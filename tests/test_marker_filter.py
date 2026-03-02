"""
Unit tests for marker filter module
"""
import unittest
from panssr import marker_filter


class TestMarkerFilter(unittest.TestCase):
    """Test marker filtering functionality"""

    def test_is_valid_ssr_polymorphism(self):
        """Test SSR polymorphism validation"""
        # Valid polymorphism: sizes differ by motif length multiples
        sizes = [200, 204, 208]  # Differ by 4bp each
        motif = "AT"  # 2bp motif
        result = marker_filter.is_valid_ssr_polymorphism(sizes, motif)
        self.assertTrue(result, "Should accept valid SSR polymorphism")

        # Invalid: not polymorphic (all same size)
        sizes = [200, 200, 200]
        result = marker_filter.is_valid_ssr_polymorphism(sizes, motif)
        self.assertFalse(result, "Should reject non-polymorphic marker")

        # Invalid: size differences not multiples of motif length
        sizes = [200, 203, 207]  # Differ by 3bp (not multiple of 2)
        motif = "AT"
        result = marker_filter.is_valid_ssr_polymorphism(sizes, motif)
        self.assertFalse(result, "Should reject invalid size differences")

    def test_has_successful_primers(self):
        """Test primer success checking"""
        # Valid primers
        marker = {
            "primers": {
                "PRIMER_LEFT_0_SEQUENCE": "ATCGATCGATCG",
                "PRIMER_RIGHT_0_SEQUENCE": "GCTAGCTAGCTA"
            }
        }
        result = marker_filter.has_successful_primers(marker)
        self.assertTrue(result, "Should accept marker with valid primers")

        # Missing primers
        marker = {"primers": {}}
        result = marker_filter.has_successful_primers(marker)
        self.assertFalse(result, "Should reject marker without primers")

        # No primers field
        marker = {}
        result = marker_filter.has_successful_primers(marker)
        self.assertFalse(result, "Should reject marker with no primer field")

    def test_filter_markers_basic(self):
        """Test basic marker filtering"""
        markers = [
            {
                "motif": "AT",
                "primers": {
                    "PRIMER_LEFT_0_SEQUENCE": "ATCG",
                    "PRIMER_RIGHT_0_SEQUENCE": "GCTA"
                },
                "amplicon_sizes": [200, 204, 208]  # Valid polymorphism
            },
            {
                "motif": "CG",
                "primers": {
                    "PRIMER_LEFT_0_SEQUENCE": "CGCG",
                    "PRIMER_RIGHT_0_SEQUENCE": "TAGC"
                },
                "amplicon_sizes": [300, 300, 300]  # Not polymorphic
            },
        ]

        filtered = marker_filter.filter_markers(markers, require_primers=True)
        self.assertEqual(len(filtered), 1, "Should filter out non-polymorphic marker")
        self.assertEqual(filtered[0]["motif"], "AT", "Should keep polymorphic marker")


class TestAnnotationPreference(unittest.TestCase):
    """Test annotation preference logic"""

    def test_is_preferred_annotation_intergenic(self):
        """Test intergenic preference"""
        # None means intergenic
        result = marker_filter.is_preferred_annotation(None, prefer_intergenic=True)
        self.assertTrue(result, "Should prefer intergenic (None)")

        # Exon annotation when prefer_intergenic=True
        annotation = {"type": "exon"}
        result = marker_filter.is_preferred_annotation(annotation, prefer_intergenic=True)
        self.assertFalse(result, "Should reject exon when prefer_intergenic=True")

    def test_is_preferred_annotation_repeat_filter(self):
        """Test that repeat regions are filtered"""
        annotation = {"type": "repeat_region"}
        result = marker_filter.is_preferred_annotation(annotation, prefer_intergenic=False)
        self.assertFalse(result, "Should reject repeat regions")


if __name__ == '__main__':
    unittest.main()
