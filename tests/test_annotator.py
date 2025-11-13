"""
Unit tests for annotator module
"""
import unittest
from panssrator import annotator


class TestAnnotator(unittest.TestCase):
    """Test annotation functionality"""

    def test_get_feature_priority(self):
        """Test feature priority scoring"""
        # High priority features
        self.assertGreater(annotator.get_feature_priority("CDS"),
                          annotator.get_feature_priority("intron"),
                          "CDS should have higher priority than intron")

        self.assertGreater(annotator.get_feature_priority("exon"),
                          annotator.get_feature_priority("gene"),
                          "Exon should have higher priority than gene")

        # Repeat regions should have low priority
        self.assertLess(annotator.get_feature_priority("repeat_region"),
                       annotator.get_feature_priority("gene"),
                       "Repeat region should have lower priority than gene")

        # Case insensitivity
        self.assertEqual(annotator.get_feature_priority("CDS"),
                        annotator.get_feature_priority("cds"),
                        "Priority should be case-insensitive")

    def test_annotate_ssr_no_chromosome(self):
        """Test annotation when chromosome not in trees"""
        ssr = {"chrom": "chr99", "start": 100, "end": 150}
        trees = {}  # Empty trees
        result = annotator.annotate_ssr(ssr, trees)
        self.assertIsNone(result, "Should return None for missing chromosome")

    def test_annotate_ssr_missing_chrom_field(self):
        """Test annotation when SSR has no chromosome field"""
        ssr = {"start": 100, "end": 150}  # No chrom field
        trees = {}
        result = annotator.annotate_ssr(ssr, trees)
        self.assertIsNone(result, "Should return None for missing chrom field")


if __name__ == '__main__':
    unittest.main()
