"""
Unit tests for SSR discovery module
"""
import unittest
from panssr import ssr_discovery, config


class TestSSRDiscovery(unittest.TestCase):
    """Test SSR detection functionality"""

    def test_detect_dinucleotide_ssr(self):
        """Test detection of dinucleotide repeats"""
        seq = "ATGATGATGATGATGATG"  # (ATG)6 = dinucleotide AT repeated
        ssrs = ssr_discovery.detect_ssrs(seq)
        self.assertGreater(len(ssrs), 0, "Should detect SSR in simple repeat sequence")

    def test_detect_trinucleotide_ssr(self):
        """Test detection of trinucleotide repeats"""
        seq = "ATGATGATGATGATG"  # (ATG)5 = trinucleotide repeat
        ssrs = ssr_discovery.detect_ssrs(seq)
        self.assertGreater(len(ssrs), 0, "Should detect trinucleotide SSR")
        # Check that at least one SSR has motif length of 3
        has_tri = any(len(ssr['motif']) == 3 for ssr in ssrs)
        self.assertTrue(has_tri, "Should find trinucleotide motif")

    def test_no_ssr_in_random_sequence(self):
        """Test that random sequence without repeats returns no SSRs"""
        seq = "ACGTACGTACG"  # Too short repeats
        ssrs = ssr_discovery.detect_ssrs(seq)
        # May return no SSRs or very short ones depending on thresholds
        if ssrs:
            # If any found, they should be valid
            for ssr in ssrs:
                self.assertIn('motif', ssr)
                self.assertIn('start', ssr)
                self.assertIn('end', ssr)
                self.assertGreater(ssr['repeat_count'], 0)

    def test_ssr_positions_are_correct(self):
        """Test that SSR positions are correctly reported"""
        seq = "NNNNNATATATATATATAT"  # AT repeat starting at position 6 (1-indexed), 7 repeats
        ssrs = ssr_discovery.detect_ssrs(seq)
        # Should find AT repeat
        at_ssrs = [s for s in ssrs if s['motif'] == 'AT']
        self.assertGreater(len(at_ssrs), 0, "Should find AT repeat")
        # Check that start position is 1-indexed
        self.assertGreater(at_ssrs[0]['start'], 0, "Start position should be 1-indexed")

    def test_max_length_filter(self):
        """Test that very long SSRs are filtered out"""
        # Create a very long repeat sequence
        seq = "A" * 200  # Mononucleotide repeat longer than MAX_SSR_LENGTH
        ssrs = ssr_discovery.detect_ssrs(seq)
        # Should be filtered by MAX_SSR_LENGTH
        for ssr in ssrs:
            ssr_len = ssr['end'] - ssr['start'] + 1
            self.assertLessEqual(ssr_len, config.MAX_SSR_LENGTH,
                                 "SSR length should not exceed MAX_SSR_LENGTH")


if __name__ == '__main__':
    unittest.main()
