"""
Unit tests for utils module
"""
import unittest
from panssr import utils


class TestUtils(unittest.TestCase):
    """Test utility functions"""

    def test_reverse_complement(self):
        """Test reverse complement function"""
        # Test simple sequence
        self.assertEqual(utils.reverse_complement("ATCG"), "CGAT")
        self.assertEqual(utils.reverse_complement("AAAA"), "TTTT")
        self.assertEqual(utils.reverse_complement("GCGC"), "GCGC")

        # Test with lowercase
        self.assertEqual(utils.reverse_complement("atcg"), "cgat")

        # Test palindrome
        self.assertEqual(utils.reverse_complement("GAATTC"), "GAATTC")

    def test_replace_ambiguity_codes(self):
        """Test IUPAC ambiguity code replacement"""
        # Test standard bases (no change)
        result = utils.replace_ambiguity_codes("ATCG")
        self.assertEqual(result, "ATCG")

        # Test ambiguity codes
        result = utils.replace_ambiguity_codes("R")  # R = A or G
        self.assertEqual(result, "[AG]")

        result = utils.replace_ambiguity_codes("Y")  # Y = C or T
        self.assertEqual(result, "[CT]")

        result = utils.replace_ambiguity_codes("N")  # N = any base
        self.assertEqual(result, "[ACGT]")

        # Test mixed sequence
        result = utils.replace_ambiguity_codes("ATNGC")
        self.assertEqual(result, "AT[ACGT]GC")

    def test_timeit_decorator(self):
        """Test that timeit decorator works"""
        @utils.timeit
        def dummy_function():
            return "test"

        result = dummy_function()
        self.assertEqual(result, "test")


if __name__ == '__main__':
    unittest.main()
