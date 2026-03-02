"""Tests for lazy module loading in package init."""

import panssr


def test_lazy_access_basic_modules():
    assert hasattr(panssr, "utils")
    assert panssr.utils is not None
    assert hasattr(panssr, "io_tools")
    assert panssr.io_tools is not None
