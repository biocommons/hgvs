#!/usr/bin/env python

import pytest

import hgvs
from hgvs.exceptions import HGVSInvalidIntervalError


def test_602(parser, am37):
    """ensure that there's no ref check for variants beyond bounds of transcript"""

    hgvs_g = "NC_000006.11:g.33415675_33422480del"
    var_g = parser.parse(hgvs_g)

    hgvs.global_config.mapping.strict_bounds = True  # default
    with pytest.raises(HGVSInvalidIntervalError):
        var_c = am37.g_to_c(var_g, "NM_006772.2")

    hgvs.global_config.mapping.strict_bounds = False
    var_c = am37.g_to_c(var_g, "NM_006772.2")  # No error


if __name__ == "__main__":
    from hgvs.easy import am37, parser

    test_602(parser, am37)
