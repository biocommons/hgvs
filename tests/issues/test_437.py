"""https://github.com/biocommons/hgvs/issues/437"""

import pytest


import hgvs
from hgvs.exceptions import HGVSInvalidVariantError


def test_437_RMRP_terminii(parser, am37):
    hgvs.global_config.mapping.strict_bounds = True

    # Generate n. and g. variants at terminal positions of tx
    s1_n = parser.parse("NR_003051.3:n.1G>T")
    s1_g = am37.n_to_g(s1_n)
    assert str(s1_g) == "NC_000009.11:g.35658015C>A"

    e1_n = parser.parse("NR_003051.3:n.267G>T")
    e1_g = am37.n_to_g(e1_n)
    assert str(e1_g) == "NC_000009.11:g.35657749C>A"


    # Ensure that an error is raised when outside bounds
    # Construct two variants that are 1 base outside the tx bounds
    # N.B. G is not the correct reference for either variant
    s2_n = parser.parse("NR_003051.3:n.-1G>C")
    with pytest.raises(HGVSInvalidVariantError):
        s2_g = am37.n_to_g(s2_n)

    e2_n = parser.parse("NR_003051.3:n.268G>C")
    with pytest.raises(HGVSInvalidVariantError):
        e2_g = am37.n_to_g(e2_n)

    # Disable bounds checking and retry
    hgvs.global_config.mapping.strict_bounds = False

    s2_g = am37.n_to_g(s2_n)
    assert str(s2_g) == "NC_000009.11:g.35658016A>G"

    e2_g = am37.n_to_g(e2_n)
    assert str(e2_g) == "NC_000009.11:g.35657748A>G"

    # ensure that the resulting genomic variants are +/- 1 from the
    # in-bounds variants (necessarily true given the strict str()
    # checks above)
    assert s2_g.posedit.pos.start.base - s1_g.posedit.pos.start.base == 1
    assert e2_g.posedit.pos.start.base - e1_g.posedit.pos.start.base == -1

    hgvs.global_config.mapping.strict_bounds = True



def test_437_RMRP_oob_dup(parser, am37):
    """Out-of-bounds genomic dup"""
    hgvs.global_config.mapping.strict_bounds = False

    var_n = parser.parse("NR_003051.3:n.-19_-18insACT")
    var_g = am37.n_to_g(var_n)
    assert str(var_g) == "NC_000009.11:g.35658038_35658040dup"

    hgvs.global_config.mapping.strict_bounds = True


if __name__ == "__main__":
    from hgvs.easy import *
    test_437(parser=parser, am=am37)
