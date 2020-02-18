"""Test support for optional strict bounds on transcript coordinates

https://github.com/biocommons/hgvs/issues/437

Tests to consider:
* Ensure that bounds remain checked for other sequence types

"""

import pytest


import hgvs
from hgvs.exceptions import HGVSInvalidVariantError


def test_437_RMRP_terminii(parser, am37):
    hgvs.global_config.mapping.strict_bounds = True

    # Generate n. and g. variants at terminal positions of tx, and
    # variants that are 1 base out of bounds
    s1_n = parser.parse("NR_003051.3:n.1G>T")
    e1_n = parser.parse("NR_003051.3:n.267G>T")

    s2_n = parser.parse("NR_003051.3:n.-1G>C")
    e2_n = parser.parse("NR_003051.3:n.268G>C")


    # Sanity check projections to genome for in-bound variations
    s1_g = am37.n_to_g(s1_n)
    assert str(s1_g) == "NC_000009.11:g.35658015C>A"
    e1_g = am37.n_to_g(e1_n)
    assert str(e1_g) == "NC_000009.11:g.35657749C>A"


    # Disable bounds checking and retry
    hgvs.global_config.mapping.strict_bounds = False

    # TODO: not yet handling lack of zeroes in HGVS counting
    #s2_g = am37.n_to_g(s2_n)
    #assert str(s2_g) == "NC_000009.11:g.35658016A>G"
    #assert s2_g.posedit.pos.start.base - s1_g.posedit.pos.start.base == 1

    e2_g = am37.n_to_g(e2_n)
    assert str(e2_g) == "NC_000009.11:g.35657748A>G"
    assert e2_g.posedit.pos.start.base - e1_g.posedit.pos.start.base == -1

    hgvs.global_config.mapping.strict_bounds = True


def test_437_RMRP_enforce_strict_bounds(parser, am37):
    """Ensure that an exception is raised when outside bounds"""

    # Construct two variants that are 1 base outside the tx bounds
    # N.B. G is not the correct reference for either variant
    s2_n = parser.parse("NR_003051.3:n.-1G>C")
    e2_n = parser.parse("NR_003051.3:n.268G>C")

    with pytest.raises(HGVSInvalidVariantError):
        s2_g = am37.n_to_g(s2_n)

    with pytest.raises(HGVSInvalidVariantError):
        e2_g = am37.n_to_g(e2_n)


# TODO: not yet handling lack of zeroes in HGVS counting
def x_test_437_RMRP_oob_dup(parser, am37):
    """Intentionally preserve dup, derived from genomic sequence, when
    projecting to out-of-bounds transcript coordinates"""
    hgvs.global_config.mapping.strict_bounds = False

    var_n = parser.parse("NR_003051.3:n.-19_-18insACT")
    var_g = am37.n_to_g(var_n)
    assert str(var_g) == "NC_000009.11:g.35658038_35658040dup"

    hgvs.global_config.mapping.strict_bounds = True


if __name__ == "__main__":
    from hgvs.easy import *
    test_437(parser=parser, am=am37)
