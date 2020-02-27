"""Test support for optional strict bounds on transcript coordinates

https://github.com/biocommons/hgvs/issues/437

RMPR (- strand, chr 9)
                        |==============================|
  NR_003051.3        -1 |    1      2   //  267    268 |  269
                      T |    G          //    G      T |    T
                        |               //             |     
 NC_000009.11  35658016 | 8015   8014   // 7749   7748 | 7748
                      A |    C          //    C      A |    A
                        |               //             |     
 NC_000009.12  35658019 | 8018   8017   // 7752   7751 | 7750
                      A |    C          //    C      A |    A
                        |==============================|


Tests to consider:
* Ensure that bounds remain checked for other sequence types

"""

import unittest

import pytest

import hgvs
from hgvs.exceptions import HGVSInvalidVariantError


@pytest.mark.usefixtures("kitchen_sink_setup")
class Test437_RMRP(unittest.TestCase):
    
    def setUp(self):
        # Generate n. and g. variants at terminal positions of tx, and
        # variants that are 1 base out of bounds

        self.s1_n = self.parser.parse("NR_003051.3:n.1G>T")
        self.e1_n = self.parser.parse("NR_003051.3:n.268T>C")
        
        self.s2_n = self.parser.parse("NR_003051.3:n.-1G>C")
        self.e2_n = self.parser.parse("NR_003051.3:n.269N>C")
        

    def test_437_terminii(self):
        hgvs.global_config.mapping.strict_bounds = True

        # Sanity check projections to genome for in-bound variations
        s1_g = self.am37.n_to_g(self.s1_n)
        assert str(s1_g) == "NC_000009.11:g.35658015C>A"
        e1_g = self.am37.n_to_g(self.e1_n)
        assert str(e1_g) == "NC_000009.11:g.35657748A>G"

        # Disable bounds checking and retry
        hgvs.global_config.mapping.strict_bounds = False

        # TODO: not yet handling lack of zeroes in HGVS counting
        s2_g = self.am37.n_to_g(self.s2_n)
        assert str(s2_g) == "NC_000009.11:g.35658016A>G"
        assert s2_g.posedit.pos.start.base - s1_g.posedit.pos.start.base == 1

        e2_g = self.am37.n_to_g(self.e2_n)
        assert str(e2_g) == "NC_000009.11:g.35657747A>G"
        assert e2_g.posedit.pos.start.base - e1_g.posedit.pos.start.base == -1

        hgvs.global_config.mapping.strict_bounds = True


    def test_437_enforce_strict_bounds(self):
        """Ensure that an exception is raised when outside bounds"""

        # TODO: These exceptions are raised, but with the wrong message
        # use `match=` arg
        with pytest.raises(HGVSInvalidVariantError):
            s2_g = self.am37.n_to_g(self.s2_n)

        with pytest.raises(HGVSInvalidVariantError):
            e2_g = self.am37.n_to_g(self.e2_n)


# TODO: not yet handling lack of zeroes in HGVS counting
def x_test_437_oob_dup(parser, am37):
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
