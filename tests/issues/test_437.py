r"""Test support for optional strict bounds on transcript coordinates

https://github.com/biocommons/hgvs/issues/437

    RMPR (- strand, chr 9)

    NC_000009.12          7751    7759  //      8018           8033      8043
    NC_000009.11          7748    7756  //      8015           8030      8040
                             |       |  //         |              |         |
               +      acaaaaaACAGCCGCG  //  CACGAACCacgtcctcagcttcacagagtagtatt
               -      tgtttttTGTCGGCGC  //  GTGCTTGGtgcaggagtcgaagtgtctcatcataa
                      aaaaaaaTGTCGGCGC  //  GTGCTTGG                 
                             |       |  //         |\        |         |
    NR_003051.3            268     260  //        <1 -1    -10       -20

    tests:                  ^^                     ^^                >-------<
                            43                     21                   dup

"""

import unittest

import pytest

import hgvs
import hgvs.alignmentmapper
import hgvs.variantmapper
from hgvs.exceptions import HGVSInvalidVariantError


@pytest.mark.usefixtures("kitchen_sink_setup")
class Test437_RMRP(unittest.TestCase):
    def setUp(self):
        # Generate n. and g. variants at terminal positions of tx, and
        # variants that are 1 base out of bounds

        # "NR_003051.3:n.1G>T", "NC_000009.11:g.35658015C>A"
        # "NR_003051.3:n.268T>C", "NC_000009.11:g.35657748A>G"

        # "NR_003051.3:n.-1G>C", "NC_000009.11:g.35658016A>G"
        # "NR_003051.3:n.269N>C", "NC_000009.11:g.35657747A>G"
        # "NR_003051.3:n.-19_-18insACT", "NC_000009.11:g.35658038_35658040dup"   # ~= 40_41insCAT; NB rotated 1

        self.s1_n = self.parser.parse("NR_003051.3:n.1G>T")
        self.e1_n = self.parser.parse("NR_003051.3:n.268T>C")

        self.s2_n = self.parser.parse("NR_003051.3:n.-1G>C")
        self.e2_n = self.parser.parse("NR_003051.3:n.269N>C")

        self.vm = hgvs.variantmapper.VariantMapper(self.hdp)
        self.alm37 = hgvs.alignmentmapper.AlignmentMapper(self.hdp, "NR_003051.3", "NC_000009.11", "splign")
        self.alm38 = hgvs.alignmentmapper.AlignmentMapper(self.hdp, "NR_003051.3", "NC_000009.12", "splign")

    def test_tx_start(self):
        hgvs_n = "NR_003051.3:n.1G>T"
        hgvs_g = "NC_000009.11:g.35658015C>A"

        var_n = self.parser.parse(hgvs_n)
        var_g = self.parser.parse(hgvs_g)

        var_ng = self.am37.n_to_g(var_n)
        assert var_ng.posedit.pos == var_g.posedit.pos
        assert str(var_ng) == hgvs_g

        var_gn = self.am37.g_to_n(var_g, var_n.ac)
        assert var_gn.posedit.pos == var_n.posedit.pos
        assert str(var_gn) == hgvs_n

    def test_terminii(self):
        hgvs.global_config.mapping.strict_bounds = True

        # Sanity check projections to genome for in-bound variations
        # n -> g -> n roundtrip
        s1_g = self.am37.n_to_g(self.s1_n)
        assert str(s1_g) == "NC_000009.11:g.35658015C>A"
        assert self.s1_n == self.am37.g_to_n(s1_g, self.s1_n.ac)

        e1_g = self.am37.n_to_g(self.e1_n)
        assert str(e1_g) == "NC_000009.11:g.35657748A>G"
        assert self.e1_n == self.am37.g_to_n(e1_g, self.e1_n.ac)

        # Disable bounds checking and retry
        hgvs.global_config.mapping.strict_bounds = False

        s2_g = self.am37.n_to_g(self.s2_n)
        assert str(s2_g) == "NC_000009.11:g.35658016A>G"
        assert s2_g.posedit.pos.start.base - s1_g.posedit.pos.start.base == 1
        # BUG: assert self.s2_n == self.am37.g_to_n(s2_g, self.s2_n.ac)

        e2_g = self.am37.n_to_g(self.e2_n)
        assert str(e2_g) == "NC_000009.11:g.35657747A>G"
        assert e2_g.posedit.pos.start.base - e1_g.posedit.pos.start.base == -1
        # BUG: assert self.e2_n == self.am37.g_to_n(e2_g, self.e2_n.ac)

        hgvs.global_config.mapping.strict_bounds = True

    def test_enforce_strict_bounds(self):
        """Ensure that an exception is raised when outside bounds"""

        # TODO: These exceptions are raised, but with the wrong message
        # use `match=` arg
        with pytest.raises(HGVSInvalidVariantError):
            s2_g = self.am37.n_to_g(self.s2_n)

        with pytest.raises(HGVSInvalidVariantError):
            e2_g = self.am37.n_to_g(self.e2_n)


def test_oob_dup(parser, am37):
    """Intentionally preserve dup, derived from genomic sequence, when
    projecting to out-of-bounds transcript coordinates

    """
    hgvs.global_config.mapping.strict_bounds = False

    # n_to_g
    var_n = parser.parse("NR_003051.3:n.-19_-18insACT")
    var_g = am37.n_to_g(var_n)
    assert str(var_g) == "NC_000009.11:g.35658038_35658040dup"

    # g_to_n
    var_n2 = am37.g_to_n(var_g, "NR_003051.3")
    assert str(var_n2) == "NR_003051.3:n.-21_-19dup"

    hgvs.global_config.mapping.strict_bounds = True


def test_invitae_examples(parser, am37):
    """bidirectional g↔t tests of out-of-bounds variants provided by Invitae"""
    invitae_examples = [
        ("NC_000009.11:g.35658020C>T", "NR_003051.3:n.-5G>A"),
        # See #592: hgvs doesn't support n.* coordinates, so rewrite
        # n.*7 as n.275.
        # ("NC_000009.11:g.35657741A>G", "NR_003051.3:n.*7T>C"),
        ("NC_000009.11:g.35657741A>G", "NR_003051.3:n.275T>C"),
        # Same g. position, ensuring correct c-n equivalence
        ("NC_000001.10:g.18807339T>C", "NM_152375.2:n.-85T>C"),
        ("NC_000001.10:g.18807339T>C", "NM_152375.2:c.-137T>C"),
        # Same transcript, 5' and 3' out-of-bounds
        # disabled these tests because NCBI responses were failing.
        # Yet another reason to create a pure REST interface
        # ("NC_000007.13:g.94953895G>A", "NM_000446.5:c.-108C>T"),
        # ("NC_000007.13:g.94927495G>A", "NM_000446.5:c.*761C>T"),
    ]

    hgvs.global_config.mapping.strict_bounds = False

    for hgvs_g, hgvs_t in invitae_examples:
        var_g = parser.parse(hgvs_g)
        var_t = parser.parse(hgvs_t)
        assert hgvs_g == str(am37.t_to_g(var_t))
        if var_t.type == "c":
            assert hgvs_t == str(am37.g_to_c(var_g, var_t.ac))
        else:
            assert hgvs_t == str(am37.g_to_n(var_g, var_t.ac))

    hgvs.global_config.mapping.strict_bounds = True


if __name__ == "__main__":
    from hgvs.easy import *

    test_invitae_examples(parser=parser, am37=am37)
