import unittest

from parameterized import parameterized

import hgvs
import hgvs.dataproviders
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig
from hgvs.repeats import RepeatAnalyser, count_repetitive_units


class TestRepeats(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.hdp = hgvs.dataproviders.uta.connect()

    @parameterized.expand(
        [
            ("abc", "abc", 1),
            ("CC", "C", 2),
            ("GAGA", "GA", 2),
            ("", "", 1),
            ("AAAAAB", "AAAAAB", 1),
            ("ABBCABBCABBC", "ABBC", 3),
        ]
    )
    def test_count_repetitive_units(self, s: str, u: str, c: int):

        observed_c, observed_u = count_repetitive_units(s)

        self.assertEqual(u, observed_u)
        self.assertEqual(c, observed_c)

    @parameterized.expand(
        [
            ("NC_000019.10:g.45770205del", True, "C", 6, 5, "C[6]>C[5]"),
            ("NC_000019.10:g.45770204_45770205del", True, "C", 6, 4, "C[6]>C[4]"),
            ("NC_000019.10:g.45770205insC", True, "C", 6, 7, "C[6]>C[7]"),
            ("NC_000019.10:g.45770206del", False, None, 0, 0, "A>"),
            ("NC_000007.13:g.36561662dup", True, "C", 2, 3, "C[2]>C[3]"),
        ]
    )
    def test_homopolymer(self, hgvs_g, is_repeat, repeat_unit, ref_count, alt_count, s):

        hp = hgvs.parser.Parser()
        var_g = hp.parse_hgvs_variant(hgvs_g)
        config = PrettyConfig(self.hdp, None, None)
        dc = DataCompiler(config)
        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        self.assertEqual(s, str(ra))
        self.assertEqual(is_repeat, ra.is_repeat)
        self.assertEqual(repeat_unit, ra.repeat_unit)
        self.assertEqual(ref_count, ra.ref_count)
        self.assertEqual(alt_count, ra.alt_count)

    @parameterized.expand(
        [
            (
                "NC_000021.8:g.46020668_46020682del",
                False,
                None,
                0,
                0,
                "CTGCTGCGCCCCCAGCTGCTGCGCCCC>CTGCTGCGCCCC",
            ),
            (
                "NC_000012.11:g.33049660_33049680dup",
                False,
                None,
                0,
                0,
                "GGGGGCTGCCATGGGGCCGGTGGGGGC>GGGGGCTGCCATGGGGCCGGTGGGGGCTGCCATGGGGCCGGTGGGGGC",
            ),
            ("NC_000005.10:g.123346517_123346518insATTA", True, "ATTA", 2, 3, "ATTA[2]>ATTA[3]"),
            ("NC_000005.10:g.123346522_123346525dup", True, "ATTA", 2, 3, "ATTA[2]>ATTA[3]"),
            ("NC_000019.10:g.45770210_45770212del", True, "GCA", 20, 19, "GCA[20]>GCA[19]"),
        ]
    )
    def test_repeats(self, hgvs_g, is_repeat, repeat_unit, ref_count, alt_count, s):
        hp = hgvs.parser.Parser()
        var_g = hp.parse_hgvs_variant(hgvs_g)
        config = PrettyConfig(self.hdp, None, None)
        dc = DataCompiler(config)
        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        self.assertEqual(is_repeat, ra.is_repeat)
        self.assertEqual(repeat_unit, ra.repeat_unit)
        self.assertEqual(ref_count, ra.ref_count)
        self.assertEqual(alt_count, ra.alt_count)
        self.assertEqual(s, str(ra))
