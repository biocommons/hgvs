# -*- coding: utf-8 -*-
import unittest

import pytest

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from hgvs.pretty.prettyprint import PrettyPrint


@pytest.mark.skip(
    reason="The pretty print tests are data hungry. If we were to add the data to the test cache, we would inflate the size of the cache. As such only running when necessary."
)
@pytest.mark.quick
@pytest.mark.models
class Test_SimplePosition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hp = hgvs.parser.Parser()
        cls.hdp = hgvs.dataproviders.uta.connect(mode=None, cache=None)
        cls.assembly_mapper37 = AssemblyMapper(cls.hdp, assembly_name="GRCh37")
        cls.assembly_mapper38 = AssemblyMapper(cls.hdp, assembly_name="GRCh38")
        cls.pp = PrettyPrint(cls.hdp, cls.assembly_mapper37, reverse_display=False)
        cls.pp.use_color = False

        cls.atta_expected_results = (
            "          :   123,346,500         123,346,520         123,346,540\n"
            + "chrom pos :   |    .    |    .    |    .    |    .    |    .\n"
            + "seq    -> : ATAAAGCTTTTCCAAATGTTATTAATTACTGGCATTGCTTTTTGCCAA\n"
            + "region    :                     |------|                    \n"
            + "tx seq <- : TATTTCGAAAAGGTTTACAATAATTAATGACCGTAACGAAAAACGGTT\n"
            + "tx pos    :  |    .    |    .   |   |    .    |    .    |   \n"
            + "          :  *20       *10      *1  2880      2870      2860\n"
            + "aa seq <- :                      TerAsnSerAlaAsnSerLysAlaLeu\n"
            + "aa pos    :                         |||            ...      \n"
            + "          :                         960                     \n"
            + "ref>alt   : ATTA[2]>ATTA[3]\n"
        )

    def test_var_c1_forward(self):
        """test c1 on -> strand"""

        hgvs_c = "NM_198689.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000021.8:g.46020522=\n"
            + "hgvs_c    : NM_198689.2:c.1=\n"
            + "hgvs_p    : NP_941962.1:p.Met1?\n"
            + "          :         46,020,510          46,020,530\n"
            + "chrom pos :    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : CCTCCAGTTCAATCCCCAGCATGGCCGCGTCCACTATGTCT\n"
            + "tx ref dif:                          X               \n"
            + "region    :                     =                    \n"
            + "tx seq -> : cctccagttcaatccccagcATGGCTGCGTCCACTATGTCT\n"
            + "tx pos    : |    .    |    .    |   .    |    .    | \n"
            + "          : -20       -10       1        10        20\n"
            + "aa seq -> :                     MetAlaAlaSerThrMetSer\n"
            + "aa pos    :                                 ...      \n"
            + "          :                     1                    \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_c1_reverse(self):
        """test c1 on <- strand"""

        hgvs_c = "NM_001177507.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36763753=\n"
            + "hgvs_c    : NM_001177507.2:c.1=\n"
            + "hgvs_p    : NP_001170978.1:p.Met1?\n"
            + "          :        36,763,740          36,763,760\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : GATTTTCCAGGGGGACTGCATCTCCGAGCTATGCACCCCAA\n"
            + "region    :                     =                    \n"
            + "tx seq <- : CTAAAAGGTCCCCCTGACGTAgaggctcgatacgtggggtt\n"
            + "tx pos    :  |    .    |    .   |    .    |    .    |\n"
            + "          :  20        10       1         -10       -20\n"
            + "aa seq <- : IleLysTrpProSerGlnMet                    \n"
            + "aa pos    :       ...                                \n"
            + "          :                   1                      \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_c1_reverse_flipped_display(self):
        """test the reversed display on <- strand"""
        hgvs_c = "NM_001177507.2:c.1="
        var_c = self.hp.parse(hgvs_c)

        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            show_reverse_strand=True,
            reverse_display=True,
        )
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36763753=\n"
            + "hgvs_c    : NM_001177507.2:c.1=\n"
            + "hgvs_p    : NP_001170978.1:p.Met1?\n"
            + "          :    36,763,770          36,763,750\n"
            + "chrom pos :    |    .    |    .    |    .    |    .  \n"
            + "seq    <- : AACCCCACGTATCGAGCCTCTACGTCAGGGGGACCTTTTAG\n"
            + "seq    -> : TTGGGGTGCATAGCTCGGAGATGCAGTCCCCCTGGAAAATC\n"
            + "region    :                     =                    \n"
            + "tx seq -> : ttggggtgcatagctcggagATGCAGTCCCCCTGGAAAATC\n"
            + "tx pos    : |    .    |    .    |   .    |    .    | \n"
            + "          : -20       -10       1        10        20\n"
            + "aa seq -> :                     MetGlnSerProTrpLysIle\n"
            + "aa pos    :                                 ...      \n"
            + "          :                     1                    \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_substitution(self):
        hgvs_g = "NC_000007.13:g.36561662C>T"
        var_g = self.hp.parse(hgvs_g)

        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper38,
            show_reverse_strand=True,
            reverse_display=False,
            use_color=False,
        )

        result = pp.display(var_g, "NM_001177507.2")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561662C>T\n"
            + "hgvs_c    : NM_001177507.2:c.1486G>A\n"
            + "hgvs_p    : NP_001170978.1:p.(Gly496Arg)\n"
            + "          :         36,561,650          36,561,670\n"
            + "chrom pos :    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "seq    <- : ATGGAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "region    :                     T                    \n"
            + "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :       |    .    |    .    |    .    |    \n"
            + "          :       1500      1490      1480      1470\n"
            + "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "aa pos    :       |||            ...            |||  \n"
            + "          :       500                           490  \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_ins(self):
        """[ATTA]x2 -> x3"""
        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001166226.1")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA\n"
            + "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT\n"
            + "hgvs_p    : NP_001159698.1:p.?\n"
            + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_atta_forward(self):
        """the ATTA[2]>ATTA[3] variant now displayed forward facing:"""

        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)
        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            show_reverse_strand=True,
            reverse_display=True,
        )
        result = pp.display(var_g, "NM_001166226.1")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA\n"
            + "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT\n"
            + "hgvs_p    : NP_001159698.1:p.?\n"
            + "          :      123,346,540         123,346,520         123,346,500\n"
            + "chrom pos : .    |    .    |    .    |    .    |    .    |  \n"
            + "seq    <- : AACCGTTTTTCGTTACGGTCATTAATTATTGTAAACCTTTTCGAAATA\n"
            + "seq    -> : TTGGCAAAAAGCAATGCCAGTAATTAATAACATTTGGAAAAGCTTTAT\n"
            + "region    :                     |------|                    \n"
            + "tx seq -> : TTGGCAAAAAGCAATGCCAGTAATTAATAACATTTGGAAAAGCTTTAT\n"
            + "tx pos    :    |    .    |    .    |   |   .    |    .    | \n"
            + "          :    2860      2870      2880         *10       *20\n"
            + "aa seq -> : LeuAlaLysSerAsnAlaSerAsnTer                     \n"
            + "aa pos    :       ...            |||                        \n"
            + "          :                      960                        \n"
            + "ref>alt   : ATTA[2]>ATTA[3]\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_var_g_dup(self):
        hgvs_g = "NC_000005.10:g.123346522_123346525dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001166226.1")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346522_123346525dup\n"
            + "hgvs_c    : NM_001166226.1:c.2880_2883dup\n"
            + "hgvs_p    : NP_001159698.1:p.(=)\n"
            + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_insertion(self):
        "A shuffleable insertion, shuffleable unit: TCGTCATC additional residues: G"
        hgvs_g = "NC_000004.11:g.1643284_1643285insTCGTCATCG"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000004.11:g.1643284_1643285insTCGTCATCG\n"
            + "hgvs_c    : NM_001174070.2:c.932_933insCGATGACGA\n"
            + "hgvs_p    : NP_001167541.1:p.(Asp309_Asp311dup)\n"
            + "          :      1,643,270 1,643,280 1,643,290 1,643,300 1,643,310\n"
            + "chrom pos : .    |    .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TCACTGGGGTGTCATCCTCATCGTCATCTTCGTAATTGAGGGAGCAAA\n"
            + "region    :                     |------|                    \n"
            + "tx seq <- : AGTGACCCCACAGTAGGAGTAGCAGTAGAAGCATTAACTCCCTCGTTT\n"
            + "tx pos    :   |    .    |    .    |    .    |    .    |    .\n"
            + "          :   950       940       930       920       910\n"
            + "aa seq <- : sValProThrAspAspGluAspAspAspGluTyrAsnLeuSerCysLe\n"
            + "aa pos    :        ...            |||            ...        \n"
            + "          :                       310                       \n"
            + "ref>alt   : TCGTCATC>TCGTCATCGTCGTCATC\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_insertion_size_1(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663insT"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561662_36561663insT\n"
            + "hgvs_c    : NM_001177507.2:c.1485_1486insA\n"
            + "hgvs_p    : NP_001170978.1:p.(Gly496ArgfsTer39)\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |  \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                    ^^                   \n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :      |    .    |    .    |    .    |    \n"
            + "          :      1500      1490      1480      1470\n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "aa pos    :      |||            ...            |||  \n"
            + "          :      500                           490  \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_2bp(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561662_36561663del\n"
            + "hgvs_c    : NM_001177507.2:c.1485_1486del\n"
            + "hgvs_p    : NP_001170978.1:p.(Asp495GlufsTer39)\n"
            + "          :         36,561,650          36,561,670\n"
            + "chrom pos :    .    |    .    |    .    |    .    |   \n"
            + "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     xx                    \n"
            + "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :       |    .    |    .    |    .    |    .\n"
            + "          :       1500      1490      1480      1470\n"
            + "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "aa pos    :       |||            ...            |||   \n"
            + "          :       500                           490   \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_first_base(self):
        hgvs_n = "NM_198689.2:n.1C>G"
        var_n = self.hp.parse(hgvs_n)

        result = self.pp.display(var_n)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000021.8:g.46020497C>G\n"
            + "hgvs_n    : NM_198689.2:n.1C>G\n"
            + "          :    46,020,480          46,020,500\n"
            + "chrom pos :    |    .    |    .    |    .    |    .  \n"
            + "seq    -> : CTCACTCACCCACTCACTCCCATCTCCTCCAGTTCAATCCC\n"
            + "region    :                     G                    \n"
            + "tx seq -> :                     CATCTCCTCCAGTTCAATCCC\n"
            + "tx pos    :                         .    |    .    | \n"
            + "aa seq -> :                                          \n"
            + "aa pos    :                                          \n"
            + "          :                                          \n"
        ).split("\n")

        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_intergenic(self):
        hgvs_g = "NC_000021.8:g.29894C>A"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000021.8:g.29894C>A\n"
            + "          :       29,880    29,890    29,900    29,910\n"
            + "chrom pos :  .    |    .    |    .    |    .    |    \n"
            + "seq    -> : NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            + "region    :                     A                    \n"
            + "tx pos    :                                          \n"
            + "aa pos    :                                          \n"
            + "          :                                          "
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561662del\n"
            + "hgvs_c    : NM_001177507.2:c.1487del\n"
            + "hgvs_p    : NP_001170978.1:p.(Gly496AspfsTer122)\n"
            + "          :          36,561,650          36,561,670\n"
            + "chrom pos :     .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                     xx                    \n"
            + "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :        |    .    |    .    |    .    |    \n"
            + "          :        1500      1490      1480      1470\n"
            + "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "aa pos    :        |||            ...            |||  \n"
            + "          :        500                           490  \n"
            + "ref>alt   : C[2]>C[1]\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_del_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561663del\n"
            + "hgvs_c    : NM_001177507.2:c.1485del\n"
            + "hgvs_p    : NP_001170978.1:p.(Asp495GlufsTer123)\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     x                    \n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "aa pos    :      |||            ...            |||   \n"
            + "          :      500                           490   \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_dup_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561662dup\n"
            + "hgvs_c    : NM_001177507.2:c.1487dup\n"
            + "hgvs_p    : NP_001170978.1:p.(Phe497IlefsTer38)\n"
            + "          :          36,561,650          36,561,670\n"
            + "chrom pos :     .    |    .    |    .    |    .    |  \n"
            + "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT\n"
            + "region    :                     ||                    \n"
            + "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA\n"
            + "tx pos    :        |    .    |    .    |    .    |    \n"
            + "          :        1500      1490      1480      1470\n"
            + "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl\n"
            + "aa pos    :        |||            ...            |||  \n"
            + "          :        500                           490  \n"
            + "ref>alt   : C[2]>C[3]\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_dup_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561663dup\n"
            + "hgvs_c    : NM_001177507.2:c.1485dup\n"
            + "hgvs_p    : NP_001170978.1:p.(Gly496TrpfsTer39)\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     |                    \n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "aa pos    :      |||            ...            |||   \n"
            + "          :      500                           490   \n"
            + "ref>alt   : A>AA\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_identity(self):
        hgvs_g = "NC_000007.13:g.36561663="
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000007.13:g.36561663=\n"
            + "hgvs_c    : NM_001177507.2:c.1485=\n"
            + "hgvs_p    : NP_001170978.1:p.(Asp495=)\n"
            + "          :        36,561,650          36,561,670\n"
            + "chrom pos :   .    |    .    |    .    |    .    |   \n"
            + "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG\n"
            + "region    :                     =                    \n"
            + "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC\n"
            + "tx pos    :      |    .    |    .    |    .    |    .\n"
            + "          :      1500      1490      1480      1470\n"
            + "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln\n"
            + "aa pos    :      |||            ...            |||   \n"
            + "          :      500                           490   \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_tiny(self):
        """Test a variant with bad input."""
        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)

        tiny_pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            padding_left=0,
            padding_right=0,
            reverse_display=False,
        )

        result = tiny_pp.display(var_g, "NM_001166226.1")
        print(result)

        result = result.split("\n")

        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA\n"
            + "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT\n"
            + "hgvs_p    : NP_001159698.1:p.?\n"
            + "          :   123,346,520\n"
            + "chrom pos :   |    .\n"
            + "seq    -> : ATTAATTA\n"
            + "region    : |------|\n"
            + "tx seq <- : TAATTAAT\n"
            + "tx pos    : |   |   \n"
            + "          : *1  2880\n"
            + "aa seq <- :  TerAsnS\n"
            + "aa pos    :     ||| \n"
            + "          :     960 \n"
            + "ref>alt   : ATTA[2]>ATTA[3]"
        ).split("\n")

        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    @pytest.mark.skip(reason="CNVs not implemented yet")
    def test_cnv(self):
        """Test a CNV variant. TODO: make display compact"""
        hgvs_g = "NC_000005.10:g.123345517_123346518del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)

        print(result)

    def test_hgvs_c(self):
        """Test a hgvs_c variant overlapping start codon on reverse strand."""
        hgvs_c = "NM_004572.3:c.-9_12dup"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            padding_left=10,
            padding_right=110,
            use_color=False,
            show_legend=False,
            reverse_display=False,
        )
        result = pp.display(var_c)

        print(result)

        result = result.split("\n")
        expected_str = (
            "NC_000012.11:g.33049660_33049680dup\n"
            + "NM_004572.3:c.-9_12dup\n"
            + "NP_004563.2:p.(Met1_Pro4dup)\n"
            + "      33,049,650          33,049,670          33,049,690          33,049,710          33,049,730          33,049,750          33,049,770          33,049,790\n"
            + " .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |\n"
            + "CTGGGGCGCCGGGGGCTGCCATGGGGCCGGTGGGGGCGACCGAGCTGCTCGCCTGCCTCTGGACTCGCGGGCGAAGCCGCCACGGAGCTGGGGGCGCTGGCGCGAGCCCCGCCCCGCTCGAGTCCGGCCCCGCCCCTGGCCCGCCCC\n"
            + "          |-------------------------|                                                                                                              \n"
            + "GACCCCGCGGCCCCCGACGGTAccccggccacccccgctggctcgacgagcggacggagacctgagcgcccgcttcggcggtgcctcgacccccgcgaccgcgctcggggcggggcgagctcaggccggggcgggga          \n"
            + "  |    .    |    .   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .          \n"
            + "  20        10       1         -10       -20       -30       -40       -50       -60       -70       -80       -90       -100      -110\n"
            + "aProAlaGlyProAlaAlaMet                                                                                                                             \n"
            + "       ...                                                                                                                                         \n"
            + "                   1                                                                                                                               \n"
            + "GGGGGCTGCCATGGGGCCGGTGGGGGC>GGGGGCTGCCATGGGGCCGGTGGGGGCTGCCATGGGGCCGGTGGGGGC"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree(self):
        """Test a tx ref disagree variant."""
        hgvs_g = "NM_001111.4:c.298G>A"
        var_g = self.hp.parse(hgvs_g)
        result = self.pp.display(var_g)

        print(result)

        result = result.split("\n")

        # note the X in the transscript sequence
        expected_str = (
            "hgvs_g    : NC_000001.10:g.154574820=\n"
            + "hgvs_c    : NM_001111.4:c.298G>A\n"
            + "hgvs_p    : NP_001102.2:p.(Gly100Arg)\n"
            + "          : 154,574,800         154,574,820         154,574,840\n"
            + "chrom pos : |    .    |    .    |    .    |    .    |\n"
            + "seq    -> : TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA\n"
            + "tx ref dif:                     X                    \n"
            + "region    :                     =                    \n"
            + "tx seq <- : AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT\n"
            + "tx pos    :    .    |    .    |    .    |    .    |  \n"
            + "          :         310       300       290       280\n"
            + "aa seq <- : ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa\n"
            + "aa pos    :    ...            |||            ...     \n"
            + "          :                   100                    \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_ref_ins(self):
        """A ref disagree with a region inserted inref, that is missing in transcript"""
        # hgvs_g = "NC_000001.10:g.154574820_154574821delinsCA"
        # var_g = self.hp.parse(hgvs_g)
        # hgvs_c = "NM_020469.2:c.188_189="
        # hgvs_c = "NM_003777.3:c.5475dup" # an I variant
        # hgvs_c =

        hgvs_c = "NM_198689.2:c.124_135="
        # this would match chromosome: "NM_198689.2:c.124_135insCTGCTGCGCCCCCAG"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            infer_hgvs_c=True,
            padding_left=30,
            padding_right=40,
        )
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000021.8:g.46020668_46020682del\n"
            "hgvs_c    : NM_198689.2:c.124_135=\n"
            "hgvs_p    : NP_941962.1:p.(Cys42=)\n"
            "          :     46,020,630          46,020,650          46,020,670          46,020,690          46,020,710\n"
            "chrom pos :     |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  \n"
            "seq    -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCCCTGCTGCGCCCCCAGCTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCCGT\n"
            "tx ref dif:                               IIIIIIIIIIIIIII                                                 XX \n"
            "region    :                               x-------------------------x                                        \n"
            "tx seq -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCC---------------CTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCTAT\n"
            "tx pos    : .    |    .    |    .    |                   .    |    .    |    .    |    .    |    .    |    . \n"
            "          :      110       120       130                      140       150       160       170       180\n"
            "aa seq -> : pAspCysProGluSerCysCysGluProPr---------------oCysCysAlaProAlaProCysLeuSerLeuValCysThrProValSerTyr\n"
            "aa pos    : .            |||            ..               .            |||            ...            |||      \n"
            "          :              40                                           50                            60       \n"
            "ref>alt   : CTGCTGCGCCCCCAGCTGCTGCGCCCC>CTGCTGCGCCCC\n"
        ).split("\n")

        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_del(self):
        # hgvs_g = "NC_000001.10:g.154574820_154574821delinsCA" - one base is a svn relative to the tx and part of the variant -> NM_001025107.2:c.-589C>T
        # var_g = self.hp.parse(hgvs_g)
        # hgvs_c = "NM_020469.2:c.188_189=" is  NC_000009.11:g.136135237_136135238delinsGC in ref

        hgvs_c = "NM_000682.6:c.901_911del"  # a del variant
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            infer_hgvs_c=True,
            padding_left=30,
            padding_right=40,
            reverse_display=False,
        )
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")

        expected_str = (
            "hgvs_g    : NC_000002.11:g.96780987_96780997del\n"
            + "hgvs_c    : NM_000682.6:c.901_911del\n"
            + "hgvs_p    : NP_000673.2:p.(Glu301GlyfsTer7)\n"
            + "          :     96,780,960          96,780,980                   96,781,000          96,781,020\n"
            + "chrom pos :     |    .    |    .    |    .    |    .  _________  |    .    |    .    |    .    |    .  \n"
            + "seq    -> : TGCCTGGGGTTCACACTCTTCCTCCTCCTCCTCCTCCTCTTC.........AGCTTCATCCTCTGGAGATGCCCCACAAACACCCTCCTTC\n"
            + "tx ref dif:                                           DDDDDDDDD                                        \n"
            + "region    :                               x----------x                                                 \n"
            + "tx seq <- : ACGGACCCCAAGTGTGAGAAGGAGGAGGAGGAGGAGGAGAAGGAGGAGAAGTCGAAGTAGGAGACCTCTACGGGGTGTTTGTGGGAGGAAG\n"
            + "tx pos    :   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .   \n"
            + "          :   940       930       920       910       900       890       880       870       860\n"
            + "aa seq <- : AlaGlnProGluCysGluGluGluGluGluGluGluGluGluGluGluGluAlaGluAspGluProSerAlaGlyCysValGlyGluLysG\n"
            + "aa pos    :             |||            ...            |||            ...            |||            ... \n"
            + "          :             310                           300                           290                \n"
            + "ref>alt   : CTCCTCCTCTTC>C\n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_ref_disagree_del_reverse(self):
        hgvs_c = "NM_000682.6:c.901_911del"  # a del variant

        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp,
            self.assembly_mapper37,
            infer_hgvs_c=True,
            padding_left=30,
            padding_right=40,
            reverse_display=True,
        )
        result = pp.display(var_c)
        print(result)
        result = result.split("\n")

        expected_str = (
            "hgvs_g    : NC_000002.11:g.96780987_96780997del\n"
            "hgvs_c    : NM_000682.6:c.901_911del\n"
            "hgvs_p    : NP_000673.2:p.(Glu301GlyfsTer7)\n"
            "          :        96,781,030          96,781,010                   96,780,990          96,780,970\n"
            "chrom pos :   .    |    .    |    .    |    .    |  _________  .    |    .    |    .    |    .    |    \n"
            "seq    <- : CTTCCTCCCACAAACACCCCGTAGAGGTCTCCTACTTCGA.........CTTCTCCTCCTCCTCCTCCTCCTTCTCACACTTGGGGTCCGT\n"
            "seq    -> : GAAGGAGGGTGTTTGTGGGGCATCTCCAGAGGATGAAGCT.........GAAGAGGAGGAGGAGGAGGAGGAAGAGTGTGAACCCCAGGCA\n"
            "tx ref dif:                                         DDDDDDDDD                                          \n"
            "region    :                                                  x----------x                              \n"
            "tx seq -> : GAAGGAGGGTGTTTGTGGGGCATCTCCAGAGGATGAAGCTGAAGAGGAGGAAGAGGAGGAGGAGGAGGAGGAAGAGTGTGAACCCCAGGCA\n"
            "tx pos    :    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  \n"
            "          :         860       870       880       890       900       910       920       930       940\n"
            "aa seq -> : nLysGluGlyValCysGlyAlaSerProGluAspGluAlaGluGluGluGluGluGluGluGluGluGluGluGluCysGluProGlnAla\n"
            "aa pos    :  ...            |||            ...            |||            ...            |||            \n"
            "          :                 290                           300                           310            \n"
            "ref>alt   : CTCCTCCTCTTC>C"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    @pytest.mark.skip(
        reason="actually not that special, but still a nice variant since there is a large shuffle-able sequence on both ends."
    )
    def test_exon_boundary_overlap_forward_strand(self):
        hgvs_c = "NM_001283009.2:c.1228_1266+39del"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp, self.assembly_mapper37, show_legend=True, use_color=False
        )

        result = pp.display(var_c)

        print(result)

    def test_ruler(self):
        """Test the ruler display option turned on."""
        hgvs_c = "NM_001111.4:c.298G>A"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(
            self.hdp, self.assembly_mapper37, show_legend=False, reverse_display=False
        )

        result = pp.display(var_c, "NM_001111.4")

        print(result)

        result = result.split("\n")

        # note the X in the transscript ref-disagree row
        expected_str = (
            "NC_000001.10:g.154574820=\n"
            + "NM_001111.4:c.298G>A\n"
            + "NP_001102.2:p.(Gly100Arg)\n"
            + "154,574,800         154,574,820         154,574,840\n"
            + "|    .    |    .    |    .    |    .    |\n"
            + "TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA\n"
            + "                    X                    \n"
            + "                    =                    \n"
            + "AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT\n"
            + "   .    |    .    |    .    |    .    |  \n"
            + "        310       300       290       280\n"
            + "ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa\n"
            + "   ...            |||            ...     \n"
            + "                  100                    \n"
        ).split("\n")
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)

    def test_rna_coding(self):
        """a rna coding transcript."""
        hgvs_n = "NR_146230.2:n.10G>A"
        var_n = self.hp.parse(hgvs_n)

        pp = PrettyPrint(self.hdp, self.assembly_mapper37, show_reverse_strand=True)
        result = pp.display(var_n)
        print(result)
        expected_str = (
            "hgvs_g    : NC_000001.10:g.167905930G>A\n"
            + "hgvs_n    : NR_146230.2:n.10G>A\n"
            + "          : 167,905,910         167,905,930         167,905,950\n"
            + "chrom pos : |    .    |    .    |    .    |    .    |\n"
            + "seq    -> : TGCTGATCTTTGGATGTTCTGGTTAGTCTAAGAAGGAGAGT\n"
            + "seq    <- : ACGACTAGAAACCTACAAGACCAATCAGATTCTTCCTCTCA\n"
            + "region    :                     A                    \n"
            + "tx seq -> :            GGATGTTCTGGTTAGTCTAAGAAGGAGAGT\n"
            + "tx pos    :                .    |    .    |    .    |\n"
            + "          :                     10        20        30\n"
            + "aa seq -> :                                          \n"
            + "aa pos    :                                          \n"
            + "          :                                          \n"
        )
        for r, e in zip(result, expected_str):
            self.assertEqual(e, r)
