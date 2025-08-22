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
            "chrom pos :   |    .    |    .    |    .    |    .    |    .\n"
            "seq    -> : ATAAAGCTTTTCCAAATGTTATTAATTACTGGCATTGCTTTTTGCCAA\n"
            "region    :                     |------|                    \n"
            "tx seq <- : TATTTCGAAAAGGTTTACAATAATTAATGACCGTAACGAAAAACGGTT\n"
            "tx pos    :  |    .    |    .   |   |    .    |    .    |   \n"
            "          :  *20       *10      *1  2880      2870      2860\n"
            "aa seq <- :                      TerAsnSerAlaAsnSerLysAlaLeu\n"
            "aa pos    :                         |||            ...      \n"
            "          :                         960                     \n"
            "ref>alt   : ATTA[2]>ATTA[3]\n"
        )

    def test_var_c1_forward(self):
        """test c1 on -> strand"""
        hgvs_c = "NM_198689.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000021.8:g.46020522=",
            "hgvs_c    : NM_198689.2:c.1=",
            "hgvs_p    : NP_941962.1:p.Met1?",
            "          :         46,020,510          46,020,530",
            "chrom pos :    .    |    .    |    .    |    .    |  ",
            "seq    -> : CCTCCAGTTCAATCCCCAGCATGGCCGCGTCCACTATGTCT",
            "tx ref dif:                          X               ",
            "region    :                     =                    ",
            "tx seq -> : cctccagttcaatccccagcATGGCTGCGTCCACTATGTCT",
            "tx pos    : |    .    |    .    |   .    |    .    | ",
            "          : -20       -10       1        10        20",
            "aa seq -> :                     MetAlaAlaSerThrMetSer",
            "aa pos    :                                 ...      ",
            "          :                     1                    ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_var_c1_reverse(self):
        """test c1 on <- strand"""
        hgvs_c = "NM_001177507.2:c.1="
        var_c = self.hp.parse(hgvs_c)
        result = self.pp.display(var_c)
        print(result)
        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36763753=",
            "hgvs_c    : NM_001177507.2:c.1=",
            "hgvs_p    : NP_001170978.1:p.Met1?",
            "          :        36,763,740          36,763,760",
            "chrom pos :   .    |    .    |    .    |    .    |   ",
            "seq    -> : GATTTTCCAGGGGGACTGCATCTCCGAGCTATGCACCCCAA",
            "region    :                     =                    ",
            "tx seq <- : CTAAAAGGTCCCCCTGACGTAgaggctcgatacgtggggtt",
            "tx pos    :  |    .    |    .   |    .    |    .    |",
            "          :  20        10       1         -10       -20",
            "aa seq <- : IleLysTrpProSerGlnMet                    ",
            "aa pos    :       ...                                ",
            "          :                   1                      ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36763753=",
            "hgvs_c    : NM_001177507.2:c.1=",
            "hgvs_p    : NP_001170978.1:p.Met1?",
            "          :    36,763,770          36,763,750",
            "chrom pos :    |    .    |    .    |    .    |    .  ",
            "seq    <- : AACCCCACGTATCGAGCCTCTACGTCAGGGGGACCTTTTAG",
            "seq    -> : TTGGGGTGCATAGCTCGGAGATGCAGTCCCCCTGGAAAATC",
            "region    :                     =                    ",
            "tx seq -> : ttggggtgcatagctcggagATGCAGTCCCCCTGGAAAATC",
            "tx pos    : |    .    |    .    |   .    |    .    | ",
            "          : -20       -10       1        10        20",
            "aa seq -> :                     MetGlnSerProTrpLysIle",
            "aa pos    :                                 ...      ",
            "          :                     1                    ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561662C>T",
            "hgvs_c    : NM_001177507.2:c.1486G>A",
            "hgvs_p    : NP_001170978.1:p.(Gly496Arg)",
            "          :         36,561,650          36,561,670",
            "chrom pos :    .    |    .    |    .    |    .    |  ",
            "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT",
            "seq    <- : ATGGAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA",
            "region    :                     T                    ",
            "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA",
            "tx pos    :       |    .    |    .    |    .    |    ",
            "          :       1500      1490      1480      1470",
            "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGl",
            "aa pos    :       |||            ...            |||  ",
            "          :       500                           490  ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_var_g_ins(self):
        """[ATTA]x2 -> x3"""
        hgvs_g = "NC_000005.10:g.123346517_123346518insATTA"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001166226.1")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA\n"
            "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT\n"
            "hgvs_p    : NP_001159698.1:p.?\n" + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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
        expected_str = [
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA",
            "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT",
            "hgvs_p    : NP_001159698.1:p.?",
            "          :      123,346,540         123,346,520         123,346,500",
            "chrom pos : .    |    .    |    .    |    .    |    .    |  ",
            "seq    <- : AACCGTTTTTCGTTACGGTCATTAATTATTGTAAACCTTTTCGAAATA",
            "seq    -> : TTGGCAAAAAGCAATGCCAGTAATTAATAACATTTGGAAAAGCTTTAT",
            "region    :                     |------|                    ",
            "tx seq -> : TTGGCAAAAAGCAATGCCAGTAATTAATAACATTTGGAAAAGCTTTAT",
            "tx pos    :    |    .    |    .    |   |   .    |    .    | ",
            "          :    2860      2870      2880         *10       *20",
            "aa seq -> : LeuAlaLysSerAsnAlaSerAsnTer                     ",
            "aa pos    :       ...            |||                        ",
            "          :                      960                        ",
            "ref>alt   : ATTA[2]>ATTA[3]",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_var_g_dup(self):
        hgvs_g = "NC_000005.10:g.123346522_123346525dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001166226.1")
        print(result)
        result = result.split("\n")
        expected_str = (
            "hgvs_g    : NC_000005.10:g.123346522_123346525dup\n"
            "hgvs_c    : NM_001166226.1:c.2880_2883dup\n"
            "hgvs_p    : NP_001159698.1:p.(=)\n" + self.atta_expected_results
        ).split("\n")
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_insertion(self):
        "A shuffleable insertion, shuffleable unit: TCGTCATC additional residues: G"
        hgvs_g = "NC_000004.11:g.1643284_1643285insTCGTCATCG"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000004.11:g.1643284_1643285insTCGTCATCG",
            "hgvs_c    : NM_001174070.2:c.932_933insCGATGACGA",
            "hgvs_p    : NP_001167541.1:p.(Asp309_Asp311dup)",
            "          :      1,643,270 1,643,280 1,643,290 1,643,300 1,643,310",
            "chrom pos : .    |    .    |    .    |    .    |    .    |  ",
            "seq    -> : TCACTGGGGTGTCATCCTCATCGTCATCTTCGTAATTGAGGGAGCAAA",
            "region    :                     |------|                    ",
            "tx seq <- : AGTGACCCCACAGTAGGAGTAGCAGTAGAAGCATTAACTCCCTCGTTT",
            "tx pos    :   |    .    |    .    |    .    |    .    |    .",
            "          :   950       940       930       920       910",
            "aa seq <- : sValProThrAspAspGluAspAspAspGluTyrAsnLeuSerCysLe",
            "aa pos    :        ...            |||            ...        ",
            "          :                       310                       ",
            "ref>alt   : TCGTCATC>TCGTCATCGTCGTCATC",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_insertion_size_1(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663insT"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561662_36561663insT",
            "hgvs_c    : NM_001177507.2:c.1485_1486insA",
            "hgvs_p    : NP_001170978.1:p.(Gly496ArgfsTer39)",
            "          :        36,561,650          36,561,670",
            "chrom pos :   .    |    .    |    .    |    .    |  ",
            "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT",
            "region    :                    ^^                   ",
            "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA",
            "tx pos    :      |    .    |    .    |    .    |    ",
            "          :      1500      1490      1480      1470",
            "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGl",
            "aa pos    :      |||            ...            |||  ",
            "          :      500                           490  ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_del_2bp(self):
        hgvs_g = "NC_000007.13:g.36561662_36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561662_36561663del",
            "hgvs_c    : NM_001177507.2:c.1485_1486del",
            "hgvs_p    : NP_001170978.1:p.(Asp495GlufsTer39)",
            "          :         36,561,650          36,561,670",
            "chrom pos :    .    |    .    |    .    |    .    |   ",
            "seq    -> : TACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG",
            "region    :                     xx                    ",
            "tx seq <- :    GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC",
            "tx pos    :       |    .    |    .    |    .    |    .",
            "          :       1500      1490      1480      1470",
            "aa seq <- :    GluAsnProHisPheGlyAspValProGluIleLeuGln",
            "aa pos    :       |||            ...            |||   ",
            "          :       500                           490   ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_first_base(self):
        hgvs_n = "NM_198689.2:n.1C>G"
        var_n = self.hp.parse(hgvs_n)

        result = self.pp.display(var_n)
        print(result)
        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000021.8:g.46020497C>G",
            "hgvs_n    : NM_198689.2:n.1C>G",
            "          :    46,020,480          46,020,500",
            "chrom pos :    |    .    |    .    |    .    |    .  ",
            "seq    -> : CTCACTCACCCACTCACTCCCATCTCCTCCAGTTCAATCCC",
            "region    :                     G                    ",
            "tx seq -> :                     CATCTCCTCCAGTTCAATCCC",
            "tx pos    :                         .    |    .    | ",
            "aa seq -> :                                          ",
            "aa pos    :                                          ",
            "          :                                          ",
            "",
        ]

        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_intergenic(self):
        hgvs_g = "NC_000021.8:g.29894C>A"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g)
        print(result)
        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000021.8:g.29894C>A",
            "          :       29,880    29,890    29,900    29,910",
            "chrom pos :  .    |    .    |    .    |    .    |    ",
            "seq    -> : NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
            "region    :                     A                    ",
            "tx pos    :                                          ",
            "aa pos    :                                          ",
            "          :                                          ",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_del_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561662del",
            "hgvs_c    : NM_001177507.2:c.1487del",
            "hgvs_p    : NP_001170978.1:p.(Gly496AspfsTer122)",
            "          :          36,561,650          36,561,670",
            "chrom pos :     .    |    .    |    .    |    .    |  ",
            "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT",
            "region    :                     xx                    ",
            "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA",
            "tx pos    :        |    .    |    .    |    .    |    ",
            "          :        1500      1490      1480      1470",
            "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl",
            "aa pos    :        |||            ...            |||  ",
            "          :        500                           490  ",
            "ref>alt   : C[2]>C[1]",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_del_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663del"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561663del",
            "hgvs_c    : NM_001177507.2:c.1485del",
            "hgvs_p    : NP_001170978.1:p.(Asp495GlufsTer123)",
            "          :        36,561,650          36,561,670",
            "chrom pos :   .    |    .    |    .    |    .    |   ",
            "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG",
            "region    :                     x                    ",
            "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC",
            "tx pos    :      |    .    |    .    |    .    |    .",
            "          :      1500      1490      1480      1470",
            "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln",
            "aa pos    :      |||            ...            |||   ",
            "          :      500                           490   ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_dup_1bp_shuffleable(self):
        hgvs_g = "NC_000007.13:g.36561662dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561662dup",
            "hgvs_c    : NM_001177507.2:c.1487dup",
            "hgvs_p    : NP_001170978.1:p.(Phe497IlefsTer38)",
            "          :          36,561,650          36,561,670",
            "chrom pos :     .    |    .    |    .    |    .    |  ",
            "seq    -> : TTACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCT",
            "region    :                     ||                    ",
            "tx seq <- :     GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGA",
            "tx pos    :        |    .    |    .    |    .    |    ",
            "          :        1500      1490      1480      1470",
            "aa seq <- :     GluAsnProHisPheGlyAspValProGluIleLeuGl",
            "aa pos    :        |||            ...            |||  ",
            "          :        500                           490  ",
            "ref>alt   : C[2]>C[3]",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_dup_1bp(self):
        hgvs_g = "NC_000007.13:g.36561663dup"
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561663dup",
            "hgvs_c    : NM_001177507.2:c.1485dup",
            "hgvs_p    : NP_001170978.1:p.(Gly496TrpfsTer39)",
            "          :        36,561,650          36,561,670",
            "chrom pos :   .    |    .    |    .    |    .    |   ",
            "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG",
            "region    :                     |                    ",
            "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC",
            "tx pos    :      |    .    |    .    |    .    |    .",
            "          :      1500      1490      1480      1470",
            "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln",
            "aa pos    :      |||            ...            |||   ",
            "          :      500                           490   ",
            "ref>alt   : A>AA",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_identity(self):
        hgvs_g = "NC_000007.13:g.36561663="
        var_g = self.hp.parse(hgvs_g)

        result = self.pp.display(var_g, "NM_001177507.2")
        print(result)

        result = result.split("\n")
        expected_str = [
            "hgvs_g    : NC_000007.13:g.36561663=",
            "hgvs_c    : NM_001177507.2:c.1485=",
            "hgvs_p    : NP_001170978.1:p.(Asp495=)",
            "          :        36,561,650          36,561,670",
            "chrom pos :   .    |    .    |    .    |    .    |   ",
            "seq    -> : ACCTCGTTGGGGTGGAATCCATCCACGGGCTCGATGAGCTG",
            "region    :                     =                    ",
            "tx seq <- :   GAGCAACCCCACCTTAGGTAGGTGCCCGAGCTACTCGAC",
            "tx pos    :      |    .    |    .    |    .    |    .",
            "          :      1500      1490      1480      1470",
            "aa seq <- :   GluAsnProHisPheGlyAspValProGluIleLeuGln",
            "aa pos    :      |||            ...            |||   ",
            "          :      500                           490   ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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

        expected_str = [
            "hgvs_g    : NC_000005.10:g.123346517_123346518insATTA",
            "hgvs_c    : NM_001166226.1:c.*1_*2insTAAT",
            "hgvs_p    : NP_001159698.1:p.?",
            "          :   123,346,520",
            "chrom pos :   |    .",
            "seq    -> : ATTAATTA",
            "region    : |------|",
            "tx seq <- : TAATTAAT",
            "tx pos    : |   |   ",
            "          : *1  2880",
            "aa seq <- :  TerAsnS",
            "aa pos    :     ||| ",
            "          :     960 ",
            "ref>alt   : ATTA[2]>ATTA[3]",
        ]

        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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
        expected_str = [
            "NC_000012.11:g.33049660_33049680dup",
            "NM_004572.3:c.-9_12dup",
            "NP_004563.2:p.(Met1_Pro4dup)",
            "      33,049,650          33,049,670          33,049,690          33,049,710          33,049,730          33,049,750          33,049,770          33,049,790",
            " .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |",
            "CTGGGGCGCCGGGGGCTGCCATGGGGCCGGTGGGGGCGACCGAGCTGCTCGCCTGCCTCTGGACTCGCGGGCGAAGCCGCCACGGAGCTGGGGGCGCTGGCGCGAGCCCCGCCCCGCTCGAGTCCGGCCCCGCCCCTGGCCCGCCCC",
            "          |-------------------------|                                                                                                              ",
            "GACCCCGCGGCCCCCGACGGTAccccggccacccccgctggctcgacgagcggacggagacctgagcgcccgcttcggcggtgcctcgacccccgcgaccgcgctcggggcggggcgagctcaggccggggcgggga          ",
            "  |    .    |    .   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .          ",
            "  20        10       1         -10       -20       -30       -40       -50       -60       -70       -80       -90       -100      -110",
            "aProAlaGlyProAlaAlaMet                                                                                                                             ",
            "       ...                                                                                                                                         ",
            "                   1                                                                                                                               ",
            "GGGGGCTGCCATGGGGCCGGTGGGGGC>GGGGGCTGCCATGGGGCCGGTGGGGGCTGCCATGGGGCCGGTGGGGGC",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_ref_disagree(self):
        """Test a tx ref disagree variant."""
        hgvs_g = "NM_001111.4:c.298G>A"
        var_g = self.hp.parse(hgvs_g)
        result = self.pp.display(var_g)

        print(result)

        result = result.split("\n")

        # note the X in the transscript sequence
        expected_str = [
            "hgvs_g    : NC_000001.10:g.154574820=",
            "hgvs_c    : NM_001111.4:c.298G>A",
            "hgvs_p    : NP_001102.2:p.(Gly100Arg)",
            "          : 154,574,800         154,574,820         154,574,840",
            "chrom pos : |    .    |    .    |    .    |    .    |",
            "seq    -> : TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA",
            "tx ref dif:                     X                    ",
            "region    :                     =                    ",
            "tx seq <- : AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT",
            "tx pos    :    .    |    .    |    .    |    .    |  ",
            "          :         310       300       290       280",
            "aa seq <- : ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa",
            "aa pos    :    ...            |||            ...     ",
            "          :                   100                    ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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
        expected_str = [
            "hgvs_g    : NC_000021.8:g.46020668_46020682del",
            "hgvs_c    : NM_198689.2:c.124_135=",
            "hgvs_p    : NP_941962.1:p.(Cys42=)",
            "          :     46,020,630          46,020,650          46,020,670          46,020,690          46,020,710",
            "chrom pos :     |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  ",
            "seq    -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCCCTGCTGCGCCCCCAGCTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCCGT",
            "tx ref dif:                               IIIIIIIIIIIIIII                                                 XX ",
            "region    :                               x-------------------------x                                        ",
            "tx seq -> : CGACTGCCCAGAGAGCTGCTGCGAGCCCCC---------------CTGCTGCGCCCCGGCCCCCTGCCTGAGCCTGGTCTGCACCCCAGTGAGCTAT",
            "tx pos    : .    |    .    |    .    |                   .    |    .    |    .    |    .    |    .    |    . ",
            "          :      110       120       130                      140       150       160       170       180",
            "aa seq -> : pAspCysProGluSerCysCysGluProPr---------------oCysCysAlaProAlaProCysLeuSerLeuValCysThrProValSerTyr",
            "aa pos    : .            |||            ..               .            |||            ...            |||      ",
            "          :              40                                           50                            60       ",
            "ref>alt   : CTGCTGCGCCCCCAGCTGCTGCGCCCC>CTGCTGCGCCCC",
            "",
        ]

        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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

        expected_str = [
            "hgvs_g    : NC_000002.11:g.96780987_96780997del",
            "hgvs_c    : NM_000682.6:c.901_911del",
            "hgvs_p    : NP_000673.2:p.(Glu301GlyfsTer7)",
            "          :     96,780,960          96,780,980                   96,781,000          96,781,020",
            "chrom pos :     |    .    |    .    |    .    |    .  _________  |    .    |    .    |    .    |    .  ",
            "seq    -> : TGCCTGGGGTTCACACTCTTCCTCCTCCTCCTCCTCCTCTTC.........AGCTTCATCCTCTGGAGATGCCCCACAAACACCCTCCTTC",
            "tx ref dif:                                           DDDDDDDDD                                        ",
            "region    :                               x----------x                                                 ",
            "tx seq <- : ACGGACCCCAAGTGTGAGAAGGAGGAGGAGGAGGAGGAGAAGGAGGAGAAGTCGAAGTAGGAGACCTCTACGGGGTGTTTGTGGGAGGAAG",
            "tx pos    :   |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .   ",
            "          :   940       930       920       910       900       890       880       870       860",
            "aa seq <- : AlaGlnProGluCysGluGluGluGluGluGluGluGluGluGluGluGluAlaGluAspGluProSerAlaGlyCysValGlyGluLysG",
            "aa pos    :             |||            ...            |||            ...            |||            ... ",
            "          :             310                           300                           290                ",
            "ref>alt   : CTCCTCCTCTTC>C",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

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

        expected_str = [
            "hgvs_g    : NC_000002.11:g.96780987_96780997del",
            "hgvs_c    : NM_000682.6:c.901_911del",
            "hgvs_p    : NP_000673.2:p.(Glu301GlyfsTer7)",
            "          :        96,781,030          96,781,010                   96,780,990          96,780,970",
            "chrom pos :   .    |    .    |    .    |    .    |  _________  .    |    .    |    .    |    .    |    ",
            "seq    <- : CTTCCTCCCACAAACACCCCGTAGAGGTCTCCTACTTCGA.........CTTCTCCTCCTCCTCCTCCTCCTTCTCACACTTGGGGTCCGT",
            "seq    -> : GAAGGAGGGTGTTTGTGGGGCATCTCCAGAGGATGAAGCT.........GAAGAGGAGGAGGAGGAGGAGGAAGAGTGTGAACCCCAGGCA",
            "tx ref dif:                                         DDDDDDDDD                                          ",
            "region    :                                                  x----------x                              ",
            "tx seq -> : GAAGGAGGGTGTTTGTGGGGCATCTCCAGAGGATGAAGCTGAAGAGGAGGAAGAGGAGGAGGAGGAGGAGGAAGAGTGTGAACCCCAGGCA",
            "tx pos    :    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  ",
            "          :         860       870       880       890       900       910       920       930       940",
            "aa seq -> : nLysGluGlyValCysGlyAlaSerProGluAspGluAlaGluGluGluGluGluGluGluGluGluGluGluGluCysGluProGlnAla",
            "aa pos    :  ...            |||            ...            |||            ...            |||            ",
            "          :                 290                           300                           310            ",
            "ref>alt   : CTCCTCCTCTTC>C",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    @pytest.mark.skip(
        reason="actually not that special, but still a nice variant since there is a large shuffle-able sequence on both ends."
    )
    def test_exon_boundary_overlap_forward_strand(self):
        hgvs_c = "NM_001283009.2:c.1228_1266+39del"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, self.assembly_mapper37, show_legend=True, use_color=False)

        result = pp.display(var_c)

        print(result)

    def test_ruler(self):
        """Test the ruler display option turned on."""
        hgvs_c = "NM_001111.4:c.298G>A"
        var_c = self.hp.parse(hgvs_c)
        pp = PrettyPrint(self.hdp, self.assembly_mapper37, show_legend=False, reverse_display=False)

        result = pp.display(var_c, "NM_001111.4")

        print(result)

        result = result.split("\n")

        # note the X in the transscript ref-disagree row
        expected_str = [
            "NC_000001.10:g.154574820=",
            "NM_001111.4:c.298G>A",
            "NP_001102.2:p.(Gly100Arg)",
            "154,574,800         154,574,820         154,574,840",
            "|    .    |    .    |    .    |    .    |",
            "TCTCTGGAGCCCCTGACTTCTGAGATGCACGCCCCTGGGGA",
            "                    X                    ",
            "                    =                    ",
            "AGAGACCTCGGGGACTGAAGGCTCTACGTGCGGGGACCCCT",
            "   .    |    .    |    .    |    .    |  ",
            "        310       300       290       280",
            "ArgGlnLeuGlyGlnSerGlyLeuHisValGlyArgProVa",
            "   ...            |||            ...     ",
            "                  100                    ",
            "",
        ]
        for r, e in zip(result, expected_str, strict=False):
            assert e == r

    def test_rna_coding(self):
        """a rna coding transcript."""
        hgvs_n = "NR_146230.2:n.10G>A"
        var_n = self.hp.parse(hgvs_n)

        pp = PrettyPrint(self.hdp, self.assembly_mapper37, show_reverse_strand=True)
        result = pp.display(var_n)
        print(result)
        expected_str = (
            "hgvs_g    : NC_000001.10:g.167905930G>A\n"
            "hgvs_n    : NR_146230.2:n.10G>A\n"
            "          : 167,905,910         167,905,930         167,905,950\n"
            "chrom pos : |    .    |    .    |    .    |    .    |\n"
            "seq    -> : TGCTGATCTTTGGATGTTCTGGTTAGTCTAAGAAGGAGAGT\n"
            "seq    <- : ACGACTAGAAACCTACAAGACCAATCAGATTCTTCCTCTCA\n"
            "region    :                     A                    \n"
            "tx seq -> :            GGATGTTCTGGTTAGTCTAAGAAGGAGAGT\n"
            "tx pos    :                .    |    .    |    .    |\n"
            "          :                     10        20        30\n"
            "aa seq -> :                                          \n"
            "aa pos    :                                          \n"
            "          :                                          \n"
        )
        for r, e in zip(result, expected_str, strict=False):
            assert e == r
