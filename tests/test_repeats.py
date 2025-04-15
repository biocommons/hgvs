import pytest
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig, VariantCoords
from hgvs.repeats import RepeatAnalyser, detect_repetitive_block_lengths, get_repeat_str


def to_hgvs_repeat(chrom_ac, fs: VariantCoords, ra: RepeatAnalyser):
    """Experimental method to see how close we can get to HGVS repeat nomenclature.
    Note: In many cases we are not yet creating valid HGVS strings here. However this works at least for some cases.
    (homoploymeric repeats). If a repeat can be shuffled this needs some modification to use the right-aligned
    representation as well as trimming off bases that can be shuffled (they come out as repeats of length 1 at the moment.)
    """
    return f"{chrom_ac}:g.{fs.start + 1}_{fs.end}{ra.alt_str}"


@pytest.mark.skip(reason="would add too much caching burden.")
class TestHGVSExamples:
    """The variants in this class are taken from the HGVS nomenclature page about repeats
    https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/
    """

    def test_intergenic_repeat(self, parser, hdp):
        # we can't parse HGVS repeat representations yet:
        # NC_000014.8:g.101179660_101179695TG[14]
        # We have to to use an alternative way to describe this variant:
        hgvs_g = "NC_000014.8:g.101179660_101179667del"

        var_g = parser.parse_hgvs_variant(hgvs_g)

        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)

        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        assert ra.is_repeat
        assert ra.ref_str == "TG[18]"
        assert ra.alt_str == "TG[14]"
        assert (
            to_hgvs_repeat(var_g.ac, fs, ra)
            == "NC_000014.8:g.101179660_101179695TG[14]"
        )  # a valid HGVS name

    def test_CACNA1_repeat(self, parser, hdp, am38):
        """NM_023035.2:c.6955_6993CAG[26] (or c.6955_6993dup )
        A repeated CAG tri-nucleotide sequence (CAG on reverse strand, CTG on forward strand).
        """
        hgvs_c = "NM_023035.2:c.6955_6993dup"
        var_c = parser.parse_hgvs_variant(hgvs_c)
        var_g = am38.c_to_g(var_c)
        assert str(var_g) == "NC_000019.10:g.13207860_13207898dup"

        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)

        fs = dc.get_shuffled_variant(var_g, 0)

        ra = RepeatAnalyser(fs)
        assert ra.is_repeat
        assert ra.ref_str == "CTG[13]C[1]"
        assert ra.alt_str == "CTG[26]C[1]"
        assert (
            to_hgvs_repeat(var_g.ac, fs, ra)
            == "NC_000019.10:g.13207859_13207898CTG[26]C[1]"
        )  # not valid. It is 5' shuffled and there's a C that can be shuffled which needs to get trimmed.

        # This repeat is shuffle-able. These are alternate representations of each other:
        # CTG[13]C[1]
        # C[1]TGC[13]

        # here we test the alternative shuffled representation:
        ra2 = RepeatAnalyser(fs, reverse=True)
        assert ra2.is_repeat
        assert ra2.ref_str == "C[1]TGC[13]"
        assert ra2.alt_str == "C[1]TGC[26]"
        assert (
            to_hgvs_repeat(var_g.ac, fs, ra2)
            == "NC_000019.10:g.13207859_13207898C[1]TGC[26]"
        )  # not valid, we should trim off the shuffle-able C for that.

    def test_ATXN7_repeat(self, parser, hdp):
        """A repeated AGC tri-nucleotide sequence in the ATXN7 gene on chromosome 3

        Note: In literature, the tri-nucleotide repeat, encoding a poly-Gln repeat on protein level,
        is known as the CAG repeat (CAG is the codon for Gln).

        This repeat is shuffle-able.

        There are three alternative representations that are all equally valid for this repeat:
        GCA[10]G[1]C[1]
        G[1]CAG[10]C[1]
        G[1]C[1]AGC[10]

        At the moment we can create the first and third representation of this repeat, but not yet the middle one. Note: the middle one is the one
        that would match the codons in the transcript, and prob make most sense from a biological perspective.
        """

        hgvs_g = "NC_000003.12:g.63912688_63912689insCAGCAGCAG"
        var_g = parser.parse_hgvs_variant(hgvs_g)

        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)

        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        assert ra.is_repeat
        assert ra.ref_str == "GCA[10]G[1]C[1]"
        assert ra.alt_str == "GCA[13]G[1]C[1]"

        # invalid HGVS, 5' aligned, and shows the two shuffle-able bases at the end.
        assert (
            to_hgvs_repeat(var_g.ac, fs, ra)
            == "NC_000003.12:g.63912685_63912716GCA[13]G[1]C[1]"
        )

        ra2 = RepeatAnalyser(fs, reverse=True)
        assert ra2.is_repeat
        assert ra2.ref_str == "G[1]C[1]AGC[10]"
        assert ra2.alt_str == "G[1]C[1]AGC[13]"

        # invalid HGVS, due to the shuffle-able bases, but comes pretty close to what HGVS recommends.
        assert (
            to_hgvs_repeat(var_g.ac, fs, ra2)
            == "NC_000003.12:g.63912685_63912716G[1]C[1]AGC[13]"
        )


@pytest.mark.skip(reason="would add too much caching burden.")
@pytest.mark.parametrize(
    "hgvs_g, is_repeat, repeat_unit, ref_count, alt_count, s",
    [
        ("NC_000019.10:g.45770205del", True, "C", 6, 5, "C[6]>C[5]"),
        ("NC_000019.10:g.45770204_45770205del", True, "C", 6, 4, "C[6]>C[4]"),
        ("NC_000019.10:g.45770205insC", True, "C", 6, 7, "C[6]>C[7]"),
        ("NC_000019.10:g.45770206del", False, None, 0, 0, "A>"),
        ("NC_000007.13:g.36561662dup", True, "C", 2, 3, "C[2]>C[3]"),
    ],
)
class TestHomopolymerRepeats:
    def test_homopolymer(
        self, parser, hdp, hgvs_g, is_repeat, repeat_unit, ref_count, alt_count, s
    ):
        var_g = parser.parse_hgvs_variant(hgvs_g)
        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)
        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        assert s == str(ra)
        assert is_repeat == ra.is_repeat

        if repeat_unit is None:
            assert len(ra.repeat_units_alt) == 0
        else:
            assert repeat_unit == ra.repeat_units_alt[0].repeat_unit

        if ref_count == 0:
            assert len(ra.repeat_units_ref) == 0
        else:
            assert ref_count == ra.repeat_units_ref[0].repeat_count

        if alt_count == 0:
            assert len(ra.repeat_units_alt) == 0
        else:
            assert alt_count == ra.repeat_units_alt[0].repeat_count


@pytest.mark.skip(reason="would add too much caching burden.")
@pytest.mark.parametrize(
    "hgvs_g, is_repeat, repeat_unit, ref_count, alt_count, str_repr",
    [
        (
            "NC_000021.8:g.46020668_46020682del",
            False,
            "CTG",
            2,
            2,
            "CTGCTGCGCCCCCAGCTGCTGCGCCCC>CTGCTGCGCCCC",
        ),
        (
            "NC_000012.11:g.33049660_33049680dup",
            False,
            "G",
            5,
            5,
            "GGGGGCTGCCATGGGGCCGGTGGGGGC>GGGGGCTGCCATGGGGCCGGTGGGGGCTGCCATGGGGCCGGTGGGGGC",
        ),
        (
            "NC_000005.10:g.123346517_123346518insATTA",
            True,
            "ATTA",
            2,
            3,
            "ATTA[2]>ATTA[3]",
        ),
        (
            "NC_000005.10:g.123346522_123346525dup",
            True,
            "ATTA",
            2,
            3,
            "ATTA[2]>ATTA[3]",
        ),
        (
            "NC_000019.10:g.45770210_45770212del",
            True,
            "CAG",
            20,
            19,
            "CAG[20]C[1]A[1]>CAG[19]C[1]A[1]",
        ),  # note this one can be shuffled CAG/GCA
        ("NC_000007.14:g.117548628_117548629insTTTT", True, "T", 7, 11, "T[7]>T[11]"),
        ("NC_000009.11:g.35079521_35079523del", True, "TGG", 2, 1, "TGG[2]>TGG[1]"),
        (
            "NC_000001.11:g.6490477_6490484del",
            True,
            "TCTAAGGC",
            2,
            1,
            "TCTAAGGC[2]T[1]C[1]>TCTAAGGC[1]T[1]C[1]",
        ),
    ],
)
class TestRepeats:
    def test_repeats(
        self,
        parser,
        hdp,
        hgvs_g,
        is_repeat,
        repeat_unit,
        ref_count,
        alt_count,
        str_repr,
    ):
        var_g = parser.parse_hgvs_variant(hgvs_g)
        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)
        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)

        assert str_repr == str(ra)

        # note: in some cases is_repeat is False, although repeat_unit contains some data.
        # however in those cases the repeat patter was not considered significant enough to name the whole variant accordingly
        assert is_repeat == ra.is_repeat
        assert repeat_unit == ra.repeat_units_alt[0].repeat_unit
        assert ref_count == ra.repeat_units_ref[0].repeat_count
        assert alt_count == ra.repeat_units_alt[0].repeat_count


@pytest.mark.skip(reason="would add too much caching burden.")
class TestLargeVariant:
    def test_large_variant(self, parser, hdp):
        """Large variants would take too much time to analyze, plus the results are not that useful."""

        hgvs_g = "NC_000001.11:g.9976249_9983617dup"
        var_g = parser.parse_hgvs_variant(hgvs_g)
        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)
        fs = dc.get_shuffled_variant(var_g, 0)
        ra = RepeatAnalyser(fs)
        assert not ra.is_repeat


class TestRepeatMethods:
    @pytest.mark.parametrize(
        "ref, alt, ref_count, ref_repeat_unit, alt_count, alt_repeat_unit, mixed_repeat_ref, mixed_repeat_alt, expected_ref_repeat, expected_alt_repeat",
        [
            ("T", "TTTTTTT", None, None, 7, "T", False, False, "T[1]", "T[7]"),
            ("TTTTTTT", "T", 7, "T", None, None, False, False, "T[7]", "T[1]"),
            ("T", "TT", None, None, 2, "T", False, False, "T[1]", "T[2]"),
            ("TT", "T", 2, "T", None, None, False, False, "T[2]", "T[1]"),
            ("TT", "TT", 2, "T", 2, "T", False, False, "T[2]", "T[2]"),
            ("", "TT", None, None, 2, "T", False, False, None, None),
            ("TT", "", 2, "T", None, None, False, False, None, None),
            ("TT", "TTTTTTT", 2, "T", 7, "T", False, False, "T[2]", "T[7]"),
            ("GA", "GAGA", None, None, 2, "GA", False, False, "GA[1]", "GA[2]"),
            ("TTTT", "TTTTTTT", 4, "T", 7, "T", False, False, "T[4]", "T[7]"),
            ("abc", "abc", None, None, None, None, False, False, None, None),
            ("abcabc", "abcab", 2, "abc", None, None, False, False, "abc[2]", None),
            ("", "", None, None, None, None, False, False, None, None),
            ("GA", "GCGC", None, None, 2, "GC", False, False, None, "GC[2]"),
            ("GCGC", "GA", 2, "GC", None, None, False, False, "GC[2]", None),
            (
                "GAGAAAA",
                "GAGAAAAA",
                2,
                "GA",
                3,
                "A",
                True,
                True,
                "GA[2]A[3]",
                "GA[2]A[4]",
            ),
            (
                "GAGAAAA",
                "GAGAA",
                4,
                "A",
                2,
                "GA",
                True,
                False,
                "GA[2]A[3]",
                "GA[2]A[1]",
            ),
            (
                "GAGAGAABCABCABCGA",
                "ABCABC",
                4,
                "A",
                2,
                "ABC",
                True,
                False,
                "GA[3]ABC[3]G[1]A[1]",
                "ABC[2]",
            ),
            (
                "GAGAGAABCABCABCTGAGA",
                "ABCABC",
                4,
                "A",
                2,
                "ABC",
                True,
                False,
                "GA[3]ABC[3]T[1]GA[2]",
                "ABC[2]",
            ),
            (
                "ABBCABBCABBC",
                "ABBCABBC",
                3,
                "ABBC",
                2,
                "ABBC",
                False,
                False,
                "ABBC[3]",
                "ABBC[2]",
            ),
        ],
    )
    def test_detect_repetitive_block_lengths(
        self,
        ref,
        alt,
        ref_count,
        ref_repeat_unit,
        alt_count,
        alt_repeat_unit,
        mixed_repeat_ref,
        mixed_repeat_alt,
        expected_ref_repeat,
        expected_alt_repeat,
    ):
        repeat_units_ref = detect_repetitive_block_lengths(ref)
        repeat_units_alt = detect_repetitive_block_lengths(alt)

        if ref_count is None:
            assert len(repeat_units_ref) == 0
        if alt_count is None:
            assert len(repeat_units_alt) == 0

        if len(repeat_units_ref) == 1:
            ref_u = repeat_units_ref[0].repeat_unit
            ref_l = repeat_units_ref[0].repeat_count
            assert ref_u == ref_repeat_unit
            assert ref_l == ref_count
            assert not mixed_repeat_ref
        elif len(repeat_units_ref) > 1:
            # more than one repeat found.
            assert mixed_repeat_ref

        if len(repeat_units_alt) == 1:
            alt_u = repeat_units_alt[0].repeat_unit
            alt_l = repeat_units_alt[0].repeat_count
            assert alt_u == alt_repeat_unit
            assert alt_l == alt_count
            assert not mixed_repeat_alt
        elif len(repeat_units_alt) > 1:
            assert mixed_repeat_alt

        # next: rebuild ref and alt strings:
        ref_repeat = get_repeat_str(ref, alt, repeat_units_ref, repeat_units_alt)
        assert ref_repeat == expected_ref_repeat

        alt_repeat = get_repeat_str(alt, ref, repeat_units_alt, repeat_units_ref)
        assert alt_repeat == expected_alt_repeat
