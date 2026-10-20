import random
from types import SimpleNamespace

import pytest
from bioutils.sequences import reverse_complement

from hgvs.extras.babelfish import SYMBOLIC_ALTS, Babelfish, VCFCoordinate, _common_prefix_len

SYMBOLIC_TYPES = {v: k for k, v in SYMBOLIC_ALTS.items()}
CHR2 = "NC_000002.12"

NORM_HGVS_VCF = [
    # Columns are: (normed-HGVS, non-normalized HGVS, VCF coordinates, non-norm VCF)
    # no-op
    (
        "NC_000006.12:g.49949407=",
        [],
        ("6", 49949407, "A", ".", "identity"),
        [
            ("6", 49949407, "A", "A", "identity"),
            # Test case insensitivity
            ("6", 49949407, "A", "a", "identity"),
            ("6", 49949407, "a", "A", "identity"),
        ],
    ),
    # Test multi-base identity
    (
        "NC_000006.12:g.49949407_49949408=",
        [],
        ("6", 49949407, "AA", ".", "identity"),
        [("6", 49949407, "AA", "AA", "identity")],
    ),
    # snv
    (
        "NC_000006.12:g.49949407A>T",
        [],
        # was ("6", 49949406, "AA", "AT", "sub") however VT parsimony rules say it should be those below
        ("6", 49949407, "A", "T", "sub"),
        [],
    ),
    # delins
    (
        "NC_000006.12:g.49949413_49949414delinsCC",
        [],
        # This was ("6", 49949412, "AAA", "ACC", "delins") - however VT parsimony rules say it should be those below
        ("6", 49949413, "AA", "CC", "delins"),
        [],
    ),
    # del, no shift
    ("NC_000006.12:g.49949415del", [], ("6", 49949414, "AT", "A", "del"), []),
    # del, w/ shift
    (
        "NC_000006.12:g.49949414del",
        ["NC_000006.12:g.49949413del"],
        ("6", 49949409, "GA", "G", "del"),
        [],
    ),
    ("NC_000006.12:g.49949413_49949414del", [], ("6", 49949409, "GAA", "G", "del"), []),
    # ins, no shift
    ("NC_000006.12:g.49949413_49949414insC", [], ("6", 49949413, "A", "AC", "ins"), []),
    ("NC_000006.12:g.49949414_49949415insCC", [], ("6", 49949414, "A", "ACC", "ins"), []),
    # ins/dup, w/shift
    (
        "NC_000006.12:g.49949414dup",
        ["NC_000006.12:g.49949413_49949414insA", "NC_000006.12:g.49949414_49949415insA"],
        ("6", 49949409, "G", "GA", "dup"),
        [],
    ),
    (
        "NC_000006.12:g.49949413_49949414dup",
        ["NC_000006.12:g.49949414_49949415insAA"],
        ("6", 49949409, "G", "GAA", "dup"),
        [],
    ),
]


@pytest.mark.extra
@pytest.mark.vcr
def test_hgvs_to_vcf(parser, babelfish38):
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """

    def _h2v(h):
        return babelfish38.hgvs_to_vcf(parser.parse(h))

    for norm_hgvs_string, alt_hgvs, expected_variant_coordinate, _ in NORM_HGVS_VCF:
        for hgvs_string in [norm_hgvs_string, *alt_hgvs]:
            variant_coordinates = _h2v(hgvs_string)
            assert variant_coordinates == expected_variant_coordinate


def test_vcf_to_hgvs(babelfish38):
    def _v2h(*v):
        return babelfish38.vcf_to_g_hgvs(*v)

    for expected_hgvs_string, _, norm_variant_coordinate, alt_variant_coordinate in NORM_HGVS_VCF:
        for variant_coordinate in [norm_variant_coordinate, *alt_variant_coordinate]:
            *v, _typ = variant_coordinate  # last column is type ie "dup"
            hgvs_g = _v2h(*v)
            hgvs_string = hgvs_g.format()
            assert hgvs_string == expected_hgvs_string


def test_vcf_to_hgvs_contig_chrom(babelfish38):
    hgvs_g = babelfish38.vcf_to_g_hgvs("NC_000006.12", 49949409, "GAA", "G")
    assert hgvs_g.format() == "NC_000006.12:g.49949413_49949414del"


class _NormalizeCalledError(Exception):
    pass


def _anchor_babelfish(hdp, symbolic_alt_min_length, normalize_max_length, **kwargs):
    babelfish = Babelfish(
        hdp,
        assembly_name="GRCh38",
        symbolic_alt_min_length=symbolic_alt_min_length,
        normalize_max_length=normalize_max_length,
        **kwargs,
    )

    def _normalize(var):
        raise _NormalizeCalledError(var)

    babelfish.hn = SimpleNamespace(normalize=_normalize)
    return babelfish


@pytest.fixture
def symbolic_babelfish(hdp):
    return _anchor_babelfish(hdp, symbolic_alt_min_length=1000, normalize_max_length=1000)


SYMBOLIC_HGVS_VCF = [
    # (HGVS, VCF coordinates, VCF END)
    ("NC_000002.12:g.1000000_224225011dup", ("2", 999999, "T", "<DUP>", "dup"), 224225011),
    ("NC_000002.12:g.1000000_224225011del", ("2", 999999, "T", "<DEL>", "del"), 224225011),
    ("NC_000002.12:g.1000000_224225011inv", ("2", 999999, "T", "<INV>", "inv"), 224225011),
]


def test_hgvs_to_vcf_symbolic(parser, symbolic_babelfish):
    for hgvs_string, expected_variant_coordinate, end in SYMBOLIC_HGVS_VCF:
        var_g = parser.parse(hgvs_string)
        assert symbolic_babelfish.hgvs_to_vcf(var_g) == expected_variant_coordinate
        vc = symbolic_babelfish.hgvs_to_vcf_coordinate(var_g)
        _chrom, position, ref, alt, _typ = expected_variant_coordinate
        assert vc.as_tuple() == (CHR2, position, ref, alt, end)
        assert symbolic_babelfish.vcf_chrom(vc) == "2"
        assert vc.svlen == symbolic_babelfish.symbolic_alt_length(var_g)
        svtype = alt[1:-1]
        svlen = -vc.svlen if svtype == "DEL" else vc.svlen  # VCF 4.3
        assert symbolic_babelfish.vcf_info(vc) == {"SVTYPE": svtype, "SVLEN": svlen, "END": end}
        assert str(vc) == f"{CHR2}:{position}-{end} {alt}"
        assert symbolic_babelfish.vcf_coordinate_to_g_hgvs(vc).format() == hgvs_string


def test_hgvs_to_vcf_chrom_style_accession(parser, hdp):
    babelfish = _anchor_babelfish(hdp, 1000, 1000, chrom_style="accession")
    for hgvs_string, (_chrom, *rest), _end in SYMBOLIC_HGVS_VCF:
        var_g = parser.parse(hgvs_string)
        assert babelfish.hgvs_to_vcf(var_g) == (CHR2, *rest)
        assert babelfish.vcf_chrom(babelfish.hgvs_to_vcf_coordinate(var_g)) == CHR2
    with pytest.raises(ValueError, match="chrom_style"):
        Babelfish(hdp, "GRCh38", chrom_style="chr")


def test_hgvs_to_vcf_requires_g(parser, symbolic_babelfish):
    with pytest.raises(RuntimeError, match="Expected g"):
        symbolic_babelfish.hgvs_to_vcf(parser.parse("NM_000059.4:c.100del"))


def test_hgvs_to_vcf_symbolic_explicit_sequence(parser, symbolic_babelfish):
    var_g = parser.parse("NC_000002.12:g.1000000_1001999del" + "A" * 2000)
    assert symbolic_babelfish.hgvs_to_vcf(var_g) == ("2", 999999, "T", "<DEL>", "del")


def test_hgvs_to_vcf_symbolic_uncertain_uses_inner_interval(parser, symbolic_babelfish):
    var_g = parser.parse("NC_000002.12:g.(999000_1000000)_(224225011_224226000)dup")
    assert symbolic_babelfish.hgvs_to_vcf(var_g) == ("2", 999999, "T", "<DUP>", "dup")


def test_hgvs_to_vcf_normalized(parser, hdp):
    normalized = [
        ("NC_000002.12:g.1000000_1000998dup", 1000, 1000),  # below symbolic_alt_min_length
        ("NC_000002.12:g.1000000_1001999dup", 1000, 10000),  # below normalize_max_length
        ("NC_000002.12:g.1000000_1009998dup", 1000, None),  # below default normalize_max_length
        ("NC_000002.12:g.1000000_1000003delGG", 1, 1),  # sequence doesn't fit interval
        ("NC_000002.12:g.1000000_1000003delinsA", 1, 1),
        ("NC_000002.12:g.1000000_1000998dup", None, None),  # not requested
    ]
    for hgvs_string, symbolic_alt_min_length, normalize_max_length in normalized:
        babelfish = _anchor_babelfish(hdp, symbolic_alt_min_length, normalize_max_length)
        var_g = parser.parse(hgvs_string)
        with pytest.raises(_NormalizeCalledError):
            babelfish.hgvs_to_vcf(var_g)


def test_normalize_max_length_requires_symbolic(hdp):
    for symbolic_alt_min_length in (None, 1000):
        with pytest.raises(ValueError, match="normalize_max_length"):
            Babelfish(hdp, "GRCh38", symbolic_alt_min_length, normalize_max_length=50)


def test_normalize_max_length_default(hdp):
    for symbolic_alt_min_length, normalize_max_length in [(50, 1_000_000), (2_000_000, 2_000_000)]:
        babelfish = Babelfish(hdp, "GRCh38", symbolic_alt_min_length)
        assert babelfish.normalize_max_length == normalize_max_length


def test_vcf_to_hgvs_symbolic(symbolic_babelfish):
    for expected_hgvs_string, variant_coordinate, end in SYMBOLIC_HGVS_VCF:
        chrom, position, ref, alt, _ = variant_coordinate
        hgvs_g = symbolic_babelfish.vcf_to_g_hgvs(chrom, position, ref, alt, end=end)
        assert hgvs_g.format() == expected_hgvs_string


def test_vcf_to_hgvs_symbolic_svlen(symbolic_babelfish):
    expected = "NC_000002.12:g.1000000_224225011del"
    for svlen in (223225012, -223225012):  # VCF 4.3 deletions have negative SVLEN
        hgvs_g = symbolic_babelfish.vcf_to_g_hgvs("2", 999999, "T", "<DEL>", svlen=svlen)
        assert hgvs_g.format() == expected


def test_vcf_to_hgvs_symbolic_requires_valid_end(symbolic_babelfish):
    with pytest.raises(ValueError, match="requires end"):
        symbolic_babelfish.vcf_to_g_hgvs("2", 999999, "T", "<DUP>")
    with pytest.raises(ValueError, match="must be after position"):
        symbolic_babelfish.vcf_to_g_hgvs("2", 999999, "T", "<DEL>", end=999999)
    with pytest.raises(ValueError, match="Unsupported"):
        symbolic_babelfish.vcf_to_g_hgvs("2", 999999, "T", "<CNV>", end=1000999)


class _SyntheticHDP:
    """A random contig, so normalization runs without UTA/SeqRepo. The only repeat is
    TANDEM_REPEAT, 2 copies of 301_360 (1-based)"""

    TANDEM_REPEAT = (300, 360, 420)  # 0-based start of copy 1, copy 2, end

    def __init__(self, length=1000):
        rng = random.Random(0)  # noqa: S311
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        start, middle, end = self.TANDEM_REPEAT
        self.seq = seq[:middle] + seq[start:middle] + seq[end:]

    def get_seq(self, ac, start_i=None, end_i=None):  # noqa: ARG002
        return self.seq[start_i:end_i]


@pytest.fixture
def synthetic_babelfish():
    return Babelfish(
        _SyntheticHDP(),
        assembly_name="GRCh38",
        symbolic_alt_min_length=50,
        normalize_max_length=100,
    )


def test_vcf_to_hgvs_small_symbolic_matches_explicit(synthetic_babelfish):
    seq = synthetic_babelfish.hdp.seq
    position, end = 100, 120  # 20bp < symbolic_alt_min_length
    pad = seq[position - 1]
    span = seq[position:end]
    explicit = {
        "<DEL>": (position, pad + span, pad),
        "<DUP>": (position, pad, pad + span),
        "<INV>": (position + 1, span, reverse_complement(span)),
    }
    for symbolic_alt, (explicit_position, ref, alt) in explicit.items():
        hgvs_g = synthetic_babelfish.vcf_to_g_hgvs("2", position, pad, symbolic_alt, end=end)
        assert hgvs_g.posedit.edit.type == SYMBOLIC_TYPES[symbolic_alt]
        assert hgvs_g == synthetic_babelfish.vcf_to_g_hgvs("2", explicit_position, ref, alt)


def test_vcf_to_hgvs_large_explicit_as_symbolic(synthetic_babelfish):
    seq = synthetic_babelfish.hdp.seq
    position, end = 500, 620  # 120bp >= normalize_max_length
    pad = seq[position - 1]
    span = seq[position:end]
    preceding = seq[position - 120 : position]
    explicit_symbolic = [
        ((position, pad + span, pad), "<DEL>"),
        ((position, pad, pad + span), "<DUP>"),  # left aligned, copy before the original
        ((end, seq[end - 1], seq[end - 1] + span), "<DUP>"),  # copy after the original
        ((position + 1, span, reverse_complement(span)), "<INV>"),
    ]
    for (explicit_position, ref, alt), symbolic_alt in explicit_symbolic:
        hgvs_g = synthetic_babelfish.vcf_to_g_hgvs("2", explicit_position, ref, alt)
        assert (
            hgvs_g.format() == f"NC_000002.12:g.{position + 1}_{end}{SYMBOLIC_TYPES[symbolic_alt]}"
        )
    # An insertion that isn't a duplication stays explicit
    hgvs_g = synthetic_babelfish.vcf_to_g_hgvs(
        "2", position, pad, pad + reverse_complement(preceding)
    )
    assert hgvs_g.posedit.edit.type == "ins"


def test_symbolic_size_normalized(parser, synthetic_babelfish):
    """Between symbolic_alt_min_length and normalize_max_length: symbolic, but shifted"""
    seq = synthetic_babelfish.hdp.seq
    start, middle, end = _SyntheticHDP.TANDEM_REPEAT
    pad = seq[start - 1]
    repeat = seq[start:middle]
    for typ, symbolic_alt, ref, alt in [
        ("del", "<DEL>", pad + repeat, pad),
        ("dup", "<DUP>", pad, pad + repeat),
    ]:
        hgvs_string = f"NC_000002.12:g.{middle + 1}_{end}{typ}"  # 3' shifted
        explicit_g = synthetic_babelfish.vcf_to_g_hgvs("2", start, ref, alt)
        symbolic_g = synthetic_babelfish.vcf_to_g_hgvs("2", start, pad, symbolic_alt, end=middle)
        assert explicit_g.format() == symbolic_g.format() == hgvs_string
        # Left aligned
        for hgvs_input in [hgvs_string, f"NC_000002.12:g.{start + 1}_{middle}{typ}"]:
            variant_coordinate = synthetic_babelfish.hgvs_to_vcf(parser.parse(hgvs_input))
            assert variant_coordinate == ("2", start, pad, symbolic_alt, typ)


def test_vcf_version_svlen(parser, hdp):
    """A deletion's SVLEN is negative up to VCF 4.3 and positive from 4.4; either is read"""
    var_g = parser.parse("NC_000002.12:g.1000000_1000999del")
    for vcf_version, svlen in (("4.3", -1000), ("4.4", 1000)):
        babelfish = _anchor_babelfish(hdp, 1000, 1000, vcf_version=vcf_version)
        vc = babelfish.hgvs_to_vcf_coordinate(var_g)
        assert vc == VCFCoordinate(CHR2, 999999, "T", "<DEL>", end=1000999)
        assert vc.svlen == 1000  # the coordinate doesn't know the VCF version
        assert babelfish.vcf_info(vc)["SVLEN"] == svlen
        assert babelfish.vcf_to_g_hgvs("2", 999999, "T", "<DEL>", svlen=svlen) == var_g
        dup = VCFCoordinate(CHR2, 999999, "T", "<DUP>", end=1000999)
        assert babelfish.vcf_info(dup)["SVLEN"] == 1000
    with pytest.raises(ValueError, match="vcf_version"):
        Babelfish(hdp, "GRCh38", vcf_version="4")


def test_vcf_coordinate_explicit_end(babelfish38):
    vc = VCFCoordinate("NC_000006.12", 49949409, "GAA", "G")
    assert vc.end == 49949411
    assert vc.svlen is None
    assert babelfish38.vcf_info(vc) == {}
    assert str(vc) == "NC_000006.12:49949409 GAA>G"
    assert babelfish38.vcf_coordinate("6", 49949409, "GAA", "G") == vc
    with pytest.raises(ValueError, match="ends at 49949411"):
        VCFCoordinate("NC_000006.12", 49949409, "GAA", "G", end=49949412)


def test_vcf_coordinate_as_symbolic_and_explicit(synthetic_babelfish):
    seq = synthetic_babelfish.hdp.seq
    position, end = 500, 560  # 60bp >= symbolic_alt_min_length
    pad = seq[position - 1]
    span = seq[position:end]
    explicit_symbolic = [
        (
            VCFCoordinate(CHR2, position, pad + span, pad),
            VCFCoordinate(CHR2, position, pad, "<DEL>", end),
        ),
        (
            VCFCoordinate(CHR2, position, pad, pad + span),
            VCFCoordinate(CHR2, position, pad, "<DUP>", end),
        ),
        (
            VCFCoordinate(CHR2, position + 1, span, reverse_complement(span)),
            VCFCoordinate(CHR2, position, pad, "<INV>", end),
        ),
    ]
    for explicit, symbolic in explicit_symbolic:
        assert synthetic_babelfish.vcf_coordinate_as_symbolic(explicit) == symbolic
        assert synthetic_babelfish.vcf_coordinate_as_explicit(symbolic) == explicit
        # Already in the requested form
        assert synthetic_babelfish.vcf_coordinate_as_symbolic(symbolic) == symbolic
        assert synthetic_babelfish.vcf_coordinate_as_explicit(explicit) == explicit
    # Copy after the original - the padding base moves to before the original
    copy_after = VCFCoordinate(CHR2, end, seq[end - 1], seq[end - 1] + span)
    assert synthetic_babelfish.vcf_coordinate_as_symbolic(copy_after) == VCFCoordinate(
        CHR2, position, pad, "<DUP>", end
    )
    # A multi-base ref isn't a dup, however long the insertion
    delins = VCFCoordinate(CHR2, position, seq[position - 1 : end], pad + "GATTACA" * 20)
    assert synthetic_babelfish.vcf_coordinate_as_symbolic(delins) == delins
    # Without symbolic_alt_min_length, a min_length is needed
    babelfish = Babelfish(synthetic_babelfish.hdp, assembly_name="GRCh38")
    with pytest.raises(ValueError, match="min_length"):
        babelfish.vcf_coordinate_as_symbolic(delins)
    assert babelfish.vcf_coordinate_as_symbolic(explicit_symbolic[0][0], min_length=50).is_symbolic


def test_hgvs_to_vcf_small_inv_explicit(parser, synthetic_babelfish):
    """An inv below symbolic_alt_min_length is explicit ref/alt, with no padding base"""
    seq = synthetic_babelfish.hdp.seq
    start, end = 101, 120
    span = seq[start - 1 : end]
    var_g = parser.parse(f"NC_000002.12:g.{start}_{end}inv")
    expected = ("2", start, span, reverse_complement(span), "inv")
    assert synthetic_babelfish.hgvs_to_vcf(var_g) == expected
    assert synthetic_babelfish.vcf_to_g_hgvs(*expected[:4]).format() == var_g.format()


def test_symbolic_alt_min_length_boundary(parser, synthetic_babelfish):
    """Exactly symbolic_alt_min_length (50) bases is symbolic, one fewer is explicit"""
    seq = synthetic_babelfish.hdp.seq
    position = 500
    pad = seq[position - 1]
    for length, is_symbolic in ((50, True), (49, False)):
        span = seq[position : position + length]
        explicit = [
            VCFCoordinate(CHR2, position, pad + span, pad),
            VCFCoordinate(CHR2, position, pad, pad + span),
            VCFCoordinate(CHR2, position + 1, span, reverse_complement(span)),
        ]
        for vc in explicit:
            assert synthetic_babelfish.vcf_coordinate_as_symbolic(vc).is_symbolic == is_symbolic
        for typ in SYMBOLIC_ALTS:
            var_g = parser.parse(f"NC_000002.12:g.{position + 1}_{position + length}{typ}")
            assert synthetic_babelfish.hgvs_to_vcf_coordinate(var_g).is_symbolic == is_symbolic


def test_large_delins_not_symbolic(parser, hdp):
    """ClinVar 869248: a 1517bp delins - its alt isn't the first ref base, so it isn't a <DEL>"""
    babelfish = Babelfish(hdp, assembly_name="GRCh37", symbolic_alt_min_length=50)
    start, end = 5247125, 5248641
    ref = hdp.get_seq("NC_000011.9", start - 1, end)
    hgvs_string = f"NC_000011.9:g.{start}_{end}delinsT"
    assert babelfish.vcf_to_g_hgvs("11", start, ref, "T").format() == hgvs_string
    vc = babelfish.hgvs_to_vcf_coordinate(parser.parse(hgvs_string))
    assert vc == VCFCoordinate("NC_000011.9", start, ref, "T")


def test_common_prefix_len_partial_full_empty():
    """Cover partial, full, and empty common prefixes for sequence strings."""
    # Partial
    assert _common_prefix_len("ACGT", "ACTT") == 2
    # Full (one is prefix of the other, or equal)
    assert _common_prefix_len("ACGT", "ACGT") == 4
    assert _common_prefix_len("AC", "ACGT") == 2
    # Empty
    assert _common_prefix_len("ACGT", "TGCA") == 0
    assert _common_prefix_len("", "ACGT") == 0
    assert _common_prefix_len("A", "") == 0
