import pytest


NORM_HGVS_VCF = [
    # Columns are: (normed-HGVS, non-normalized HGVS, VCF coordinates, non-norm VCF)
    # no-op
    (
        "NC_000006.12:g.49949407=",
        [],
        ("6", 49949407, "A", ".", "identity"),
        [("6", 49949407, "A", "A", "identity"),
         # Test case insensitivity
         ("6", 49949407, "A", "a", "identity"),
         ("6", 49949407, "a", "A", "identity"),]
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
        for hgvs_string in [norm_hgvs_string] + alt_hgvs:
            variant_coordinates = _h2v(hgvs_string)
            assert variant_coordinates == expected_variant_coordinate


def test_vcf_to_hgvs(parser, babelfish38):
    def _v2h(*v):
        return babelfish38.vcf_to_g_hgvs(*v)

    for expected_hgvs_string, _, norm_variant_coordinate, alt_variant_coordinate in NORM_HGVS_VCF:
        for variant_coordinate in [norm_variant_coordinate] + alt_variant_coordinate:
            *v, typ = variant_coordinate  # last column is type ie "dup"
            hgvs_g = _v2h(*v)
            hgvs_string = hgvs_g.format()
            assert hgvs_string == expected_hgvs_string


def test_vcf_to_hgvs_contig_chrom(parser, babelfish38):
    hgvs_g = babelfish38.vcf_to_g_hgvs("NC_000006.12", 49949409, "GAA", "G")
    assert hgvs_g.format() == "NC_000006.12:g.49949413_49949414del"
