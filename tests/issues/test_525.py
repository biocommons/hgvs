"""https://github.com/biocommons/hgvs/issues/525"""


def test_525(parser, am38):
    """https://github.com/biocommons/hgvs/issues/525"""

    # simple test case
    hgvs = "NM_000551.3:c.3_4insTAG"    # insert stop in phase at AA 2
    var_c = parser.parse_hgvs_variant(hgvs)
    var_p = am38.c_to_p(var_c)
    assert str(var_p) == "NP_000542.1:p.(Pro2Ter)"

    # variant reported in issue
    hgvs = "NM_001256125.1:c.1015_1016insAGGGACTGGGCGGGGCCATGGTCT"
    var_c = parser.parse_hgvs_variant(hgvs)
    var_p = am38.c_to_p(var_c)
    assert str(var_p) == "NP_001243054.2:p.(Trp339Ter)"
