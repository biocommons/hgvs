import pytest

# https://github.com/biocommons/hgvs/issues/801
# Uncertain c./n. intervals should parse the same combinations that were
# already supported for g. intervals.

testdata = [
    "NM_000333.3:c.(4_246)del",
    "NM_000333.3:c.(4_6)_246del",
    "NM_000333.3:c.4_(245_246)del",
    "NM_000333.3:c.(4_6)_(245_246)del",
    "NM_000333.3:c.(?_6)_245del",
    "NM_000333.3:c.4_(245_?)del",
    "NM_000333.3:c.(?_6)_(245_?)del",
]


@pytest.mark.parametrize("hgvs_str", testdata)
def test_parse_uncertain_c_interval(hgvs_str, parser):
    var = parser.parse(hgvs_str)
    assert str(var) == hgvs_str


@pytest.mark.parametrize(
    "hgvs_str",
    [s.replace("NM_000333.3:c.", "NM_000333.3:n.") for s in testdata],
)
def test_parse_uncertain_n_interval(hgvs_str, parser):
    var = parser.parse(hgvs_str)
    assert str(var) == hgvs_str
