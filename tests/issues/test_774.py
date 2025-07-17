import pytest

testdata = [
    ('NM_000203.5:c.1957_*3del', 'NP_000194.2:p.(Pro653CysfsTer?)'),
    ('NM_004985.5:c.564_*2del', 'NP_004976.2:p.(Met188IlefsTer8)'),
    ('NM_004333.6:c.2298_*4del', 'NP_004324.2:p.(His766GlnfsTer26)'),
    ('NM_153223.4:c.2955_*2del', 'NP_694955.2:p.(Ser985ArgfsTer15)'),
]

@pytest.mark.parametrize("var_c_str,var_p_str", testdata)
def test_c_to_p(var_c_str, var_p_str, parser, am37):
    var_c = parser.parse(var_c_str)
    var_p = am37.c_to_p(var_c)
    assert str(var_p) == var_p_str
