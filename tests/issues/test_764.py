def test_vm_g_to_c(parser, vm):
    var_g = parser.parse("NC_000007.13:g.106545832_106545833insAG")
    var_c = vm.g_to_c(var_g, "NM_002649.2")
    assert str(var_c) == "NM_002649.2:c.3309_*1insAG"


def test_am37_g_to_c(parser, am37):
    var_g = parser.parse("NC_000007.13:g.106545832_106545833insAG")
    var_c = am37.g_to_c(var_g, "NM_002649.2")
    assert str(var_c) == "NM_002649.2:c.3309_*1insAG"
