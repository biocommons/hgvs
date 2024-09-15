import os

import hgvs
from hgvs.exceptions import HGVSError
import pytest
from support import CACHE
import support.mock_input_source as mock_input_data_source

sanity_cases = [
    {
        "name": "ins at intron/exon boundary",
        "var_c": "NM_999999.1:c.9-1_9insG",
        "var_p": "MOCK:p.(Lys4GlufsTer?)"
    },
    {
        "name": "ins at exon/intron boundary",
        "var_c": "NM_999999.1:c.39_39+1insC",
        "var_p": "MOCK:p.(=)"
    },
    {
        "name": "dup at intron/exon boundary",
        "var_c": "NM_999999.1:c.9-1dup",
        "var_p": "MOCK:p.?",
        "exc": HGVSError

    },
    {
        "name": "dup at exon/intron boundary",
        "var_c": "NM_999999.1:c.39+1dup",
        "var_p": "MOCK:p.?",
        "exc": HGVSError
    },
]

real_cases = [
    {
        "name": "ins at intron/exon boundary",
        "var_c": "NM_004380.2:c.3251-1_3251insA",
        "var_p": "NP_004371.2:p.(Ile1084AsnfsTer3)"
    },
    {
        "name": "ins at exon/intron boundary",
        "var_c": "NM_004380.2:c.3250_3250+1insT",
        "var_p": "NP_004371.2:p.(Phe1085LeufsTer2)"
    },
    {
        "name": "dup at intron/exon boundary",
        "var_c": "NM_024529.4:c.132-2_132-1dup",
        "var_p": "NP_078805.3:p.(Thr45GlyfsTer65)",
        "exc": HGVSError
    },
    {
        "name": "dup at intron/exon boundary",
        "var_c": "NM_004985.4:c.112-8_112-1dup",
        "var_p": "NP_004976.2:p.(Asp38LeufsTer10)",
        "exc": HGVSError
    },
    {
        "name": "dup at intron/exon boundary on negative strand",
        "var_c": "NM_004380.2:c.3251-1dup",
        "var_p": "NP_004371.2:p.(Ile1084SerfsTer3)",
        "exc": HGVSError
    },
    {
        "name": "dup at exon/intron boundary",
        "var_c": "NM_024529.4:c.131+1_131+3dup",
        "var_p": "NP_078805.3:p.(Thr45Ter)",
        "exc": HGVSError
    },
    {
        "name": "dup at exon/intron boundary on negative strand",
        "var_c": "NM_004985.4:c.111+1_111+4dup",
        "var_p": "NP_004976.2:p.(Asp38ValfsTer11)",
        "exc": HGVSError
    },
    {
        "name": "snv with broken cigar mapping",
        "var_c": "NM_014638.2:c.515-99G>T",
        "var_p": "NP_055453.2:p.?"
    },
    {
        "name": "ins with broken cigar mapping",
        "var_c": "NM_004799.2:c.71-7048_71-7047insATAT",
        "var_p": "NP_004790.2:p.?"
    },
    {
        "name": "del with broken cigar mapping",
        "var_c": "NM_007324.2:c.71-7826_71-7802del",
        "var_p": "NP_015563.2:p.?"
    },
]


@pytest.fixture(scope="module")
def hp():
    return hgvs.parser.Parser()


@pytest.fixture(scope="module")
def hdp():
    return hgvs.dataproviders.uta.connect(
        mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
    )


@pytest.fixture(scope="module")
def vm(hdp):
    return hgvs.variantmapper.VariantMapper(hdp)


@pytest.fixture(scope="module")
def am37(hdp):
    return hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")


@pytest.fixture(scope="module")
def mock_hdp():
    fn = os.path.join(os.path.dirname(__file__), "..", "data", "sanity_cp.tsv")
    return mock_input_data_source.MockInputSource(fn)


@pytest.fixture(scope="module")
def mock_vm(mock_hdp):
    return hgvs.variantmapper.VariantMapper(mock_hdp, prevalidation_level="INTRINSIC")


@pytest.mark.parametrize("case", sanity_cases)
def test_sanity_c_to_p(case, hp, mock_vm):
    var_c = hp.parse(case["var_c"])
    if "exc" in case:
        with pytest.raises(case["exc"]):
            mock_vm.c_to_p(var_c, "MOCK")
    else:
        assert str(mock_vm.c_to_p(var_c, "MOCK")) == case["var_p"]


@pytest.mark.parametrize("case", real_cases)
def test_real_c_to_p(case, hp, vm, am37):
    var_c = hp.parse(case["var_c"])
    if "exc" in case:
        with pytest.raises(case["exc"]):
            vm.c_to_p(var_c)
    else:
        assert str(vm.c_to_p(var_c)) == case["var_p"]
    assert str(am37.c_to_p(var_c)) == case["var_p"]
