import os
from contextlib import contextmanager

import hgvs
from hgvs.exceptions import HGVSError
import pytest
import support.mock_input_source as mock_input_data_source

sanity_cases = [
    {
        "name": "ins at intron/exon boundary",
        "var_c": "NM_999999.1:c.2-1_2insG",
        "exonic":{
            "var_p": "MOCK:p.Met1?",
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },
    },
    {
        "name": "ins at exon/intron boundary",
        "var_c": "NM_999999.1:c.26_26+1insC",
        "exonic":{
            "var_p": "MOCK:p.(Lys9AsnfsTer?)",
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary",
        "var_c": "NM_999999.1:c.2dup",
        "exonic":{
            "var_p": "MOCK:p.Met1?",
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary",
        "var_c": "NM_999999.1:c.2-1dup",
        "exonic":{
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },

    },
    {
        "name": "dup at exon/intron boundary",
        "var_c": "NM_999999.1:c.26dup",
        "exonic":{
            "var_p": "MOCK:p.(Ter10IleextTer?)",
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary",
        "var_c": "NM_999999.1:c.26+1dup",
        "exonic":{
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "MOCK:p.?",
        },
    },
]

real_cases = [
    {
        "name": "ins at intron/exon boundary",
        "var_c": "NM_004380.2:c.3251-1_3251insA",
        "exonic":{
            "var_p": "NP_004371.2:p.(Ile1084AsnfsTer3)",
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "ins at exon/intron boundary",
        "var_c": "NM_004380.2:c.3250_3250+1insT",
        "exonic":{
            "var_p": "NP_004371.2:p.(Phe1085LeufsTer2)",
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in exon, single-base)",
        "var_c": "NM_024529.4:c.132dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45AspfsTer21)",
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in exon, multi-base)",
        "var_c": "NM_024529.4:c.132_133dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45ArgfsTer65)",
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in intron, single-base)",
        "var_c": "NM_024529.4:c.132-1dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45AspfsTer21)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in intron, multi-base)",
        "var_c": "NM_024529.4:c.132-2_132-1dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45GlyfsTer65)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in exon, single-base, negative strand)",
        "var_c": "NM_004380.2:c.3251dup",
        "exonic":{
            "var_p": "NP_004371.2:p.(Phe1085LeufsTer2)",
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in exon, multi-base, negative strand)",
        "var_c": "NM_004380.2:c.3251_3252dup",
        "exonic":{
            "var_p": "NP_004371.2:p.(Phe1085SerfsTer15)",
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in intron, single-base, negative strand)",
        "var_c": "NM_004380.2:c.3251-1dup",
        "exonic":{
            "var_p": "NP_004371.2:p.(Ile1084SerfsTer3)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "dup at intron/exon boundary (in intron, multi-base, negative strand)",
        "var_c": "NM_004380.2:c.3251-2_3251-1dup",
        "exonic":{
            "var_p": "NP_004371.2:p.(Ile1084LysfsTer16)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_004371.2:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in exon, single-base)",
        "var_c": "NM_024529.4:c.131dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45AspfsTer21)",
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in exon, multi-base)",
        "var_c": "NM_024529.4:c.130_131dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45GlyfsTer65)",
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in intron, single-base)",
        "var_c": "NM_024529.4:c.131+1dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45AspfsTer21)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in intron, multi-base)",
        "var_c": "NM_024529.4:c.131+1_131+3dup",
        "exonic":{
            "var_p": "NP_078805.3:p.(Thr45Ter)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_078805.3:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in exon, single-base, negative strand)",
        "var_c": "NM_004985.4:c.111dup",
        "exonic":{
            "var_p": "NP_004976.2:p.(Asp38GlyfsTer10)",
        },
        "intronic":{
            "var_p": "NP_004976.2:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in exon, multi-base, negative strand)",
        "var_c": "NM_004985.4:c.110_111dup",
        "exonic":{
            "var_p": "NP_004976.2:p.(Asp38ArgfsTer8)",
        },
        "intronic":{
            "var_p": "NP_004976.2:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in intron, single-base, negative strand)",
        "var_c": "NM_004985.4:c.111+1dup",
        "exonic":{
            "var_p": "NP_004976.2:p.(Asp38GlyfsTer10)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_004976.2:p.?",
        },
    },
    {
        "name": "dup at exon/intron boundary (in intron, multi-base, negative strand)",
        "var_c": "NM_004985.4:c.111+1_111+4dup",
        "exonic":{
            "var_p": "NP_004976.2:p.(Asp38ValfsTer11)",
            "exc": HGVSError,
        },
        "intronic":{
            "var_p": "NP_004976.2:p.?",
        },
    },
    {
        "name": "snv with broken cigar mapping",
        "var_c": "NM_014638.2:c.515-99G>T",
        "exonic":{
            "var_p": "NP_055453.2:p.?",
        },
        "intronic":{
            "var_p": "NP_055453.2:p.?",
        },
    },
    {
        "name": "ins with broken cigar mapping",
        "var_c": "NM_004799.2:c.71-7048_71-7047insATAT",
        "exonic":{
            "var_p": "NP_004790.2:p.?",
        },
        "intronic":{
            "var_p": "NP_004790.2:p.?",
        },
    },
    {
        "name": "del with broken cigar mapping",
        "var_c": "NM_007324.2:c.71-7826_71-7802del",
        "exonic":{
            "var_p": "NP_015563.2:p.?",
        },
        "intronic":{
            "var_p": "NP_015563.2:p.?",
        },
    },
]


@pytest.fixture(scope="module")
def mock_hdp():
    fn = os.path.join(os.path.dirname(__file__), "..", "data", "sanity_cp.tsv")
    return mock_input_data_source.MockInputSource(fn)


@pytest.fixture(scope="module")
def mock_vm(mock_hdp):
    return hgvs.variantmapper.VariantMapper(mock_hdp, prevalidation_level="INTRINSIC")


@contextmanager
def config_setting(module, setting, temp_value):
    old_value = getattr(module, setting)
    try:
        setattr(module, setting, temp_value)
        yield
    finally:
        setattr(module, setting, old_value)


@pytest.mark.parametrize("case", sanity_cases)
def test_sanity_c_to_p(case, parser, mock_vm):
    with config_setting(hgvs.global_config.mapping, 'ins_at_boundary_is_intronic', False):
        var_c = parser.parse(case["var_c"])
        if "exc" in case["exonic"]:
            with pytest.raises(case["exonic"]["exc"]):
                mock_vm.c_to_p(var_c, "MOCK")
        else:
            assert str(mock_vm.c_to_p(var_c, "MOCK")) == case["exonic"]["var_p"]
    with config_setting(hgvs.global_config.mapping, 'ins_at_boundary_is_intronic', True):
        var_c = parser.parse(case["var_c"])
        if "exc" in case["intronic"]:
            with pytest.raises(case["intronic"]["exc"]):
                mock_vm.c_to_p(var_c, "MOCK")
        else:
            assert str(mock_vm.c_to_p(var_c, "MOCK")) == case["intronic"]["var_p"]


@pytest.mark.parametrize("case", real_cases)
def test_real_c_to_p(case, parser, vm, am37):
    with config_setting(hgvs.global_config.mapping, 'ins_at_boundary_is_intronic', False):
        var_c = parser.parse(case["var_c"])
        if "exc" in case["exonic"]:
            with pytest.raises(case["exonic"]["exc"]):
                vm.c_to_p(var_c)
        else:
            assert str(vm.c_to_p(var_c)) == case["exonic"]["var_p"]
        assert str(am37.c_to_p(var_c)) == case["exonic"]["var_p"]
    with config_setting(hgvs.global_config.mapping, 'ins_at_boundary_is_intronic', True):
        var_c = parser.parse(case["var_c"])
        if "exc" in case["intronic"]:
            with pytest.raises(case["intronic"]["exc"]):
                vm.c_to_p(var_c)
        else:
            assert str(vm.c_to_p(var_c)) == case["intronic"]["var_p"]
        assert str(am37.c_to_p(var_c)) == case["intronic"]["var_p"]
