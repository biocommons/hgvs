from contextlib import contextmanager

import hgvs
import pytest

cases = [
    {
        "name": "ins with splice region preserved",
        "var_c": "NM_004119.2:c.1837+21_1837+22insCGAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAA",
        "exonic":{
            "var_p": "NP_004110.2:p.(Lys614_Val615insAsnGlyMetCysGlnThrArgGluTyrGluTyrAspLeuLysTrpGluPheProArgGluAsnLeuGluPheGlyLys)",
            "is_shifted": True,
        },
        "intronic": {
            "var_p": "NP_004110.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_004119.2:c.1835_1837+3dup",
        "exonic":{
            "var_p": "NP_004110.2:p.(Gly613_Lys614insIleGly)",
            "is_shifted": True,
        },
        "intronic": {
            "var_p": "NP_004110.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_005228.4:c.2284-5_2290dup",
        "exonic":{
            "var_p": "NP_005219.2:p.(Ala763_Tyr764insPheGlnGluAla)",
            "is_shifted": True,
        },
        "intronic": {
            "var_p": "NP_005219.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_004456.4:c.2196-1_2196dup",
        "exonic":{
            "var_p": "NP_004447.2:p.(Tyr733AspfsTer8)",
            "is_shifted": True,
        },
        "intronic": {
            "var_p": "NP_004447.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_016222.3:c.27+2_27+5dup",
        "exonic":{
            "var_p": "NP_057306.2:p.(Arg10ValfsTer20)",
            "is_shifted": True,
        },
        "intronic": {
            "var_p": "NP_057306.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_182758.2:c.2953-31_2953-26dup",
        "exonic":{
            "var_p": "NP_877435.2:p.?",
            "is_shifted": False,
        },
        "intronic": {
            "var_p": "NP_877435.2:p.?",
            "is_shifted": False,
        },
    },
    {
        "name": "dup with broken cigar mapping",
        "var_c": "NM_000267.3:c.8315-290_8457dup",
        "exonic":{
            "var_p": "NP_000258.1:p.?",
            "is_shifted": False,
        },
        "intronic": {
            "var_p": "NP_000258.1:p.?",
            "is_shifted": False,
        },
    }
]


@pytest.mark.parametrize("case", cases)
def test_real_c_to_p(case, parser, am37, monkeypatch):
    monkeypatch.setattr(hgvs.global_config.mapping, 'shift_over_boundary_is_intronic', False)

    var_c = parser.parse(case["var_c"])
    var_p = am37.c_to_p(var_c)
    assert str(var_p) == case["exonic"]["var_p"]
    if var_p.posedit:
        assert var_p.posedit.is_shifted == case["exonic"]["is_shifted"]

    monkeypatch.setattr(hgvs.global_config.mapping, 'shift_over_boundary_is_intronic', True)

    var_c = parser.parse(case["var_c"])
    var_p = am37.c_to_p(var_c)
    assert str(var_p) == case["intronic"]["var_p"]
    if var_p.posedit:
        assert var_p.posedit.is_shifted == case["intronic"]["is_shifted"]
