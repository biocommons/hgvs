import os

import hgvs
import pytest
from support import CACHE

cases = [
    {
        "name": "ins with splice region preserved",
        "var_c": "NM_004119.2:c.1837+21_1837+22insCGAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAA",
        "var_p": "NP_004110.2:p.(Lys614_Val615insAsnGlyMetCysGlnThrArgGluTyrGluTyrAspLeuLysTrpGluPheProArgGluAsnLeuGluPheGlyLys)"
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_004119.2:c.1835_1837+3dup",
        "var_p": "NP_004110.2:p.(Gly613_Lys614insIleGly)"
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_005228.4:c.2284-5_2290dup",
        "var_p": "NP_005219.2:p.(Ala763_Tyr764insPheGlnGluAla)"
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_004456.4:c.2196-1_2196dup",
        "var_p": "NP_004447.2:p.(Tyr733AspfsTer8)"
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_016222.3:c.27+2_27+5dup",
        "var_p": "NP_057306.2:p.(Arg10ValfsTer20)"
    },
    {
        "name": "dup with splice region preserved",
        "var_c": "NM_182758.2:c.2953-31_2953-26dup",
        "var_p": "NP_877435.2:p.?"
    },
    {
        "name": "dup with broken cigar mapping",
        "var_c": "NM_000267.3:c.8315-290_8457dup",
        "var_p": "NP_000258.1:p.?"
    }
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
def am37(hdp):
    return hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")


@pytest.mark.parametrize("case", cases)
def test_real_c_to_p(case, hp, am37):
    var_c = hp.parse(case["var_c"])
    assert str(am37.c_to_p(var_c)) == case["var_p"]
