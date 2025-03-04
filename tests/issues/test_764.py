import os

import hgvs
import pytest
from support import CACHE


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


def test_vm_g_to_c(hp, vm):
    var_g = hp.parse("NC_000007.13:g.106545832_106545833insAG")
    var_c = vm.g_to_c(var_g, "NM_002649.2")
    assert str(var_c) == "NM_002649.2:c.3309_*1insAG"


def test_am37_g_to_c(hp, am37):
    var_g = hp.parse("NC_000007.13:g.106545832_106545833insAG")
    var_c = am37.g_to_c(var_g, "NM_002649.2")
    assert str(var_c) == "NM_002649.2:c.3309_*1insAG"
