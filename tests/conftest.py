import pytest

import hgvs.easy
from hgvs.extras.babelfish import Babelfish


@pytest.fixture(scope="session")
def parser():
    return hgvs.easy.parser


@pytest.fixture(scope="session")
def am38():
    return hgvs.easy.am38


@pytest.fixture(scope="session")
def hdp():
    return hgvs.easy.hdp


@pytest.fixture(scope="session")
def babelfish38(hdp):
    return Babelfish(hdp, assembly_name="GRCh38")
