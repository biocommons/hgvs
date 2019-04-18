import pytest

import hgvs.easy


@pytest.fixture(scope="session")
def parser():
    return hgvs.easy.parser


@pytest.fixture(scope="session")
def am38():
    return hgvs.easy.am38
