import pytest

import hgvs.easy
from hgvs.extras.babelfish import Babelfish


@pytest.fixture(scope="session", autouse=True)
def parser():
    return hgvs.easy.parser


@pytest.fixture(scope="session")
def am37():
    return hgvs.easy.am37


@pytest.fixture(scope="session")
def am38():
    return hgvs.easy.am38


@pytest.fixture(scope="session")
def hdp():
    return hgvs.easy.hdp


@pytest.fixture(scope="session")
def babelfish38(hdp):
    return Babelfish(hdp, assembly_name="GRCh38")



@pytest.fixture(scope="class")
def kitchen_sink_setup(request, hdp, parser, am37, am38):
    """
    Adds common fixtures to an instance passed via request.cls. 

    For example:

    @pytst.mark.usefixtures("kitchen_sink_setup")
    class MyTest(unittest.TestCase):
        def test_foo(self):
             self.parser.parse(...)

    The class must be subclass of unittest.TestCase

    See https://docs.pytest.org/en/5.3.5/unittest.html#mixing-pytest-fixtures-into-unittest-testcase-subclasses-using-marks
    
    """

    request.cls.hdp = hdp
    request.cls.parser = parser
    request.cls.am37 = am37
    request.cls.am38 = am38
