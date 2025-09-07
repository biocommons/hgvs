import os
from pathlib import Path

import pytest
from support import CACHE

from hgvs.assemblymapper import AssemblyMapper
import hgvs.easy
from hgvs.extras.babelfish import Babelfish
from hgvs.variantmapper import VariantMapper


def remove_request_headers(request):
    headers_to_remove = ["user-agent"]
    headers = getattr(request, "headers", request)
    for header in headers_to_remove:
        headers.pop(header, None)
        headers.pop(header.title(), None)
    return request


def remove_response_headers(response):
    headers_to_remove = ["date", "server"]
    for header in headers_to_remove:
        response["headers"].pop(header, None)
        response["headers"].pop(header.title(), None)
    return response


@pytest.fixture(scope="function")
def vcr_config(request):
    """See https://pytest-vcr.readthedocs.io/en/latest/configuration/"""
    test_file_path = Path(request.node.fspath)
    return {
        "cassette_library_dir": str(
            test_file_path.with_name("cassettes") / test_file_path.stem
        ),
        "record_mode": os.environ.get("VCR_RECORD_MODE", "new_episodes"),
        "cassette_name": f"{request.node.name}.yaml",
        "before_record_request": remove_request_headers,
        "before_record_response": remove_response_headers,
    }


hgvs.easy.hdp = hgvs.dataproviders.uta.connect(
    mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
)


@pytest.fixture(scope="session", autouse=True)
def parser():
    return hgvs.easy.parser


@pytest.fixture(scope="session")
def vm(hdp):
    return VariantMapper(hdp)


@pytest.fixture(scope="session")
def am37(hdp):
    return AssemblyMapper(hdp, assembly_name="GRCh37")


@pytest.fixture(scope="session")
def am38(hdp):
    return AssemblyMapper(hdp, assembly_name="GRCh38")


@pytest.fixture(scope="session")
def hdp():
    return hgvs.easy.hdp


@pytest.fixture(scope="session")
def babelfish38(hdp):
    return Babelfish(hdp, assembly_name="GRCh38")


def pytest_report_header(config):
    env_vars = ["UTA_DB_URL", "HGVS_SEQREPO_URL", "HGVS_CACHE_MODE"]
    rv = [f"{ev}: {os.environ.get(ev)}" for ev in sorted(env_vars)]
    rv += [f"hgvs.easy.hdp={hgvs.easy.hdp.url}"]
    rv += [f"{hgvs.easy.hdp.seqfetcher.source=}"]
    return "\n".join(rv)


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
