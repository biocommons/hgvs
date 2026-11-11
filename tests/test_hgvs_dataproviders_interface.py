"""Tests for the abstract data provider interface and the version 4 bridge"""

import pytest

from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.models import TranscriptInfo

TX_INFO_DICT = {
    "hgnc": "ATM",
    "cds_start_i": 385,
    "cds_end_i": 9556,
    "tx_ac": "NM_000051.3",
    "alt_ac": "AC_000143.1",
    "alt_aln_method": "splign",
}
TX_INFO = TranscriptInfo(**TX_INFO_DICT)
TX_INFO_ARGS = ("NM_000051.3", "AC_000143.1", "splign")


class V4Provider(Interface):
    """A provider written against interface version 4: get_tx_info returns a dict."""

    required_version = "1.1"
    url = "fake://v4"

    def __init__(self, tx_info_result):
        self._tx_info_result = tx_info_result
        self.calls = []
        super().__init__()

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self.calls.append(("get_tx_info", tx_ac, alt_ac, alt_aln_method))
        return self._tx_info_result

    def data_version(self):
        return "fake"

    def schema_version(self):
        return "1.1"

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        self.calls.append(("get_tx_exons", tx_ac, alt_ac, alt_aln_method))
        return [{"tx_ac": tx_ac, "alt_ac": alt_ac, "alt_aln_method": alt_aln_method}]

    # the remaining abstract methods are not exercised here
    def get_acs_for_protein_seq(self, seq):
        raise NotImplementedError

    def get_assembly_map(self, assembly_name):
        raise NotImplementedError

    def get_gene_info(self, gene):
        raise NotImplementedError

    def get_pro_ac_for_tx_ac(self, tx_ac):
        raise NotImplementedError

    def get_seq(self, ac, start_i=None, end_i=None):
        raise NotImplementedError

    def get_similar_transcripts(self, tx_ac):
        raise NotImplementedError

    def get_tx_for_gene(self, gene):
        raise NotImplementedError

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError

    def get_tx_identity_info(self, tx_ac):
        raise NotImplementedError

    def get_tx_mapping_options(self, tx_ac):
        raise NotImplementedError


def test_bridge_is_a_version_5_interface():
    hdp = Interface.from_v4(V4Provider(dict(TX_INFO_DICT)))
    assert isinstance(hdp, Interface)
    assert hdp.interface_version() == 5


def test_bridge_converts_dict_tx_info_to_model():
    provider = V4Provider(dict(TX_INFO_DICT))
    hdp = Interface.from_v4(provider)
    assert hdp.get_tx_info(*TX_INFO_ARGS) == TX_INFO
    assert provider.calls == [("get_tx_info", *TX_INFO_ARGS)]


def test_bridge_passes_none_tx_info_through():
    hdp = Interface.from_v4(V4Provider(None))
    assert hdp.get_tx_info(*TX_INFO_ARGS) is None


def test_bridge_rejects_dict_with_unknown_keys():
    hdp = Interface.from_v4(V4Provider({**TX_INFO_DICT, "bogus": 1}))
    with pytest.raises(TypeError, match="bogus"):
        hdp.get_tx_info(*TX_INFO_ARGS)


def test_bridge_delegates_unchanged_methods_as_is():
    provider = V4Provider(dict(TX_INFO_DICT))
    hdp = Interface.from_v4(provider)
    assert hdp.get_tx_exons(*TX_INFO_ARGS) == provider.get_tx_exons(*TX_INFO_ARGS)
    assert hdp.schema_version() == "1.1"
    assert hdp.data_version() == "fake"


def test_bridge_exposes_provider_attributes():
    provider = V4Provider(dict(TX_INFO_DICT))
    hdp = Interface.from_v4(provider)
    assert hdp.url == "fake://v4"
    assert hdp.required_version == "1.1"
    with pytest.raises(AttributeError):
        _ = hdp.no_such_attribute
