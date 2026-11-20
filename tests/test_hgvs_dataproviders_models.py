"""Tests for the data provider models"""

import pickle

import pytest

from hgvs.dataproviders.models import TranscriptInfo


@pytest.fixture
def tx_info():
    return TranscriptInfo(
        tx_ac="NM_000051.3",
        alt_ac="AC_000143.1",
        alt_aln_method="splign",
        cds_start_i=385,
        cds_end_i=9556,
        hgnc="ATM",
    )


def test_transcript_info_is_frozen(tx_info):
    with pytest.raises(AttributeError):
        tx_info.cds_start_i = 0


def test_transcript_info_round_trips_through_pickle(tx_info):
    assert pickle.loads(pickle.dumps(tx_info)) == tx_info  # noqa: S301 - trusted, built above


def test_transcript_info_has_no_dict_access(tx_info):
    with pytest.raises(TypeError):
        tx_info["cds_start_i"]
