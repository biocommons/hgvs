from os import path
from unittest.mock import patch

from biocommons.seqrepo.dataproxy import SeqRepoRESTDataProxy
from biocommons.seqrepo.seqrepo import SeqRepo
from bioutils import seqfetcher

from hgvs.dataproviders.seqfetcher import SeqFetcher


@patch.object(path, "exists")
def test_seqfetcher_initialized_with_seqrepo_dir(mock_path_exists, monkeypatch):
    mock_path_exists.return_value = True
    monkeypatch.setenv("HGVS_SEQREPO_DIR", "/tmp/my-seqrepo-location")
    with patch("sqlite3.connect") as mock_sqlite_connect:
        sf = SeqFetcher()
    mock_path_exists.assert_called_once_with("/tmp/my-seqrepo-location")
    assert mock_sqlite_connect.call_args[0][0] == "/tmp/my-seqrepo-location/aliases.sqlite3"
    assert isinstance(sf.sr, SeqRepo)
    assert sf.source == "SeqRepo (/tmp/my-seqrepo-location)"


def test_seqfetcher_initialized_with_seqrepo_url(monkeypatch):
    monkeypatch.setenv("HGVS_SEQREPO_URL", "http://localhost:5000/seqrepo")
    sf = SeqFetcher()
    assert isinstance(sf.sr, SeqRepoRESTDataProxy)
    # the URL should be versioned automatically for us
    assert sf.sr.base_url == "http://localhost:5000/seqrepo/1/"
    assert sf.source == "SeqRepo REST (http://localhost:5000/seqrepo)"


def test_seqfetcher_initialized_with_public_seqrepo_sources(monkeypatch):
    monkeypatch.setenv("HGVS_SEQREPO_URL", "")
    sf = SeqFetcher()
    assert sf.sr is None
    assert sf.fetcher == seqfetcher.fetch_seq
    assert sf.source == "bioutils.seqfetcher (network fetching)"
