from unittest.mock import MagicMock

import psycopg2
import pytest

import hgvs.dataproviders.ncbi
from hgvs.exceptions import HGVSError


def _make_ncbi_stub(pooling):
    """NCBI_postgresql instance wired with mocks, without connecting."""
    hdp = object.__new__(hgvs.dataproviders.ncbi.NCBI_postgresql)
    hdp.pooling = pooling
    hdp.url = MagicMock()
    hdp.url.schema = "ncbi"
    hdp._connect = MagicMock()

    cur = MagicMock(name="cursor")
    conn = MagicMock(name="connection")
    conn.cursor.return_value = cur

    if pooling:
        hdp._pool = MagicMock(name="pool")
        hdp._pool.getconn.return_value = conn
    else:
        hdp._conn = conn

    return hdp, conn, cur


class TestNCBIGetCursorContextManager:
    def test_success_yields_cursor_and_cleans_up_pooling(self):
        hdp, conn, cur = _make_ncbi_stub(pooling=True)

        with hdp._get_cursor() as c:
            assert c is cur

        cur.close.assert_called_once_with()
        hdp._pool.putconn.assert_called_once_with(conn)

    def test_success_closes_cursor_without_pooling(self):
        hdp, _conn, cur = _make_ncbi_stub(pooling=False)

        with hdp._get_cursor() as c:
            assert c is cur

        cur.close.assert_called_once_with()

    def test_consumer_exception_propagates_and_cleans_up(self):
        hdp, conn, cur = _make_ncbi_stub(pooling=True)

        with pytest.raises(ValueError, match="boom"), hdp._get_cursor():
            raise ValueError("boom")

        cur.close.assert_called_once_with()
        hdp._pool.putconn.assert_called_once_with(conn)

    def test_acquisition_operational_error_retries_then_raises(self):
        """OperationalError while acquiring the cursor triggers reconnect
        and, once retries are exhausted, raises HGVSError."""
        hdp, conn, _cur = _make_ncbi_stub(pooling=False)
        conn.cursor.side_effect = psycopg2.OperationalError("lost")

        with pytest.raises(HGVSError), hdp._get_cursor(n_retries=1):
            pass

        assert hdp._connect.call_count == 2

    def test_consumer_operational_error_does_not_try_to_loop_in_context_manager(self):
        """An OperationalError raised *inside* the ``with`` block must
        propagate rather than being caught by the acquisition retry logic
        (which previously re-entered the generator and raised
        "generator didn't stop")."""
        hdp, _conn, cur = _make_ncbi_stub(pooling=False)

        with pytest.raises(psycopg2.OperationalError), hdp._get_cursor():
            raise psycopg2.OperationalError("query failed")

        hdp._connect.assert_not_called()
        cur.close.assert_called_once_with()

# <LICENSE>
# Copyright 2026 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
