"""Provides access to the NCBI table that allows to access transcripts by NCBI gene ids

This file is specific to Invitae and is nearly a copy of uta.py.  It
is not tested (see pytest.ini for exclusion).

"""

import os
import re
import six
import hgvs
import logging
import inspect
import psycopg2
import psycopg2.extras
import psycopg2.pool
import contextlib

from six.moves.urllib import parse as urlparse
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError

_logger = logging.getLogger(__name__)


def _stage_from_version(version):
    """return "prd", "stg", or "dev" for the given version string.  A value is always returned"""
    if version:
        m = re.match(r"^(?P<xyz>\d+\.\d+\.\d+)(?P<extra>.*)", version)
        if m:
            return "stg" if m.group("extra") else "prd"
    return "dev"


def _get_ncbi_db_url():
    """returns NCBI DB URL based on environment variables and code version

    * if NCBI_DB_URL is set, use that
    * Otherwise, if _NCBI_URL_KEY is set, use that as the name of a
      config file entry and use the corresponding URL
    * Otherwise,

    """

    if "NCBI_DB_URL" in os.environ:
        return os.environ["NCBI_DB_URL"]

    if "_NCBI_URL_KEY" in os.environ:
        url_key = os.environ["_NCBI_URL_KEY"]
    else:
        sdlc = _stage_from_version(hgvs.__version__)
        url_key = "public_{sdlc}".format(sdlc=sdlc)
    return hgvs.global_config['NCBI'][url_key]


def connect(db_url=None,
            pooling=hgvs.global_config.uta.pooling,
            application_name=None,
            mode=None,
            cache=None):
    """Connect to a uta/ncbi database instance.

    :param db_url: URL for database connection
    :type db_url: string
    :param pooling: whether to use connection pooling (postgresql only)
    :type pooling: bool
    :param application_name: log application name in connection (useful for debugging; PostgreSQL only)
    :type application_name: str

    When called with an explicit db_url argument, that db_url is used for connecting.

    When called without an explicit argument, the function default is
    determined by the environment variable UTA_DB_URL if it exists, or
    hgvs.datainterface.uta.public_db_url otherwise.

    >>> hdp = connect()
    >>> hdp.schema_version()
    '1.1'

    The format of the db_url is driver://user:pass@host/database (the same
    as that used by SQLAlchemy).  Examples:

    A remote public postgresql database:
        postgresql://anonymous:anonymous@uta.biocommons.org/uta'

    A local postgresql database:
        postgresql://localhost/uta

    A local SQLite database:
      sqlite:////tmp/uta-0.0.6.db

    For postgresql db_urls, pooling=True causes connect to use a
    psycopg2.pool.ThreadedConnectionPool.
    """

    _logger.debug('connecting to ' + str(db_url) + '...')

    if db_url is None:
        db_url = _get_ncbi_db_url()

    url = _parse_url(db_url)
    if url.scheme == 'postgresql':
        conn = NCBI_postgresql(
            url=url, pooling=pooling, application_name=application_name, mode=mode, cache=cache)
    else:
        # fell through connection scheme cases
        raise RuntimeError("{url.scheme} in {url} is not currently supported".format(url=url))
    _logger.info('connected to ' + str(db_url) + '...')
    return conn


class NCBIBase(object):
    required_version = "1.1"

    _queries = {
        "gene_id_for_hgnc":
        """
            select distinct(gene_id)
            from assocacs
            where hgnc=?
            """,
        "gene_id_for_tx":
        """
            select gene_id
            from assocacs
            where tx_ac=?
            """,
        "tx_for_gene_id":
        """
            select tx_ac 
            from assocacs
            where gene_id=?
            """,
        "hgnc_for_gene_id":
        """
            select distinct(hgnc) 
            from assocacs
            where gene_id=?
            """,
        "gene_info_for_gene_id":
        """
            select gene_id, tax_id, hgnc, maploc, aliases, type, summary, descr, xrefs   
            from geneinfo
            where gene_id=?
            """,
        "gene_info_for_hgnc":
        """
            select gene_id, tax_id, hgnc, maploc, aliases, type, summary, descr, xrefs   
            from geneinfo
            where hgnc=?
            """,
        "all_transcripts":
        """
                select distinct(tx_ac) 
                from assocacs                 
            """
    }

    def __init__(self, url, mode=None, cache=None):
        self.url = url
        if mode != 'run':
            self._connect()

    def __str__(self):
        return (
            "{n} <data_version:{dv}; schema_version:{sv}; application_name={self.application_name};"
            " url={self.url}; sequences-from={sf}>").format(
                n=type(self).__name__,
                self=self,
                dv=self.data_version(),
                sv=self.schema_version(),
                sf=os.environ.get("HGVS_SEQREPO_DIR", "seqfetcher"))

    def _fetchone(self, sql, *args):
        with self._get_cursor() as cur:
            cur.execute(sql, *args)
            return cur.fetchone()

    def _fetchall(self, sql, *args):
        with self._get_cursor() as cur:
            cur.execute(sql, *args)
            return cur.fetchall()

    def _update(self, sql, *args):
        with self._get_cursor() as cur:
            cur.execute(sql, *args)
            return cur.fetchall()

    ############################################################################
    # Queries

    def data_version(self):
        return self.url.schema

    def schema_version(self):
        return self._fetchone("select * from meta where key = 'schema_version'")['value']

    def get_ncbi_gene_id_for_hgnc(self, hgnc):
        rows = self._fetchall(self._queries['gene_id_for_hgnc'], [hgnc])
        return [r['gene_id'] for r in rows]

    def get_ncbi_gene_id_for_tx(self, tx_ac):
        rows = self._fetchall(self._queries['gene_id_for_tx'], [tx_ac])
        return [r['gene_id'] for r in rows]

    def get_tx_for_ncbi_gene_id(self, gene_id):
        rows = self._fetchall(self._queries['tx_for_gene_id'], [gene_id])
        return [r['tx_ac'] for r in rows]

    def get_hgnc_for_ncbi_gene_id(self, gene_id):
        rows = self._fetchall(self._queries['hgnc_for_gene_id'], [gene_id])
        return [r['hgnc'] for r in rows]

    def get_gene_info_for_ncbi_gene_id(self, gene_id):
        rows = self._fetchall(self._queries['gene_info_for_gene_id'], [gene_id])
        return rows

    def get_gene_info_for_hgnc(self, hgnc):
        rows = self._fetchall(self._queries['gene_info_for_hgnc'], [hgnc])
        return rows

    def get_all_transcripts(self):
        rows = self._fetchall(self._queries['all_transcripts'], [])
        return [r['tx_ac'] for r in rows]

    def store_assocacs(self, hgnc, tx_ac, gene_id, pro_ac, origin):
        sql = """
                insert into assocacs (hgnc, tx_ac, gene_id, pro_ac, origin)
                values (%s,%s,%s,%s,%s)
                
            """
        self._update(sql, [hgnc, tx_ac, gene_id, pro_ac, origin])


class NCBI_postgresql(NCBIBase):
    def __init__(self,
                 url,
                 pooling=hgvs.global_config.uta.pooling,
                 application_name=None,
                 mode=None,
                 cache=None):
        if url.schema is None:
            raise Exception("No schema name provided in {url}".format(url=url))
        self.application_name = application_name
        self.pooling = pooling
        self._conn = None
        super(NCBI_postgresql, self).__init__(url, mode, cache)

    def __del__(self):
        self.close()

    def close(self):
        if self.pooling:
            _logger.warning("Closing pool; future mapping and validation will fail.")
            self._pool.closeall()
        else:
            _logger.warning("Closing connection; future mapping and validation will fail.")
            if self._conn is not None:
                self._conn.close()

    def _connect(self):
        if self.application_name is None:
            st = inspect.stack()
            self.application_name = os.path.basename(st[-1][1])
        conn_args = dict(
            host=self.url.hostname,
            port=self.url.port,
            database=self.url.database,
            user=self.url.username,
            password=self.url.password,
            application_name=self.application_name + "/" + hgvs.__version__,
        )
        if self.pooling:
            _logger.info("Using UTA ThreadedConnectionPool")
            self._pool = psycopg2.pool.ThreadedConnectionPool(
                hgvs.global_config.uta.pool_min, hgvs.global_config.uta.pool_max, **conn_args)
        else:
            self._conn = psycopg2.connect(**conn_args)
            self._conn.autocommit = True

        self._ensure_schema_exists()

        # remap sqlite's ? placeholders to psycopg2's %s
        self._queries = {k: v.replace('?', '%s') for k, v in six.iteritems(self._queries)}

    def _ensure_schema_exists(self):
        # N.B. On AWS RDS, information_schema.schemata always returns zero rows
        r = self._fetchone("select exists(SELECT 1 FROM pg_namespace WHERE nspname = %s)",
                           [self.url.schema])
        if r[0]:
            return
        raise HGVSDataNotAvailableError("specified schema ({}) does not exist (url={})".format(
            self.url.schema, self.url))

    @contextlib.contextmanager
    def _get_cursor(self, n_retries=1):
        """Returns a context manager for obtained from a single or pooled
        connection, and sets the PostgreSQL search_path to the schema
        specified in the connection URL.

        Although *connections* are threadsafe, *cursors* are bound to
        connections and are *not* threadsafe. Do not share cursors
        across threads.

        Use this funciton like this::

            with hdp._get_cursor() as cur:
                # your code

        Do not call this function outside a contextmanager.

        """

        n_tries_rem = n_retries + 1
        while n_tries_rem > 0:
            try:

                conn = self._pool.getconn() if self.pooling else self._conn

                # autocommit=True obviates closing explicitly
                conn.autocommit = True

                cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
                cur.execute("set search_path = {self.url.schema};".format(self=self))

                yield cur

                # contextmanager executes these when context exits
                cur.close()
                if self.pooling:
                    self._pool.putconn(conn)

                break

            except psycopg2.OperationalError:

                _logger.warning(
                    "Lost connection to {url}; attempting reconnect".format(url=self.url))
                if self.pooling:
                    self._pool.closeall()
                self._connect()
                _logger.warning("Reconnected to {url}".format(url=self.url))

            n_tries_rem -= 1

        else:

            # N.B. Probably never reached
            raise HGVSError("Permanently lost connection to {url} ({n} retries)".format(
                url=self.url, n=n_retries))


class ParseResult(urlparse.ParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.

    """

    def __new__(cls, pr):
        return super(ParseResult, cls).__new__(cls, *pr)

    @property
    def database(self):
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self):
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None

    def __str__(self):
        return self.geturl()


def _parse_url(db_url):
    """parse database connection urls into components

    UTA database connection URLs follow that of SQLAlchemy, except
    that a schema may be optionally specified after the database. The
    skeleton format is:

       driver://user:pass@host/database/schema

    >>> params = _parse_url("driver://user:pass@host:9876/database/schema")

    >>> params.scheme
    u'driver'

    >>> params.hostname
    u'host'

    >>> params.username
    u'user'

    >>> params.password
    u'pass'

    >>> params.database
    u'database'

    >>> params.schema
    u'schema'

    """

    return ParseResult(urlparse.urlparse(db_url))
