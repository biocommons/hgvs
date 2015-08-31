# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import inspect
import logging
import os
import sqlite3
import types
import urlparse

import psycopg2
import psycopg2.extras
import psycopg2.pool

from bioutils.digests import seq_md5

import hgvs
from ..dataproviders.interface import Interface
from ..decorators import deprecated
from ..decorators.lru_cache import lru_cache
from ..exceptions import HGVSError, HGVSDataNotAvailableError
from .seqfetcher import SeqFetcher

_url_key = os.environ.get('_UTA_URL_KEY', 'public' if hgvs._is_released_version else 'public_dev')
_default_db_url = os.environ.get('UTA_DB_URL', hgvs.global_config['uta'][_url_key])

_logger = logging.getLogger(__name__)


def connect(db_url=_default_db_url, pooling=False, application_name=None):
    """Connect to a UTA database instance and return a UTA interface instance.

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
    url = _parse_url(db_url)
    if url.scheme == 'sqlite':
        conn = UTA_sqlite(url)
    elif url.scheme == 'postgresql':
        conn = UTA_postgresql(url=url, pooling=pooling, application_name=application_name)
    else:
        # fell through connection scheme cases
        raise RuntimeError("{url.scheme} in {url} is not currently supported".format(url=url))
    return conn


class UTABase(Interface, SeqFetcher):
    required_version = "1.1"

    sql = {
        "acs_for_protein_md5": """
            select ac
            from seq_anno
            where seq_id=?
            """,
        "gene_info": """
            select *
            from gene
            where hgnc=?
            """,
        "tx_exons": """
            select *
            from tx_exon_aln_v
            where tx_ac=? and alt_ac=? and alt_aln_method=?
            order by alt_start_i
            """,
        "tx_for_gene": """
            select hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            from transcript T
            join exon_set ES on T.ac=ES.tx_ac where alt_aln_method != 'transcript' and hgnc=?
            """,
        "tx_for_region": """
            select tx_ac,alt_ac,alt_strand,alt_aln_method,min(start_i) as start_i,max(end_i) as end_i
            from exon_set ES
            join transcript T on ES.tx_ac = t.ac
            join exon E on ES.exon_set_id=E.exon_set_id 
            where alt_ac=? and alt_aln_method=?
            group by tx_ac,alt_ac,alt_strand,alt_aln_method
            having max(end_i)>? and min(start_i)<?
            """,
        "tx_identity_info": """
            select distinct(tx_ac), alt_ac, alt_aln_method, cds_start_i, cds_end_i, lengths, hgnc
            from tx_def_summary_v
            where tx_ac=?
            """,
        "tx_info": """
            select hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            from transcript T
            join exon_set ES on T.ac=ES.tx_ac
            where tx_ac=? and alt_ac=? and alt_aln_method=?
            """,
        "tx_mapping_options": """
            select distinct tx_ac,alt_ac,alt_aln_method
            from tx_exon_aln_v where tx_ac=? and exon_aln_id is not NULL
            """,
        "tx_seq": """
            select seq
            from seq S
            join seq_anno SA on S.seq_id=SA.seq_id
            where ac=?
            """,
        "tx_similar": """
            select *
            from tx_similarity_v
            where tx_ac1 = ?
            """,
        "tx_to_pro": """
            select * from associated_accessions where tx_ac = ?
            """,
    }

    def __init__(self, url):
        self.url = url
        self._connect()
        super(UTABase, self).__init__()
        _logger.info('connected to ' + str(self.url))

    def _execute(self, sql, *args):
        cur = self._get_cursor()
        cur.execute(sql, *args)
        return cur

    ############################################################################
    ## Queries

    @lru_cache(maxsize=1)
    def data_version(self):
        return self.url.schema

    @lru_cache(maxsize=1)
    def schema_version(self):
        cur = self._execute("select * from meta where key = 'schema_version'")
        return cur.fetchone()['value']

    @lru_cache(maxsize=128)
    def get_acs_for_protein_seq(self, seq):
        """
        returns a list of protein accessions for a given sequence.  The
        list is guaranteed to contain at least one element with the
        MD5-based accession (MD5_01234abc...def56789) at the end of the
        list.
        """
        md5 = seq_md5(seq)
        cur = self._execute(self.sql['acs_for_protein_md5'], [md5])
        return [r['ac'] for r in cur.fetchall()] + ['MD5_' + md5]

    @lru_cache(maxsize=128)
    def get_gene_info(self, gene):
        """
        returns basic information about the gene.

        :param gene: HGNC gene name
        :type gene: str

        # database results
        hgnc    | ATM
        maploc  | 11q22-q23
        descr   | ataxia telangiectasia mutated
        summary | The protein encoded by this gene belongs to the PI3/PI4-kinase family. This...
        aliases | AT1,ATA,ATC,ATD,ATE,ATDC,TEL1,TELO1
        added   | 2014-02-04 21:39:32.57125

        """
        cur = self._execute(self.sql['gene_info'], [gene])
        return cur.fetchone()

    @lru_cache(maxsize=128)
    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        """
        return transcript exon info for supplied accession (tx_ac, alt_ac, alt_aln_method), or None if not found

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        :param alt_ac: specific genomic sequence (e.g., NC_000011.4)
        :type alt_ac: str

        :param alt_aln_method: sequence alignment method (e.g., splign, blat)
        :type alt_aln_method: str

        # tx_exons = db.get_tx_exons('NM_199425.2', 'NC_000020.10', 'splign')
        # len(tx_exons)
        3

        tx_exons have the following attributes::

            {
                'tes_exon_set_id' : 98390
                'aes_exon_set_id' : 298679
                'tx_ac'           : 'NM_199425.2'
                'alt_ac'          : 'NC_000020.10'
                'alt_strand'      : -1
                'alt_aln_method'  : 'splign'
                'ord'             : 2
                'tx_exon_id'      : 936834
                'alt_exon_id'     : 2999028
                'tx_start_i'      : 786
                'tx_end_i'        : 1196
                'alt_start_i'     : 25059178
                'alt_end_i'       : 25059588
                'cigar'           : '410='
            }

        For example:

        # tx_exons[0]['tx_ac']
        'NM_199425.2'

        """
        cur = self._execute(self.sql['tx_exons'], [tx_ac, alt_ac, alt_aln_method])
        rows = cur.fetchall()
        if len(rows) == 0:
            raise HGVSDataNotAvailableError(
                "No tx_exons for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                    tx_ac=tx_ac,
                    alt_ac=alt_ac,
                    alt_aln_method=alt_aln_method))
        return rows

    @lru_cache(maxsize=128)
    def get_tx_for_gene(self, gene):
        """
        return transcript info records for supplied gene, in order of decreasing length

        :param gene: HGNC gene name
        :type gene: str
        """
        cur = self._execute(self.sql['tx_for_gene'], [gene])
        return cur.fetchall()

    @lru_cache(maxsize=128)
    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        """
        return transcripts that overlap given region

        :param str alt_ac: reference sequence (e.g., NC_000007.13)
        :param str alt_aln_method: alignment method (e.g., splign)
        :param int start_i: 5' bound of region
        :param int end_i: 3' bound of region
        """
        cur = self._execute(self.sql['tx_for_region'], [alt_ac, alt_aln_method, start_i, end_i])
        return cur.fetchall()

    @lru_cache(maxsize=128)
    def get_tx_identity_info(self, tx_ac):
        """returns features associated with a single transcript.

        :param tx_ac: transcript accession with version (e.g., 'NM_199425.2')
        :type tx_ac: str

        # database output
        -[ RECORD 1 ]--+-------------
        tx_ac          | NM_199425.2
        alt_ac         | NM_199425.2
        alt_aln_method | transcript
        cds_start_i    | 283
        cds_end_i      | 1003
        lengths        | {707,79,410}
        hgnc           | VSX1

        """
        cur = self._execute(self.sql['tx_identity_info'], [tx_ac])
        rows = cur.fetchall()
        if len(rows) == 0:
            raise HGVSDataNotAvailableError(
                "No transcript definition for (tx_ac={tx_ac})".format(
                    tx_ac=tx_ac))
        return rows[0]


    @lru_cache(maxsize=128)
    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        """return a single transcript info for supplied accession (tx_ac, alt_ac, alt_aln_method), or None if not found

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        :param alt_ac: specific genomic sequence (e.g., NC_000011.4)
        :type alt_ac: str

        :param alt_aln_method: sequence alignment method (e.g., splign, blat)
        :type alt_aln_method: str

        # database output
        -[ RECORD 1 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | AC_000143.1
        alt_aln_method | splign

        """
        cur = self._execute(self.sql['tx_info'], [tx_ac, alt_ac, alt_aln_method])
        rows = cur.fetchall()
        if len(rows) == 0:
            raise HGVSDataNotAvailableError(
                "No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                    tx_ac=tx_ac,
                    alt_ac=alt_ac,
                    alt_aln_method=alt_aln_method))
        elif len(rows) == 1:
            return rows[0]
        else:
            raise HGVSError(
                "Multiple ({n}) replies for tx_info(tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                    n=len(rows),
                    tx_ac=tx_ac,
                    alt_ac=alt_ac,
                    alt_aln_method=alt_aln_method))

    @lru_cache(maxsize=128)
    def get_tx_mapping_options(self, tx_ac):
        """Return all transcript alignment sets for a given transcript
        accession (tx_ac); returns empty list if transcript does not
        exist.  Use this method to discovery possible mapping options
        supported in the database

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        # database output
        -[ RECORD 1 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | AC_000143.1
        alt_aln_method | splign
        -[ RECORD 2 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | NC_000011.9
        alt_aln_method | blat

        """
        cur = self._execute(self.sql['tx_mapping_options'], [tx_ac])
        rows = cur.fetchall()
        return rows

    @lru_cache(maxsize=128)
    def get_similar_transcripts(self, tx_ac):
        """Return a list of transcripts that are similar to the given
        transcript, with relevant similarity criteria.

        >> sim_tx = hdp.get_similar_transcripts('NM_001285829.1')
        >> dict(sim_tx[0])
        { 'cds_eq': False,
        'cds_es_fp_eq': False,
        'es_fp_eq': True,
        'tx_ac1': 'NM_001285829.1',
        'tx_ac2': 'ENST00000498907' }

        where:

        * cds_eq means that the CDS sequences are identical
        * es_fp_eq means that the full exon structures are identical
          (i.e., incl. UTR)
        * cds_es_fp_eq means that the cds-clipped portions of the exon
          structures are identical (i.e., ecluding. UTR)
        * Hint: "es" = "exon set", "fp" = "fingerprint", "eq" = "equal"

        "exon structure" refers to the start and end coordinates on a
        specified reference sequence. Thus, having the same exon
        structure means that the transcripts are defined on the same
        reference sequence and have the same exon spans on that
        sequence.

        """

        cur = self._execute(self.sql['tx_similar'], [tx_ac])
        rows = cur.fetchall()
        return rows

    @lru_cache(maxsize=128)
    def get_pro_ac_for_tx_ac(self, tx_ac):
        """Return the (single) associated protein accession for a given transcript
        accession, or None if not found."""

        cur = self._execute(self.sql['tx_to_pro'], [tx_ac])
        rows = cur.fetchall()
        try:
            return rows[0]['pro_ac']
        except IndexError:
            return None

    # Sequence fetching
    # -----------------
    # UTA stored a subset of relevant sequences in
    # the postgresql database, but that's impractical for large
    # sequences (genome scale) and is difficult to maintain.
    # The design goal is to migrate from the get_tx_seq() method to
    # externalize all sequence fetching.  See
    # TODO: Externalize sequence fetching (https://bitbucket.org/biocommons/hgvs/issue/236/)
    #
    # For the 0.4.0 release, we'll enable a a new method, fetch_seq(),
    # and deprecate get_tx_seq(). get_tx_seq() will be removed in a
    # subsequent major inteface update.  fetch_seq() itself will wrap
    # get_tx_seq() and SeqFetcher.fetch_seq() now, but is expected to
    # be entirely replaced by a more complete sequence database.
    # See https://bitbucket.org/biocommons/hgvs/issue/240/

    @lru_cache(maxsize=128)
    def fetch_seq(self, ac, start_i=None, end_i=None):
        """Fetches sequence by accession, optionally bounded by [start_i,end_i).
        See SeqFetcher.fetch_seq() for details and examples.

        This function tries _get_tx_seq() (because it's usually
        faster), and then SeqFetcher.fetch_seq().
        """

        if any(ac.startswith(pfx) for pfx in ['NM_', 'NR_', 'ENST']):
            try:
                seq = self._get_tx_seq(ac)[start_i:end_i]
                _logger.debug("fetched {ac} from UTA".format(ac=ac))
                return seq
            except HGVSDataNotAvailableError:
                pass
        # if ac not matching or on HGVSDataNotAvailableError...
        seq = super(UTABase, self).fetch_seq(ac, start_i, end_i)
        _logger.debug("fetched {ac} with SeqFetcher".format(ac=ac))
        assert seq is not None
        return seq

    # TODO: Remove get_tx_seq() in 0.5.0
    @deprecated(use_instead="fetch_seq(...)")
    def get_tx_seq(self, ac):
        """DEPRECATED: will be removed in 0.5.0"""
        return self._get_tx_seq(ac)

    def _get_tx_seq(self, ac):
        """return transcript sequence for supplied accession (ac), or None if not found

        :param ac: transcript accession with version (e.g., 'NM_000051.3')
        :type ac: str
        """
        cur = self._execute(self.sql['tx_seq'], [ac])
        row = cur.fetchone()
        if row and row['seq'] is not None:
            return row['seq']
        raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))


class UTA_postgresql(UTABase):
    def __init__(self, url, pooling=False, application_name=None):
        self.pooling = pooling
        self.application_name = application_name
        if url.schema is None:
            raise Exception("No schema name provided in {url}".format(url=url))
        super(UTA_postgresql, self).__init__(url)

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
            self._pool = psycopg2.pool.ThreadedConnectionPool(1, 10, **conn_args)
        else:
            self._conn = psycopg2.connect(**conn_args)
            self._conn.autocommit = True

        self._ensure_schema_exists()

        # remap sqlite's ? placeholders to psycopg2's %s
        self.sql = {k: v.replace('?', '%s') for k, v in self.sql.iteritems()}

    def _ensure_schema_exists(self):
        # N.B. On AWS RDS, information_schema.schemata always returns zero rows
        cur = self._execute("select exists(SELECT 1 FROM pg_namespace WHERE nspname = %s)", [self.url.schema])
        if cur.fetchone()[0]:
            return
        raise HGVSDataNotAvailableError("specified schema ({}) does not exist (url={})".format(
            self.url.schema, self.url))

    def _get_cursor(self):
        """returns a cursor obtained from a single or pooled connection, and
        sets the postgresql search_path appropriately

        Although connections are threadsafe, cursors are bound to
        connections and are not threadsafe. Therefore, be sure to not
        share cursors across threads. Cursors are reaped when they go
        out of scope, but they may be cleaned up manually with
        cursor.close(). See the psycopg2 docs for cursor lifetime

        """

        try:
            conn = self._pool.getconn() if self.pooling else self._conn
            conn.autocommit = True
            cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            return cur
        finally:
            if self.pooling:
                self._pool.putconn(conn)

    def _execute(self, sql, *args):
        cur = self._get_cursor()
        # N.B. The search_path is reset when pooling is used and the
        # connection is returned with putconn. Therefore, we need to
        # set it for every execution.
        cur.execute("set search_path = " + self.url.schema + ";")
        cur.execute(sql, *args)
        return cur


class UTA_sqlite(UTABase):
    # TODO: implement mocks (issue #237) The current sqlite db was
    # based on schema v1. No tests currently use 1.1 features from
    # sqlite, so we override the required_version here.
    required_version = "1"

    def _connect(self):
        def _sqlite3_row_dict_factory(cur, row):
            "convert sqlite row to dict"
            return dict((d[0], row[i]) for i, d in enumerate(cur.description))

        if not os.path.exists(self.url.path):
            raise IOError(self.url.path + ': Non-existent database file')
        self._conn = sqlite3.connect(self.url.path)
        self._conn.row_factory = _sqlite3_row_dict_factory

    def _get_cursor(self):
        return self._conn.cursor()


class ParseResult(urlparse.ParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.

    """

    def __new__(cls, pr):
        pr.__class__ = cls
        return pr

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
