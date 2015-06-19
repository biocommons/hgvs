# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import os
import sqlite3
import types
import urlparse

#TODO: make dynamic import with importlib
import psycopg2
import psycopg2.extras
import psycopg2.pool

from bioutils.digests import seq_md5

from ..dataproviders.interface import Interface
from ..decorators.lru_cache import lru_cache
from .ncbi import EutilsSeqSlicerMixin


_uta_urls = {
    # these are provided for developer convenience
    "local": "postgresql://localhost/uta/uta_20140210",
    "local-dev": "postgresql://localhost/uta_dev/uta1",
    "public": "postgresql://uta_public:uta_public@uta.invitae.com/uta/uta_20140210",
    "public-dev": "postgresql://uta_public:uta_public@uta.invitae.com/uta_dev/uta1",
    "sqlite-dev": "sqlite:/home/reece/projects/biocommons/hgvs/tests/db/uta-test-1.db",
}

default_db_url = os.environ.get('UTA_DB_URL', _uta_urls["public"])


def connect(db_url=default_db_url, pooling=False):
    """Connect to a UTA database instance and return a UTA interface instance.

    :param db_url: URL for database connection
    :type db_url: string
    :param pooling: whether to use connection pooling (postgresql only)
    :type pooling: bool

    When called with an explicit db_url argument, that db_url is used for connecting.

    When called without an explicit argument, the function default is
    determined by the environment variable UTA_DB_URL if it exists, or
    hgvs.datainterface.uta.public_db_url otherwise.

    >>> hdp = connect()
    >>> hdp.schema_version()
    '1'

    The format of the db_url is driver://user:pass@host/database (the same
    as that used by SQLAlchemy).  Examples:

    A remote public postgresql database:
        postgresql://uta_public:uta_public@uta.invitae.com/uta'

    A local postgresql database:
        postgresql://localhost/uta

    A local SQLite database:
      sqlite:////tmp/uta-0.0.6.db

    For postgresql db_urls, pooling=True causes connect to use a
    psycopg2.pool.ThreadedConnectionPool.
    """

    url = _parse_url(db_url)
    if url.scheme == 'sqlite':
        conn = UTA_sqlite(url)
    elif url.scheme == 'postgresql':
        conn = UTA_postgresql(url, pooling)
    else:
        # fell through connection scheme cases
        raise RuntimeError("{url.scheme} in {url} is not currently supported".format(
            url=url))
    return conn


class UTABase(Interface, EutilsSeqSlicerMixin):
    required_version = "1.0"

    _logger = logging.getLogger(__name__)

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
    }

    def __init__(self, url):
        self.url = url
        self._connect()
        super(UTABase,self).__init__()
        self._logger.info('connected to ' + str(self.url))

    def _execute(self, sql, *args):
        cur = self._get_cursor()
        cur.execute(sql, *args)
        return cur


    ############################################################################
    ## Queries

    @lru_cache(maxsize=1)
    def data_version(self):
        cur = self._execute("select * from meta where key = 'schema_version'")
        return cur.fetchone()['value']


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
        r = cur.fetchall()
        if len(r) == 0:
            return None
        else:
            return r


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
        return cur.fetchone()


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
        return cur.fetchone()


    @lru_cache(maxsize=128)
    def get_tx_seq(self, ac):
        """return transcript sequence for supplied accession (ac), or None if not found

        :param ac: transcript accession with version (e.g., 'NM_000051.3')
        :type ac: str
        """
        cur = self._execute(self.sql['tx_seq'], [ac])
        try:
            return cur.fetchone()['seq']
        except TypeError:
            return None


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
        r = cur.fetchall()
        return r



class UTA_postgresql(UTABase):
    def __init__(self, url, pooling=False):
        self.pooling = pooling
        if url.schema is None:
            raise Exception("No schema name provided in {url}".format(url=url))
        super(UTA_postgresql, self).__init__(url)

    def _connect(self):
        if self.pooling:
            self._pool = psycopg2.pool.ThreadedConnectionPool(1, 10,
                                                              host=self.url.hostname,
                                                              port=self.url.port,
                                                              database=self.url.database,
                                                              user=self.url.username,
                                                              password=self.url.password)
        else:
            self._conn = psycopg2.connect(host=self.url.hostname,
                                          port=self.url.port,
                                          database=self.url.database,
                                          user=self.url.username,
                                          password=self.url.password)

        # remap sqlite's ? placeholders to psycopg2's %s
        self.sql = {k: v.replace('?', '%s') for k, v in self.sql.iteritems()}


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

    >>> from pprint import pprint

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
