import logging
import os
import sqlite3
import urlparse

#TODO: make dynamic import with importlib
import psycopg2
import psycopg2.extras

from .interface import Interface
from ..utils.aminoacids import seq_md5


localhost_db_url = 'postgresql://localhost/uta'
public_db_url = 'postgresql://uta_public:uta_public@uta.invitae.com/uta'
default_db_url = os.environ.get('UTA_DB_URL',public_db_url)
pg_schema = 'uta_20140210'

logger = logging.getLogger(__name__)


def connect(db_url=default_db_url):
    """
    Connect to a UTA database instance and return a UTA interface instance.

    When called with an explicit db_url argument, that db_url is used for connecting.

    When called without an explicit argument, the function default is
    determined by the environment variable UTA_DB_URL if it exists, or
    hgvs.datainterface.uta.public_db_url otherwise.

    The format of the db_url is driver://user:pass@host/database (the same
    as that used by SQLAlchemy).  Examples:

    A remote public postgresql database:
        postgresql://uta_public:uta_public@uta.invitae.com/uta'

    A local postgresql database:
        postgresql://localhost/uta
    
    A local SQLite database:
      sqlite:////tmp/uta-0.0.6.db
    """

    url = urlparse.urlparse(db_url)
    
    if url.scheme == 'sqlite':
        conn = UTA_sqlite(url)
    elif url.scheme == 'postgresql':
        conn = UTA_postgresql(url)
    else:
        # fell through connection scheme cases
        raise RuntimeError("{url.scheme} in {db_url} is not supported".format(url=url, db_url=db_url))

    conn.db_url = db_url
    return conn


class UTABase(Interface):
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
            join exon E on ES.exon_set_id=E.exon_set_id 
            where alt_ac=? and alt_aln_method=?
            group by tx_ac,alt_ac,alt_strand,alt_aln_method
            having max(end_i)>=? and min(start_i)<=?
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

    def __init__(self):
        self._check_schema_version('1')


    def data_version(self):
        cur = self._get_cursor()
        cur.execute("select * from meta where key = 'schema_version'")
        return cur.fetchone()['value']


    def schema_version(self):
        return self.data_version()


    ############################################################################
    ## Queries
    def get_acs_for_protein_seq(self,seq):
        """
        returns a list of protein accessions for a given sequence.  The
        list is guaranteed to contain at least one element with the
        MD5-based accession (MD5_01234abc...def56789) at the end of the
        list.
        """
        md5 = seq_md5(seq)
        cur = self._get_cursor()
        cur.execute(self.sql['acs_for_protein_md5'],[md5])
        return [ r['ac'] for r in cur.fetchall() ] + [ 'MD5_'+md5 ]


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
        cur = self._get_cursor()
        cur.execute(self.sql['gene_info'],[gene])
        return cur.fetchone()


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
        cur = self._get_cursor()
        cur.execute(self.sql['tx_exons'],[tx_ac, alt_ac, alt_aln_method])
        r = cur.fetchall()
        if len(r) == 0:
            return None
        else:
            return r


    def get_tx_for_gene(self,gene):
        """
        return transcript info records for supplied gene, in order of decreasing length

        :param gene: HGNC gene name
        :type gene: str
        """
        cur = self._get_cursor()
        cur.execute(self.sql['tx_for_gene'],[gene])
        return cur.fetchall()

    def get_tx_for_region(self,alt_ac,alt_aln_method,start_i,end_i):
        """
        return transcripts that overlap given region

        :param str alt_ac: reference sequence (e.g., NC_000007.13)
        :param str alt_aln_method: alignment method (e.g., splign)
        :param int start_i: 5' bound of region
        :param int end_i: 3' bound of region
        """
        cur = self._get_cursor()
        cur.execute(self.sql['tx_for_region'],[alt_ac,alt_aln_method,start_i,end_i])
        return cur.fetchall()

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
        cur = self._get_cursor()
        cur.execute(self.sql['tx_identity_info'],[tx_ac])
        return cur.fetchone()


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
        cur = self._get_cursor()
        cur.execute(self.sql['tx_info'],[tx_ac, alt_ac, alt_aln_method])
        return cur.fetchone()


    def get_tx_seq(self,ac):
        """return transcript sequence for supplied accession (ac), or None if not found

        :param ac: transcript accession with version (e.g., 'NM_000051.3')
        :type ac: str
        """
        cur = self._get_cursor()
        cur.execute(self.sql['tx_seq'],[ac])
        try:
            return cur.fetchone()['seq']
        except TypeError:
            return None


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
        cur = self._get_cursor()
        cur.execute(self.sql['tx_mapping_options'],[tx_ac])
        r = cur.fetchall()
        return r




    ############################################################################
    ## INTERNAL FUNCTIONS
    def _check_schema_version(self,required_version):
        obs_Mm = self.schema_version().split('.')[:2]
        req_Mm = required_version.split('.')[:2]
        if ( obs_Mm != req_Mm ):
            raise RuntimeError("Version mismatch: Version {req_Mm} required, but {self.db_url} is version {obs_Mm}".format(
                req_Mm = '.'.join(req_Mm), self=self, obs_Mm = '.'.join(obs_Mm)))
        # else no error



class UTA_sqlite(UTABase):
    def __init__(self,url):
        def _sqlite3_row_dict_factory(cur, row):
            "convert sqlite row to dict"
            return dict( (d[0],row[i]) for i,d in enumerate(cur.description) )

        self.db_url = url.geturl()

        db_path = url.path
        if not os.path.exists(db_path):
            raise IOError(db_path + ': Non-existent database file')
        self._con = sqlite3.connect(db_path)
        self._con.row_factory = _sqlite3_row_dict_factory
        logger.info("connected to "+db_path)
        super(UTA_sqlite,self).__init__()
        
    def _get_cursor(self):
        return self._con.cursor()


class UTA_postgresql(UTABase):
    def __init__(self,url):
        self.db_url = url.geturl()

        host = url.hostname
        port = url.port or 5432
        database = url.path[1:]
        user = url.username
        password= url.password

        logger.info('connecting to '+self.db_url)
        logger.debug('connection params: host={host}, port={port}, database={database}, user={user}, password={pw}'.format(
            host=host,port=port,database=database,user=user,pw='' if password is None else '***'))
        self._con = psycopg2.connect(host=host, port=port, database=database, user=user, password=password)
        super(UTA_postgresql,self).__init__()

        # remap sqlite's ? placeholders to psycopg2's %s
        self.sql = { k: v.replace('?','%s') for k,v in self.sql.iteritems() }
        
    def _get_cursor(self):
        cur = self._con.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cur.execute('set search_path = '+pg_schema)
        return cur




## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
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
