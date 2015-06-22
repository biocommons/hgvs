import re

import requests


def _fetch_seq_ensembl(ac, start_i=None, end_i=None):
    url_fmt = "http://rest.ensembl.org/sequence/id/{ac}"
    url = url_fmt.format(ac=ac, start=start_i+1, stop=end_i)
    r = requests.get(url, headers={"Content-Type" : "application/json"})
    r.raise_for_status()
    return r.json()['seq'][start_i:end_i]


def _fetch_seq_ncbi(ac, start_i=None, end_i=None):
    url_fmt = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={ac}&rettype=fasta&seq_start={start}&seq_stop={stop}"
    url = url_fmt.format(ac=ac, start=start_i+1, stop=end_i)
    resp = requests.get(url)
    resp.raise_for_status()
    return ''.join(resp.content.splitlines()[1:])


class SeqFetcher(object):
    """Provides fetch_seq() method to fetch sequences from NCBI eutils and
    Ensembl REST interfaces

    >>> ssm = SeqFetcher()

    >>> ssm.fetch_seq('NP_056374.2',0,10)
    'MESRETLSSS'

    >>> ssm.fetch_seq('NG_032072.1',0,10)
    'AAAATTAAAT'

    >>> ssm.fetch_seq('NC_000001.10',2000000,2000030)
    'ATCACACGTGCAGGAACCCTTTTCCAAAGG'

    >>> ssm.fetch_seq('NT_113901.1',20,30)
    'TTCTTAAGCT'

    >>> ssm.fetch_seq('NW_003571030.1',20,30)
    'AGCAATATGG'

    >>> ssm.fetch_seq('ENSP00000288602',0,10)
    u'MAALSGGGGG'

    >>> ssm.fetch_seq('ENST00000288602',0,10)
    u'CGCCTCCCTT'
    
    """

    _ac_dispatch = [
        {'re': re.compile('^ENS[TP]\d+'), 'fetcher': _fetch_seq_ensembl},
        {'re': re.compile('^(?:N[CGPRTW]|AC|GL)_|^U\d+'), 'fetcher': _fetch_seq_ncbi},
        ]

    def fetch_seq(self, ac, start_i=None, end_i=None):
        """return a subsequence of the given accession for the
        interbase interval [start_i,end_i)"""
        for drec in self._ac_dispatch:
            if drec['re'].match(ac):
                return drec['fetcher'](ac, start_i, end_i)
        raise Exception("No sequence fetcher identified for " + ac)
