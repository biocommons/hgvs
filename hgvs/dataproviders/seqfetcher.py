import re

import requests

from ..decorators.lru_cache import lru_cache


# TODO: Move this fetching to the eutils package, which provides
# throttling.
@lru_cache(maxsize=100)
def _fetch_seq_ensembl(ac, start_i=None, end_i=None):
    """Fetch the specified sequence slice from Ensembl using the public
    REST interface.

    An interbase interval may be optionally provided with start_i and
    end_i. However, the Ensembl REST interface does not currently
    accept intervals, so the entire sequence is returned and sliced
    locally.

    >>> len(_fetch_seq_ensembl('ENSP00000288602'))
    766

    >>> _fetch_seq_ensembl('ENSP00000288602',0,10)
    u'MAALSGGGGG'

    >>> _fetch_seq_ensembl('ENSP00000288602')[0:10]
    u'MAALSGGGGG'

    >>> _fetch_seq_ensembl('ENSP00000288602',0,10) == _fetch_seq_ensembl('ENSP00000288602')[0:10]
    True

    """

    url_fmt = "http://rest.ensembl.org/sequence/id/{ac}"
    url = url_fmt.format(ac=ac)
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    seq = r.json()['seq']
    return seq if (start_i is None or end_i is None) else seq[start_i:end_i]


@lru_cache(maxsize=100)
def _fetch_seq_ncbi(ac, start_i=None, end_i=None):
    """Fetch sequences from NCBI using the eutils interface.  

    An interbase interval may be optionally provided with start_i and
    end_i. NCBI eutils will return just the requested subsequence,
    which might greatly reduce payload sizes (especially with
    chromosome-scale sequences).

    >>> len(_fetch_seq_ncbi('NP_056374.2'))
    1596

    Pass the desired interval rather than using Python's [] slice
    operator.

    >>> _fetch_seq_ncbi('NP_056374.2',0,10)
    'MESRETLSSS'

    >>> _fetch_seq_ncbi('NP_056374.2')[0:10]
    'MESRETLSSS'

    >>> _fetch_seq_ncbi('NP_056374.2',0,10) == _fetch_seq_ncbi('NP_056374.2')[0:10]
    True

    """

    url_fmt = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={ac}&rettype=fasta"
    if (start_i is None or end_i is None):
        url = url_fmt.format(ac=ac)
    else:
        url_fmt += "&seq_start={start}&seq_stop={stop}"
        url = url_fmt.format(ac=ac, start=start_i + 1, stop=end_i)
    resp = requests.get(url)
    resp.raise_for_status()
    return ''.join(resp.content.splitlines()[1:])


class SeqFetcher(object):
    """Provides fetch_seq() method to fetch sequences from NCBI eutils and
    Ensembl REST interfaces. This class is primarily intended as a
    mixin for HGVS data providers that doen't otherwise have access to
    sequence data.

    >>> sf = SeqFetcher()

    Without an interval, the full sequence is returned::

    >>> len(sf.fetch_seq('NP_056374.2'))
    1596

    Therefore, it's preferable to provide the interval rather than
    using Python slicing sequence on the delivered sequence::

    >>> sf.fetch_seq('NP_056374.2',0,10)   # This!
    'MESRETLSSS'

    >>> sf.fetch_seq('NP_056374.2')[0:10]  # Not this!
    'MESRETLSSS'

    >>> sf.fetch_seq('NP_056374.2',0,10) == sf.fetch_seq('NP_056374.2')[0:10]
    True

    Providing intervals is especially important for large sequences::

    >>> sf.fetch_seq('NC_000001.10',2000000,2000030)
    'ATCACACGTGCAGGAACCCTTTTCCAAAGG'

    Essentially any RefSeq, Genbank, BIC, or Ensembl sequence may be
    fetched:

    >> [(ac,sf.fetch_seq(ac,0,25))
    ... for ac in ['NG_032072.1', 'NW_003571030.1', 'NT_113901.1', 'NC_000001.10', 'NP_056374.2', 'GL000191.1', 'KB663603.1',
    ...            'ENST00000288602', 'ENSP00000288602']]
    [('NG_032072.1', 'AAAATTAAATTAAAATAAATAAAAA'),
     ('NW_003571030.1', 'TTGTGTGTTAGGGTGCTCTAAGCAA'),
     ('NT_113901.1', 'GAATTCCTCGTTCACACAGTTTCTT'),
     ('NC_000001.10', 'NNNNNNNNNNNNNNNNNNNNNNNNN'),
     ('NP_056374.2', 'MESRETLSSSRQRGGESDFLPVSSA'),
     ('GL000191.1', 'GATCCACCTGCCTCAGCCTCCCAGA'),
     ('KB663603.1', 'TTTATTTATTTTAGATACTTATCTC'),
     ('ENST00000288602', u'CGCCTCCCTTCCCCCTCCCCGCCCG'),
     ('ENSP00000288602', u'MAALSGGGGGGAEPGQALFNGDMEP')]

    """

    _ac_dispatch = [
        {'re': re.compile('^(?:AC|N[CGMPRTW])_|^[A-L]\w\d|^U\d'), 'fetcher': _fetch_seq_ncbi},
        {'re': re.compile('^ENS[TP]\d+'), 'fetcher': _fetch_seq_ensembl},
        ]

    def fetch_seq(self, ac, start_i=None, end_i=None):
        """return a subsequence of the given accession for the
        interbase interval [start_i,end_i)"""
        for drec in self._ac_dispatch:
            if drec['re'].match(ac):
                return drec['fetcher'](ac, start_i, end_i)
        raise Exception("No sequence fetcher identified for " + ac)
