import requests

class EutilsSeqSlicerMixin(object):
    """Provides seq_slice() method to fetch sequences; this is a placeholder until seqdb becomes available

    >>> es = EutilsSeqSlicerMixin()

    >>> es.seq_slice('NP_056374.2',0,10)
    'MESRETLSSS'

    >>> es.seq_slice('NG_032072.1',0,10)
    'AAAATTAAAT'

    >>> es.seq_slice('NC_000001.10',2000000,2000030)
    'ATCACACGTGCAGGAACCCTTTTCCAAAGG'

    >>> es.seq_slice('NT_113901.1',20,30)
    'TTCTTAAGCT'

    >>> es.seq_slice('NW_003571030.1',20,30)
    'AGCAATATGG'

    """

    def seq_slice(self, ac, start_i, end_i):
        """return a subsequence of the given accession for the
        interbase interval [start_i,end_i)"""

        url_fmt = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={ac}&rettype=fasta&seq_start={start}&seq_stop={stop}"
        url = url_fmt.format(ac=ac, start=start_i+1, stop=end_i)
        resp = requests.get(url)
        resp.raise_for_status()
        return ''.join(resp.content.splitlines()[1:])
        
