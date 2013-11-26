import unittest

import uta.db.transcriptdb

import hgvs.edti.uta0

class Test_EDTI_UTA0(unittest.TestCase):
    def setUp(self):
        self.edti = hgvs.edti.uta0.UTA0( uta.db.transcriptdb.TranscriptDB() )
        return self

    def test_fetch_gene_info(self):
        gi = self.edti.fetch_gene_info('VHL')

    def test_fetch_gene_transcripts(self):
        g_txs = self.edti.fetch_gene_transcripts('VHL')
        self.assertEqual( len(g_txs), 2)
        self.assertEqual( g_txs[0]['chr'], '3' )


if __name__ == '__main__':
    unittest.main()
        
