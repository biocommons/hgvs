import unittest

import bdi.sources.uta0

from hgvs.projector import Projector
from hgvs.exceptions import *
import hgvs.location

class TestHgvsProjector(unittest.TestCase):
    def setUp(self):
        self._bdi = bdi.sources.uta0.connect()
        self.ref = 'GRCh37.p10'

    # Test combinations of these, both ways
    # MCL1, multiple transcripts, SNPs mapped by NCBI
    # http://tinyurl.com/len34jc
    # http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=12036617
    # NM_001197320.1:c.514G>A
    # NM_021960.4:c.973G>A
    # NM_182763.2:c.725G>A

    def _to_and_fro( self, pj, a, b ):
        ai = hgvs.location.Interval(
            start = hgvs.location.BaseOffsetPosition( base=a[0], datum=hgvs.location.CDS_START ),
            end   = hgvs.location.BaseOffsetPosition( base=a[1], datum=hgvs.location.CDS_START ),
            )
        bi = hgvs.location.Interval(
            start = hgvs.location.BaseOffsetPosition( base=b[0], datum=hgvs.location.CDS_START ),
            end   = hgvs.location.BaseOffsetPosition( base=b[1], datum=hgvs.location.CDS_START ),
            )
        
        self.assertEquals( pj.project_interval_forward(ai), bi )
        self.assertEquals( pj.project_interval_backward(bi), ai )
        

    def test_20_60(self):
        self._to_and_fro( Projector(self._bdi,self.ref,'NM_001197320.1','NM_021960.4'), (513,514), (972,973) )

    def test_20_63(self):
        self._to_and_fro( Projector(self._bdi,self.ref,'NM_001197320.1','NM_182763.2'), (513,514), (724,725) )

    def test_60_63(self):
        self._to_and_fro( Projector(self._bdi,self.ref,'NM_021960.4','NM_182763.2'), (972,973), (724,725) )

    def test_failures(self):
        self.assertRaises( HGVSError, Projector, self._bdi,self.ref,'NM_bogus','NM_bogus' )
        self.assertRaises( HGVSError, Projector, self._bdi,'bogus','NM_001197320.1','NM_021960.4' )
        

if __name__ == '__main__':
    unittest.main()
