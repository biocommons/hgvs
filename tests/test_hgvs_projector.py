import unittest

import bdi.sources.uta0

from hgvs.exceptions import *
import hgvs.location
import hgvs.parser
import hgvs.projector

class TestHgvsProjector(unittest.TestCase):
    def setUp(self):
        self.bdi = bdi.sources.uta0.connect()
        self.ref = 'GRCh37.p10'
        self.hp = hgvs.parser.Parser()
        
    def tst_forward_and_backward(self,v1,v2):
        pj = hgvs.projector.Projector(self.bdi,self.ref,v1.ac,v2.ac)
        self.assertEqual( pj.project_variant_forward(v1), v2 )
        self.assertEqual( pj.project_variant_backward(v2), v1 )

    def test_rs201430561(self):
        # rs201430561 http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=201430561
        hgvs_c = ['NM_001197320.1:c.281C>T','NM_021960.4:c.740C>T','NM_182763.2:c.688+403C>T']
        var_c = [ self.hp.parse_hgvs_variant(h) for h in hgvs_c ]
        self.tst_forward_and_backward( var_c[0], var_c[1] )
        self.tst_forward_and_backward( var_c[0], var_c[2] )
        self.tst_forward_and_backward( var_c[1], var_c[2] )

    def test_bad_acs(self):
        hgvs_c = ['NM_001197320.1:c.281C>T','NM_021960.4:c.740C>T','NM_182763.2:c.688+403C>T']
        var_c = [ self.hp.parse_hgvs_variant(h) for h in hgvs_c ]
        pj = hgvs.projector.Projector(self.bdi,self.ref,var_c[0].ac,var_c[1].ac)
        # intentionally call p_v_f with variant on *destination* transcript, and vice versa
        self.assertRaises( RuntimeError, pj.project_variant_forward, var_c[1] )
        self.assertRaises( RuntimeError, pj.project_variant_backward, var_c[0] )


if __name__ == '__main__':
    unittest.main()
