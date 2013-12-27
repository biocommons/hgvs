import unittest

from bdi.sources.uta0_pg import UTA0

import hgvs.parser
import hgvs.hgvsmapper

class Test_HGVSMapper(unittest.TestCase):
    def setUp(self):
        self.bdi = UTA0()
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()

    def test_1(self):
        # From garcia.tsv:
        # AOAH    NM_001177507.1:c.1486G>A      
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'
    
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hm.hgvsg_to_hgvsc( var_g, 'NM_001637.3' )

        self.assertEqual( str(var_c) , hgvs_c )


if __name__ == '__main__':
    unittest.main()
