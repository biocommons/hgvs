import os
import unittest

from parsley import ParseError
import hgvs.parser

class Test_Position(unittest.TestCase):
    longMessage = True
    
    def setUp(self):
        self.parser = hgvs.parser.Parser()

    def test_Parser_roundtrip(self):
        pos_tests = """
NC_000001.10:g.150549917_150549920delinsTTTGCAA
NC_000010.10:g.99511133_99511137del5
NC_000016.9:g.28489143_28489146del
NC_000023.10:g.106890924_106890925insAA
NC_000023.10:g.107939544_107939551delCGCTGAAA
NM_000495.4:c.4994_5001delCGCTGAAA
NM_000694.2:c.752T>A
NM_001002261.3:c.805_809del5
NM_001031681.2:c.62-1G>A
NM_001032393.2:c.*515C>G
NM_001032393.2:c.-515C>G
NM_001042432.1:c.1109_1112del
NM_001197225.2:c.-183C>T
NM_002764.3:c.793_794insAA
NM_021960.4:c.984_987delinsTTGCAAA
NM_024944.2:c.389+1G>A
NM_144588.6:c.805-15_805-11del5
NP_003768.2:p.Ala4518_Ala5423delinsGln
NP_003768.2:p.Ile4452Asnfs*4
NP_003768.2:p.Leu1511Serfs*10
NP_003768.2:p.Leu4363del
NP_003768.2:p.Lys3045Thrfs*14
NP_003768.2:p.Arg4412Alafs*2
NP_003768.2:p.Thr4509Argfs*15
NP_003768.2:p.Leu2490_Thr2496del
NP_003768.2:p.Val4487Cysfs*3
NP_003768.2:p.Ser2176*
NP_003768.2:p.Arg2852*
NP_003768.2:p.Arg2907*
""".strip().split('\n')
        for t in pos_tests:
            self.assertEqual( t, str(self.parser.parse(t)) )


    def test_Parser_reject(self):
        neg_tests = """
NC_000001.10:g.155208383_155208384dup2
NM_001005741.2:c.512_513dup2
""".strip().split('\n')
        for t in neg_tests:
            with self.assertRaises(Exception, msg=t):
                self.parser.parse(t)

    def test_parser_hgvs_gauntlet(self):
        fn = os.path.join( os.path.dirname(__file__), 'data', 'hgvs-gauntlet' )
        for var in open(fn,'r'):
            var = var.strip()
            v = self.parser._grammar(var).hgvs_variant()
            self.assertEqual( str(v), var )

if __name__ == '__main__':
    unittest.main()
