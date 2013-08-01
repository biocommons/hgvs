import unittest

import hgvs.parser

class Test_Position(unittest.TestCase):
    tests = """NC_000001.10:g.150549917_150549920delinsTTTGCAA
NC_000001.10:g.155208383_155208384dup2
NC_000010.10:g.99511133_99511137del5
NC_000016.9:g.28489143_28489146del
NC_000023.10:g.106890924_106890925insAA
NC_000023.10:g.107939544_107939551delCGCTGAAA
NM_000495.4:c.4994_5001delCGCTGAAA
NM_000694.2:c.752T>A
NM_001002261.3:c.805_809del5
NM_001005741.2:c.512_513dup2
NM_001031681.2:c.62-1G>A
NM_001032393.2:c.*515C>G
NM_001032393.2:c.-515C>G
NM_001042432.1:c.1109_1112del
NM_001197225.2:c.-183C>T
NM_002764.3:c.793_794insAA
NM_021960.4:c.984_987delinsTTGCAAA
NM_024944.2:c.389+1G>A
NM_144588.6:c.805-15_805-11del5""".split('\n')

    def test_Parser_roundtrip(self):
        p = hgvs.parser.Parser()
        for t in self.tests:
            self.assertEqual( t, str(p.parse(t)) )


if __name__ == '__main__':
    unittest.main()
