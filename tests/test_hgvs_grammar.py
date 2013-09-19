import unittest

from parsley import ParseError
import hgvs.parser

class Test_Grammar(unittest.TestCase):
    longMessage = True

    def setUp(self):
        self.p = hgvs.parser.Parser()

    def test_grammar_basic(self):
        # num & snum
        self.assertEqual( self.p._grammar('55').num(), 55 )
        self.assertEqual( self.p._grammar('55').snum(), 55 )
        self.assertEqual( self.p._grammar('-55').snum(), -55 )

        # accn
        self.assertEqual( self.p._grammar('NM_01234.5').accn(), 'NM_01234.5' )
        self.assertEqual( self.p._grammar('NM_01234').accn(), 'NM_01234' )
        self.assertEqual( self.p._grammar('ENST01234').accn(), 'ENST01234' )
        self.assertEqual( self.p._grammar('LRG_1234').accn(), 'LRG_1234' )
        self.assertEqual( self.p._grammar('LRG_1234t1').accn(), 'LRG_1234t1' )


if __name__ == '__main__':
    unittest.main()
