# we're not testing hgvs.parser, but rather merely using it to load the
# grammar.  See test_hgvs_parser.py for the parser tests

import os
import unittest

from parsley import ParseError

import hgvs.parser

class Test_Parser(unittest.TestCase):
    longMessage = True

    def setUp(self):
        self.p = hgvs.parser.Parser()
        self.grammar = self.p._grammar

    def test_parser_basic(self):
        # num & snum
        self.assertEqual( self.grammar('55').num(), 55 )
        self.assertEqual( self.grammar('55').snum(), 55 )
        self.assertEqual( self.grammar('-55').snum(), -55 )

        # accn
        self.assertEqual( self.grammar('NM_01234.5').accn(), 'NM_01234.5' )
        self.assertEqual( self.grammar('NM_01234').accn(), 'NM_01234' )
        self.assertEqual( self.grammar('ENST01234').accn(), 'ENST01234' )
        self.assertEqual( self.grammar('LRG_1234').accn(), 'LRG_1234' )
        self.assertEqual( self.grammar('LRG_1234t1').accn(), 'LRG_1234t1' )


if __name__ == '__main__':
    unittest.main()
