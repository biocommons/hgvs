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

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
