# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

# we're not testing hgvs.parser, but rather merely using it to load the
# grammar.  See test_hgvs_parser.py for the parser tests

import unittest

from nose.plugins.attrib import attr

import hgvs.parser


@attr(tags=["quick"])
class Test_Parser(unittest.TestCase):
    longMessage = True

    def setUp(self):
        self.p = hgvs.parser.Parser()
        self.grammar = self.p._grammar

    def test_parser_basic(self):
        # num & snum
        self.assertEqual(self.grammar('55').num(), 55)
        self.assertEqual(self.grammar('55').snum(), 55)
        self.assertEqual(self.grammar('-55').snum(), -55)

        # accn
        self.assertEqual(self.grammar('NM_01234.5').accn(), 'NM_01234.5')
        self.assertEqual(self.grammar('NM_01234').accn(), 'NM_01234')
        self.assertEqual(self.grammar('ENST01234').accn(), 'ENST01234')
        self.assertEqual(self.grammar('LRG_1234').accn(), 'LRG_1234')
        self.assertEqual(self.grammar('LRG_1234t1').accn(), 'LRG_1234t1')

    def test_parser_aa(self):
        # terms
        self.assertEqual(self.p.parse_term1('*'), '*')
        self.assertEqual(self.p.parse_term3('Ter'), 'Ter')
        self.assertEqual(self.p.parse_term13('*'), '*')
        self.assertEqual(self.p.parse_term13('Ter'), 'Ter')

        # AAs
        self.assertEqual(self.p.parse_aa1('A'), 'A')
        self.assertEqual(self.p.parse_aa3('Ala'), 'Ala')
        self.assertEqual(self.p.parse_aa13('A'), 'A')
        self.assertEqual(self.p.parse_aa13('Ala'), 'Ala')

        # mixto
        self.assertEqual(self.p.parse_aat1('*'), '*')
        self.assertEqual(self.p.parse_aat1('A'), 'A')
        self.assertEqual(self.p.parse_aat3('Ter'), 'Ter')
        self.assertEqual(self.p.parse_aat3('Ala'), 'Ala')
        self.assertEqual(self.p.parse_aat13('*'), '*')
        self.assertEqual(self.p.parse_aat13('A'), 'A')
        self.assertEqual(self.p.parse_aat13('Ter'), 'Ter')
        self.assertEqual(self.p.parse_aat13('Ala'), 'Ala')

    def test_parser_aa_seq(self):
        # mixto
        self.assertEqual(self.p.parse_aat1_seq('*'), '*')
        self.assertEqual(self.p.parse_aat1_seq('A'), 'A')
        self.assertEqual(self.p.parse_aat3_seq('Ter'), 'Ter')
        self.assertEqual(self.p.parse_aat3_seq('Ala'), 'Ala')
        self.assertEqual(self.p.parse_aat13_seq('*'), '*')
        self.assertEqual(self.p.parse_aat13_seq('A'), 'A')
        self.assertEqual(self.p.parse_aat13_seq('Ter'), 'Ter')
        self.assertEqual(self.p.parse_aat13_seq('Ala'), 'Ala')


if __name__ == '__main__':
    unittest.main()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
