# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import pprint
import unittest

from nose.plugins.attrib import attr

from hgvs.exceptions import HGVSParseError
import hgvs.parser


class Test_Position(unittest.TestCase):
    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.parser = hgvs.parser.Parser()

    def test_parser_gauntlet(self):
        fn = os.path.join(os.path.dirname(__file__), 'data', 'gauntlet')
        for var in open(fn, 'r'):
            var = var.strip()
            if var.startswith('#') or var == '':
                continue
            v = self.parser.parse_hgvs_variant(var)
            self.assertEqual(var, str(v), 'parse-format roundtrip failed:' + pprint.pformat(v.posedit))

    @attr(tags=["quick"])
    def test_parser_reject(self):
        fn = os.path.join(os.path.dirname(__file__), 'data', 'reject')
        for var in open(fn, 'r'):
            var, msg = var.strip().split('\t')
            if var.startswith('#') or var == '':
                continue
            with self.assertRaises(HGVSParseError):
                self.parser.parse_hgvs_variant(var)
                self.assertTrue(False, msg='expected HGVSParseError: %s (%s)' % (var, msg))

    @attr(tags=["quick"])
    def test_parser_posedit_special(self):
        # See note in grammar about parsing p.=, p.?, and p.0
        self.assertEqual(str(self.parser.parse_p_posedit('0')), '0')
        self.assertEqual(str(self.parser.parse_p_posedit('0?')), '0?')
        #self.assertEqual( str(self.parser.parse_p_posedit('(0)')), '0?' )

        self.assertEqual(str(self.parser.parse_p_posedit('?')), '?')
        #self.assertEqual( str(self.parser.parse_p_posedit('??')), '?' )
        #self.assertEqual( str(self.parser.parse_p_posedit('(?)')), '?' )

        self.assertEqual(str(self.parser.parse_p_posedit('=')), '=')
        #self.assertEqual( str(self.parser.parse_p_posedit('=?')), '(=)' )
        self.assertEqual(str(self.parser.parse_p_posedit('(=)')), '(=)')


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
