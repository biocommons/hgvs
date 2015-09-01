# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.parser
import hgvs.formatter
import hgvs.dataproviders.uta


@attr(tags=["quick", "models"])
class Test_Formatter(unittest.TestCase):
    def test_formatter(self):
        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect()
        hf1  = hgvs.formatter.Formatter(hdp)
        hf2  = hgvs.formatter.Formatter(hdp, fill_ref=False)
        
        # fill reference for sequence variants
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.31_32del'))), 'NM_001166478.1:c.31_32delTT')
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.31_32del2'))), 'NM_001166478.1:c.31_32delTT')
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.2_7delinsTTTAGA'))), 'NM_001166478.1:c.2_7delTGAAGAinsTTTAGA')
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.35_36dup'))), 'NM_001166478.1:c.35_36dupTC')
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.18_24inv'))), 'NM_001166478.1:c.18_24invTCTCTTT')
        self.assertEqual(str(hf1.format(hp.parse_hgvs_variant('NM_001166478.1:c.18_24inv7'))), 'NM_001166478.1:c.18_24invTCTCTTT')
        
        # remove reference for sequence variants
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.31_32delTT'))), 'NM_001166478.1:c.31_32del')
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.31_32del2'))), 'NM_001166478.1:c.31_32del')
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.2_7delTGAAGAinsTTTAGA'))), 'NM_001166478.1:c.2_7delinsTTTAGA')
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.35_36dupTC'))), 'NM_001166478.1:c.35_36dup')
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.18_24invTCTCTTT'))), 'NM_001166478.1:c.18_24inv')
        self.assertEqual(str(hf2.format(hp.parse_hgvs_variant('NM_001166478.1:c.18_24inv7'))), 'NM_001166478.1:c.18_24inv')



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
