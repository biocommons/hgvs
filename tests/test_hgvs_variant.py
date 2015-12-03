# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.variant
import hgvs.parser
import hgvs.dataproviders.uta


@attr(tags=["quick", "models"])
class Test_SequenceVariant(unittest.TestCase):
    def test_SequenceVariant(self):
        var = hgvs.variant.SequenceVariant(ac='AC', type='B', posedit='1234DE>FG')
        self.assertEqual(str(var), 'AC:B.1234DE>FG')
    
    def test_fill_ref(self):
        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect()
        
        # fill reference for sequence variants
        var = hp.parse_hgvs_variant('NM_001166478.1:c.31_32del').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.31_32delTT')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.31_32del2').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.31_32delTT')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.2_7delinsTTTAGA').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.2_7delTGAAGAinsTTTAGA')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.35_36dup').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.35_36dupTC')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.18_24inv').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.18_24invTCTCTTT')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.18_24inv7').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.18_24invTCTCTTT')
        
        var = hp.parse_hgvs_variant('NM_001166478.1:c.18_19insACT').fill_ref(hdp)
        self.assertEqual(str(var), 'NM_001166478.1:c.18_19insACT')


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
