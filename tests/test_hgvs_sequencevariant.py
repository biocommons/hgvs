# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

from nose.plugins.attrib import attr

import hgvs
import hgvs.sequencevariant
import hgvs.parser
import hgvs.dataproviders.uta


@attr(tags=["quick", "models"])
class Test_SequenceVariant(unittest.TestCase):
    def test_SequenceVariant(self):
        var = hgvs.sequencevariant.SequenceVariant(ac="AC", type="B", posedit="1234DE>FG")
        self.assertEqual(str(var), "AC:B.1234DE>FG")
    
    def test_fill_ref(self):
        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE","run"), cache="tests/data/cache.hdp")
        
        # fill reference for sequence variants
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.31_32delTT")
        
        var = hp.parse_hgvs_variant("NM_001166478.1:c.31_32del2").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.31_32delTT")
        
        var = hp.parse_hgvs_variant("NM_001166478.1:c.2_7delinsTTTAGA").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.2_7delTGAAGAinsTTTAGA")
        
        var = hp.parse_hgvs_variant("NM_001166478.1:c.35_36dup").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.35_36dupTC")
        
        var = hp.parse_hgvs_variant("NM_001166478.1:c.18_19insACT").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.18_19insACT")

        var = hp.parse_hgvs_variant("NM_001166478.1:c.31=").fill_ref(hdp)
        self.assertEqual(str(var), "NM_001166478.1:c.31T=")

    def test_format(self):
        hp = hgvs.parser.Parser()
        
        # Global default settings
        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Arg")
        self.assertEqual(str(var), "NP_001628.1:p.Gly528Arg")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Arg")
        
        # Change global settings
        hgvs.global_config.formatting.p_3_letter = False
        self.assertEqual(str(var), "NP_001628.1:p.G528R")
        
        # Custom settings
        hgvs.global_config.formatting.p_3_letter = True
        conf = {"p_3_letter" : False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528R")
        
        var = hp.parse_hgvs_variant("NP_001628.1:p.Gly528Ter")
        conf = {"p_term_asterisk" : True}
        self.assertEqual(var.format(conf), "NP_001628.1:p.Gly528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")
        conf = {"p_3_letter" : False}
        self.assertEqual(var.format(conf), "NP_001628.1:p.G528*")
        self.assertEqual(var.format(), "NP_001628.1:p.Gly528Ter")
        


if __name__ == "__main__":
    unittest.main()

# <LICENSE>
# Copyright 2013-2015 HGVS Contributors (https://github.com/biocommons/hgvs)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
