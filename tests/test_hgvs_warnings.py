# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper


@attr(tags=["quick"])
class Test_Warnings(unittest.TestCase):
    def setUp(self):
        self.primary_assembly='GRCh37'
        self.alt_aln_method='splign'
        self.hdp = hgvs.dataproviders.uta.connect()
        self.evm = hgvs.variantmapper.EasyVariantMapper(hdp=self.hdp,
                                                        primary_assembly=self.primary_assembly,
                                                        alt_aln_method=self.alt_aln_method)
        self.hp = hgvs.parser.Parser()
        
    def test_ReferenceIndelDiscrepancy(self):
        var_c = self.hp.parse_hgvs_variant('NM_001637.3:c.1582G>A')
        var_cg = self.evm.c_to_g(var_c)
        self.assertTrue(
            any(w for w in var_cg.warnings if w.__class__ == hgvs.warnings.ReferenceIndelDiscrepancy),
            "{var_c.ac}~{self.primary_assembly} ({self.alt_aln_method}): expected a ReferenceIndelDiscrepancy".format(
                var_c=var_c,self=self))




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
