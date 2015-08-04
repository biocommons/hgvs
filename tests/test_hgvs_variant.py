# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.variant


@attr(tags=["quick", "models"])
class Test_SequenceVariant(unittest.TestCase):
    def test_SequenceVariant(self):
        var = hgvs.variant.SequenceVariant(ac='AC', type='B', posedit='1234DE>FG')
        self.assertEqual(str(var), 'AC:B.1234DE>FG')
    
    def test_CompoundVariant(self):
        var1 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12A>B')
        var2 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='23C>D')
        var3 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='34E>F')
        var4 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='45G>H')
        var5 = hgvs.variant.SequenceVariant(ac='BC', type='c', posedit='56I>J')
        
        vars = [var1, var2, var3, var4]
        phas = [hgvs.variant.PHASE1, hgvs.variant.PHASE1, hgvs.variant.PHASE2, hgvs.variant.PHASE2]
        
        var = hgvs.variant.CompoundVariant(vars, phas)
        self.assertEqual(str(var), 'AC:c.[12A>B;23C>D];[34E>F;45G>H]')
        
        varn = hgvs.variant.CompoundVariant(vars, phas, uncertain=True)
        self.assertEqual(str(varn), 'AC:c.[12A>B(;)23C>D];[34E>F(;)45G>H]')
        
        var = hgvs.variant.CompoundVariant(vars, None)
        self.assertEqual(str(var), 'AC:c.[12A>B;23C>D;34E>F;45G>H]')
        
        vars = [var1, var2, var3, var5]
        
        var = hgvs.variant.CompoundVariant(vars, phas)
        self.assertEqual(str(var), '[AC:c.12A>B;AC:c.23C>D];[AC:c.34E>F;BC:c.56I>J]')
    
    def test_MosaicVariant(self):
        var1 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12=')
        var2 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12A>C')
        var3 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12A>G')
        
        vars = [var1, var2, var3]
        
        var = hgvs.variant.MosaicVariant(vars)
        self.assertEqual(str(var), 'AC:c.[12=/12A>C/12A>G]')
        
    
    def test_ChimericVariant(self):
        var1 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12=')
        var2 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12A>C')
        var3 = hgvs.variant.SequenceVariant(ac='AC', type='c', posedit='12A>G')
        
        vars = [var1, var2, var3]
        
        var = hgvs.variant.ChimericVariant(vars)
        self.assertEqual(str(var), 'AC:c.[12=//12A>C//12A>G]')


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
