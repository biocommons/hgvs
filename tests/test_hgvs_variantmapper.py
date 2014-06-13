import unittest

from nose.plugins.attrib import attr

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper


@attr(tags=["quick"])
class Test_VariantMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()

    def test_VariantMapper_quick(self):
        # From garcia.tsv:
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'
        hgvs_p = 'NP_001628.1:p.(Gly528Arg)'  # from Mutalyzer

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.vm.g_to_c( var_g, 'NM_001637.3' )
        var_p = self.vm.c_to_p( var_c )

        self.assertEqual( str(var_c) , hgvs_c )
        self.assertEqual( str(var_p) , hgvs_p )


    def test_gcrp_invalid_input_type(self):
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hp.parse_hgvs_variant(hgvs_c)

        cases = {'gc': (self.vm.g_to_c, (var_c, 'NM_001637.3')),
                 'gr': (self.vm.g_to_r, (var_c, 'NM_001637.3')),
                 'rg': (self.vm.r_to_g, (var_c, 'NM_001637.3')),
                 'cg': (self.vm.c_to_g, (var_g, 'NM_001637.3')),
                 'cr': (self.vm.c_to_r, (var_g,)),
                 'rc': (self.vm.r_to_c, (var_g,)),
                 'cp': (self.vm.c_to_p, (var_g, None)),
                 }

        failures = []
        for key in cases:
            try:
                func, args = cases[key]
                var_result = func(*args)
                failures.append(key)
            except hgvs.exceptions.HGVSInvalidVariantError:
                pass

        self.assertFalse(failures, "conversions not failing: {}".format(failures))

    def test_gc_invalid_input_nm_accession(self):
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        with self.assertRaises(hgvs.exceptions.HGVSError):
            var_p = self.vm.c_to_p(var_g, 'NM_999999.1')



@attr(tags=["quick"])
class Test_EasyVariantMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.evm = hgvs.variantmapper.EasyVariantMapper(self.hdp,primary_assembly='GRCh37', alt_aln_method='splign')
        self.hgvs = {
            'g': 'NC_000007.13:g.36561662C>T',
            'c': 'NM_001637.3:c.1582G>A',
            'r': 'NM_001637.3:r.1983G>A',
            'p': 'NP_001628.1:p.(Gly528Arg)'
        }
        hp = hgvs.parser.Parser()
        self.var = { k:hp.parse_hgvs_variant(v) for k,v in self.hgvs.iteritems() }


    def test_c_to_g(self):
        self.assertEqual(self.hgvs['g'], str(self.evm.c_to_g(self.var['c'])))

    def test_c_to_r(self):
        self.assertEqual(self.hgvs['r'], str(self.evm.c_to_r(self.var['c'])))

    def test_g_to_c(self):
        self.assertEqual(self.hgvs['c'], str(self.evm.g_to_c(self.var['g'], self.var['c'].ac)))

    def test_g_to_r(self):
        self.assertEqual(self.hgvs['r'], str(self.evm.g_to_r(self.var['g'], self.var['r'].ac)))

    def test_c_to_p(self):
        self.assertEqual(self.hgvs['p'], str(self.evm.c_to_p(self.var['c'])))



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
