import unittest

import hgvs.dataproviders.uta

import hgvs.parser
import hgvs.hgvsmapper

class Test_HGVSMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.hdp, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()

    def test_1(self):
        # From garcia.tsv:
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'
        hgvs_p = 'NP_001628.1:p.(Gly528Arg)'  # from Mutalyzer

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hm.g_to_c( var_g, 'NM_001637.3' )
        var_p = self.hm.c_to_p( var_c )

        self.assertEqual( str(var_c) , hgvs_c )
        self.assertEqual( str(var_p) , hgvs_p )


    def test_gcrp_invalid_input_type(self):
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hp.parse_hgvs_variant(hgvs_c)

        cases = {'gc': (self.hm.g_to_c, (var_c, 'NM_001637.3')),
                 'gr': (self.hm.g_to_r, (var_c, 'NM_001637.3')),
                 'rg': (self.hm.r_to_g, (var_c, 'NM_001637.3')),
                 'cg': (self.hm.c_to_g, (var_g, 'NM_001637.3')),
                 'cr': (self.hm.c_to_r, (var_g,)),
                 'rc': (self.hm.r_to_c, (var_g,)),
                 'cp': (self.hm.c_to_p, (var_g, None)),
                 }

        failures = []
        for key in cases:
            try:
                func, args = cases[key]
                var_result = func(*args)
                failures.append(key)
            except hgvs.exceptions.InvalidHGVSVariantError:
                pass

        self.assertFalse(failures, "conversions not failing: {}".format(failures))

    def test_gc_invalid_input_nm_accession(self):
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        with self.assertRaises(hgvs.exceptions.HGVSError):
            var_p = self.hm.c_to_p(var_g, 'NM_999999.1')




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
