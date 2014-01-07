import unittest

import bdi.sources.uta0

import hgvs.parser
import hgvs.hgvsmapper

class Test_HGVSMapper(unittest.TestCase):
    def setUp(self):
        self.bdi = bdi.sources.uta0.connect()
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()

    def test_1(self):
        # From garcia.tsv:
        # AOAH    NM_001177507.1:c.1486G>A      
        hgvs_g = 'NC_000007.13:g.36561662C>T'
        hgvs_c = 'NM_001637.3:c.1582G>A'
    
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.hm.hgvsg_to_hgvsc( var_g, 'NM_001637.3' )

        self.assertEqual( str(var_c) , hgvs_c )


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
