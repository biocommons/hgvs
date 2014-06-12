# tests running from a local sqlite instance

import csv
import os
import unittest

import hgvs.dataproviders.uta

import hgvs.hgvsmapper as hgvsmapper
import hgvs.parser


def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter='\t')
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


class TestHGVSMapperFast(unittest.TestCase):
    def setUp(self):
        sqlite_path = 'sqlite://{}'.format(os.path.join(os.path.dirname(__file__), 'db', 'uta-test-1.db'))
        self.bdi = hgvs.dataproviders.uta.connect(sqlite_path)
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self.hp = hgvs.parser.Parser()
        self.failed = []

    def test_fast(self):
        for rec in gcp_file_reader('tests/data/sqlite_test_gcp.tsv'):
            self._test_gcp_mapping(rec)

    def _test_gcp_mapping(self, rec):
        var_g = self.hp.parse_hgvs_variant(rec['HGVSg'])
        var_c = self.hp.parse_hgvs_variant(rec['HGVSc'])
        var_p = self.hp.parse_hgvs_variant(rec['HGVSp']) if rec['HGVSp'] is not None and rec['HGVSp'] != '' else None

        # g -> c
        var_c_test = self.hm.hgvsg_to_hgvsc(var_g, var_c.ac)
        self.assertEquals(str(var_c_test), str(var_c))

        # c -> g
        var_g_test = self.hm.hgvsc_to_hgvsg(var_c, var_g.ac)
        self.assertEquals(str(var_g_test), str(var_g))

        if var_p is not None:
            # c -> p
            var_p_test = self.hm.hgvsc_to_hgvsp(var_c, var_p.ac)
            hgvs_p_test = str(var_p_test)
            self.assertEquals(hgvs_p_test, str(var_p))


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
