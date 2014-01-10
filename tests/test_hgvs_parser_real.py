import csv
import os
import re
import unittest

import hgvs.parser


def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter='\t')
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


class TestParseReal(unittest.TestCase):
    longMessage = True
    
    def setUp(self):
        self.parser = hgvs.parser.Parser()
        self._failed = []

        infilename = 'real_gcp.tsv'
        self._infile = os.path.join(os.path.dirname(__file__), 'data', infilename)

    def test_parse_hgvsg(self):
        key = 'HGVSg'
        omissions = [r'\(\d+\)$', r'copy\d+$', r'dup\d+$']
        self._parse_all(key, omissions)

    def test_parse_hgvsc(self):
        key = 'HGVSc'
        omissions = [r'copy\d+$', r'dup\d+$']
        self._parse_all(key, omissions)

    def test_parse_hgvsp(self):
        key = 'HGVSp'
        omissions = [r'dup', r'ext', r'Met1\?']
        self._parse_all(key, omissions)

    #
    # internal methods
    #

    def _parse_all(self, key, omissions = []):
        """Parse all tags of a given type

        :param key: record key corresponding to hgvs tag type (g, c ot p)
        :type key: str
        :param omissions: list of raw string literals corresponding to unsupported types
        :type omissions: list of raw string literals
        """

        # skip anything unsupported
        regexes = [re.compile(x) for x in omissions]

        outfilename = 'gcp_real_parsed_{}.out'.format(key)
        outfile = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        with open(outfile, 'w') as out:
            out.write("id\t{}\n".format(key))
            for rec in gcp_file_reader(self._infile):
                if rec[key]:
                    hits = [True for x in regexes if len(x.findall(rec[key])) > 0]
                    if len(hits) == 0:
                        try:
                            var_x = self.parser.parse_hgvs_variant(rec[key])
                        except:
                            self._append_fail(out, rec['id'], rec[key])
            msg = "# failed: {}".format(len(self._failed))
            self.assertTrue(len(self._failed) == 0, msg)


    def _append_fail(self, out, row_id, hgvs_tag):
        self._failed.append((row_id, hgvs_tag))
        out.write("{}\t{}\n".format(row_id, hgvs_tag))

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
