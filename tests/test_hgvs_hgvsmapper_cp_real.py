#
# Real data - cp tests
#
import csv
import os
import re
import unittest

import bdi.sources.uta0

import hgvs.hgvsmapper
import hgvs.parser

def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, 'r'), delimiter='\t')
    for rec in rdr:
        if rec['id'].startswith('#'):
            continue
        yield rec


class TestHgvsCToPReal(unittest.TestCase):

    def setUp(self):
        self.bdi = bdi.sources.uta0.connect()
        self._hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)
        self._hp = hgvs.parser.Parser()
        self._failed = []

    def test_hgvsc_to_hgvsp_real(self):
        infilename = 'real_gcp.tsv'
        outfilename = 'real_gcp.out'
        infile = os.path.join(os.path.dirname(__file__), 'data', infilename)
        outfile = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        self._run_cp_test(infile, outfile)

    def test_hgvsc_to_hgvsp_ext(self):
        infilename = 'ext.tsv'
        outfilename = 'ext.out'
        infile = os.path.join(os.path.dirname(__file__), 'data', infilename)
        outfile = os.path.join(os.path.dirname(__file__), 'data', outfilename)
        self._run_cp_test(infile, outfile)

    #
    # internal methods
    #

    def _run_cp_test(self, infile, outfile):

        with open(outfile, 'w') as out:
            out.write("id\tHGVSg\tHGVSc\tHGVSp\tConverterResult\tError\n")
            self._dup_regex = re.compile(r'dup[0-9]+$')
            for rec in gcp_file_reader(infile):
                self._run_comparison(rec, out)
        msg = "# failed: {}".format(len(self._failed))
        self.assertTrue(len(self._failed) == 0, msg)

    def _run_comparison(self, rec, out):
        (row_id, hgvsg, hgvsc, hgvsp_expected) = (rec['id'], rec['HGVSg'], rec['HGVSc'], rec['HGVSp'])

        hgvsc = self._dup_regex.sub('dup', hgvsc) # cleanup dupN

        if not row_id.startswith("#") and hgvsc and hgvsp_expected:
            try:
                var_c = self._hp.parse_hgvs_variant(hgvsc)
                var_p = self._hm.hgvsc_to_hgvsp(var_c,  hgvsp_expected.split(':')[0]) # hack until p.?, p.= etc parse
                hgvsp_actual = str(var_p)

                if hgvsp_expected != hgvsp_actual:
                    if hgvsp_expected.replace("*", "Ter") == hgvsp_actual:
                        pass
                    else:
                        self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, "MISMATCH")
            except Exception as e:
                self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, "NO_OUTPUT_EXCEPTION", e.message)


    def _append_fail(self, out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg):
        self._failed.append((row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg))
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg))

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
