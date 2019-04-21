# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

# Real data - cp tests

import os
import re
import unittest
from sys import version_info

if version_info < (3, ):
    import unicodecsv as csv
else:
    import csv

import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.parser
from support import CACHE


def gcp_file_reader(fn):
    rdr = csv.DictReader(open(fn, "r"), delimiter=str("\t"))
    for rec in rdr:
        if rec["id"].startswith("#"):
            continue
        yield rec

class TestHgvsCToPReal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "learn"), cache=CACHE)
        cls._hm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls._hp = hgvs.parser.Parser()
        cls._failed = []

    def test_c_to_p_ext(self):
        infilename = "ext.tsv"
        outfilename = "ext.out"
        infile = os.path.join(os.path.dirname(__file__), "data", infilename)
        outfile = os.path.join(os.path.dirname(__file__), "data", outfilename)
        self._run_cp_test(infile, outfile)

    #
    # internal methods
    #

    def _run_cp_test(self, infile, outfile):

        with open(outfile, "w") as out:
            out.write("id\tHGVSg\tHGVSc\tHGVSp\tConverterResult\tError\n")
            self._dup_regex = re.compile(r"dup[0-9]+$")
            for rec in gcp_file_reader(infile):
                self._run_comparison(rec, out)
        msg = "# failed: {}".format(len(self._failed))
        self.assertTrue(len(self._failed) == 0, msg)

    def _run_comparison(self, rec, out):
        (row_id, hgvsg, hgvsc, hgvsp_expected) = (rec["id"], rec["HGVSg"], rec["HGVSc"],
                                                  rec["HGVSp"])

        hgvsc = self._dup_regex.sub("dup", hgvsc)    # cleanup dupN

        if not row_id.startswith("#") and hgvsc and hgvsp_expected:
            try:
                var_c = self._hp.parse_hgvs_variant(hgvsc)
                var_p = self._hm.c_to_p(
                    var_c,
                    hgvsp_expected.split(":")[0])    # hack until p.?, p.= etc parse
                hgvsp_actual = str(var_p)

                if hgvsp_expected != hgvsp_actual:
                    if hgvsp_expected.replace("*", "Ter") == hgvsp_actual:
                        pass
                    else:
                        self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual,
                                          "MISMATCH")
            except Exception as e:
                self._append_fail(out, row_id, hgvsg, hgvsc, hgvsp_expected, "NO_OUTPUT_EXCEPTION",
                                  e.message)

    def _append_fail(self, out, row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg):
        self._failed.append((row_id, hgvsg, hgvsc, hgvsp_expected, hgvsp_actual, msg))
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(row_id, hgvsg, hgvsc, hgvsp_expected,
                                                    hgvsp_actual, msg))

    def test_c_to_p_format(self):
        hgvsc = 'NM_022464.4:c.3G>A'
        hgvsp_expected = 'NP_071909.1:p.Met1?'
        hgvsp_expected_alternative = 'NP_071909.1:p.?'

        var_c = self._hp.parse_hgvs_variant(hgvsc)
        var_p = self._hm.c_to_p(
            var_c, hgvsp_expected.split(":")[0])

        default_format_p = var_p.format()
        self.assertEqual(hgvsp_expected, default_format_p)

        alt_format_p = var_p.format(conf={'p_init_met': False})
        self.assertEqual(hgvsp_expected_alternative, alt_format_p)


if __name__ == "__main__":
    unittest.main()

# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
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
