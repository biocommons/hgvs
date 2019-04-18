# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import csv
import fileinput
import gzip
import io
import os
import sys
import types
import unittest

import pytest

import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.sequencevariant
import hgvs.variantmapper
from hgvs.exceptions import HGVSError

from support.crosschecker import CrossChecker, LineIterator
from support import CACHE

data_fn = os.path.join(os.path.dirname(__file__), "data", "clinvar.gz")


class Test_Clinvar(unittest.TestCase, CrossChecker):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp)
        self.hp = hgvs.parser.Parser()

    @pytest.mark.extra
    def test_clinvar(self, fn=data_fn, mod=None):
        """Test genome-transcript projections for 7498 clinvar variants in 4676 against genes
        for both GRCh37 and GRCh38 (when both are available).

        Approximate timings in various configurations:
        uta         sequences        time
        local       remote          1 min
        local       remote         40 min
        remote      local          10 min
        remote      remote         50 min

        In other words, you really want to run this with local sequences and UTA.
        """
        if sys.version_info < (3, ):
            fh = LineIterator(
                fh=gzip.open(fn) if fn.endswith(".gz") else io.open(fn), skip_comments=True)
        else:
            fh = LineIterator(
                fh=(io.TextIOWrapper(gzip.open(fn), encoding='utf-8')
                    if fn.endswith(".gz") else io.open(fn)),
                skip_comments=True)
        for rec in csv.DictReader(fh, delimiter=str('\t')):
            if mod and fh.lines_read % mod != 0:
                continue
            print(rec["gene"])
            variants = (self.hp.parse_hgvs_variant(vs) for vs in rec["hgvs_variants"].split())
            variants = [v for v in variants if v.type in "cg" and not v.ac.startswith("NG")]
            try:
                msg = self.crosscheck_variant_group(list(variants))
            except Exception as e:
                e.args = (e.args[0] + "\n  at line {fh.lines_read}: {fh.last_line}".format(fh=fh), )
                raise e
            self.assertIsNone(
                msg, lambda: msg + "\n   on line {fh.lines_read}: {fh.last_line}".format(fh=fh))


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
