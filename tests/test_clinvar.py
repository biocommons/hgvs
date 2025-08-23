import csv
import gzip
import io
import os
import unittest
from pathlib import Path

import pytest

import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.sequencevariant
import hgvs.variantmapper
from support import CACHE
from support.crosschecker import CrossChecker, LineIterator

data_path = Path(__file__).parent / "data" / "clinvar.gz"


class Test_Clinvar(unittest.TestCase, CrossChecker):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls.hp = hgvs.parser.Parser()

    @pytest.mark.extra
    def test_clinvar(self, fn: Path = data_path, mod=None):
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
        fh = LineIterator(
            fh=(
                io.TextIOWrapper(gzip.open(fn), encoding="utf-8")  # noqa: SIM115
                if fn.name.endswith(".gz")
                else fn.open()
            ),
            skip_comments=True,
        )
        for rec in csv.DictReader(fh, delimiter="\t"):
            if mod and fh.lines_read % mod != 0:
                continue
            print(rec["gene"])  # noqa: T201
            variants = (self.hp.parse_hgvs_variant(vs) for vs in rec["hgvs_variants"].split())
            variants = [v for v in variants if v.type in "cg" and not v.ac.startswith("NG")]
            try:
                msg = self.crosscheck_variant_group(list(variants))
            except Exception as e:
                e.args = (e.args[0] + f"\n  at line {fh.lines_read}: {fh.last_line}",)
                raise
            self.assertIsNone(msg, lambda: msg + f"\n   on line {fh.lines_read}: {fh.last_line}")  # noqa: B023


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
