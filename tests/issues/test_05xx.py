import os
import unittest

import pytest

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.validator
import hgvs.variantmapper
from support import CACHE


@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp, replace_reference=False)
        cls.vm_rr = hgvs.variantmapper.VariantMapper(cls.hdp, replace_reference=True)
        cls.hp = hgvs.parser.Parser()
        cls.hn = hgvs.normalizer.Normalizer(cls.hdp)
        cls.hv = hgvs.validator.IntrinsicValidator()
        cls.am37 = hgvs.assemblymapper.AssemblyMapper(
            cls.hdp, replace_reference=True, assembly_name="GRCh37", alt_aln_method="splign"
        )
        cls.am38 = hgvs.assemblymapper.AssemblyMapper(
            cls.hdp, replace_reference=True, assembly_name="GRCh38", alt_aln_method="splign"
        )
        cls.vn = hgvs.normalizer.Normalizer(cls.hdp, shuffle_direction=3, cross_boundaries=True)

    def test_525(self):
        """https://github.com/biocommons/hgvs/issues/525"""
        # simple test case
        hgvs = "NM_001637.3:c.3_4insTAG"  # insert stop in phase at AA 2
        var_c = self.hp.parse_hgvs_variant(hgvs)
        var_p = self.am38.c_to_p(var_c)
        assert str(var_p) == "NP_001628.1:p.(Gln2Ter)"
