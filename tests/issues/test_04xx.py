import os
import unittest

import pytest

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.sequencevariant
import hgvs.validator
import hgvs.variantmapper
from hgvs.exceptions import HGVSDataNotAvailableError
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

    def test_424_430_nochange_parse_and_format(self):
        h = "NM_012.3:c.1="
        v = self.hp.parse_hgvs_variant(h)
        self.assertEqual("NM_012.3:c.1=", str(v))
        self.assertEqual("", v.posedit.edit.ref)
        self.assertEqual("", v.posedit.edit.alt)

        h = "NM_012.3:c.1A="
        v = self.hp.parse_hgvs_variant(h)
        self.assertEqual("NM_012.3:c.1=", str(v))
        self.assertEqual("A", v.posedit.edit.ref)
        self.assertEqual("A", v.posedit.edit.alt)

    def test_459_exception_when_ac_nonexistent(self):
        bogus_ac = "NM_000000.99"
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.seqfetcher.fetch_seq(bogus_ac)

    def test_499_whole_gene_dup(self):
        # Verify that 1_*1dup works
        self.am37.c_to_p(self.hp.parse_hgvs_variant("NM_001637.3:c.1_*1dup"))

        # Now try -1_*1dup (essentially, this is issue #499)
        self.am37.c_to_p(self.hp.parse_hgvs_variant("NM_001637.3:c.-1_*1dup"))
