# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import os

import unittest

import pytest

from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSParseError, HGVSInvalidVariantError
from hgvs.enums import Datum
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.sequencevariant
import hgvs.alignmentmapper
import hgvs.validator
import hgvs.variantmapper
from support import CACHE


@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False)
        self.vm_rr = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=True)
        self.hp = hgvs.parser.Parser()
        self.hn = hgvs.normalizer.Normalizer(self.hdp)
        self.hv = hgvs.validator.IntrinsicValidator()
        self.am37 = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh37', alt_aln_method='splign')
        self.am38 = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh38', alt_aln_method='splign')
        self.vn = hgvs.normalizer.Normalizer(self.hdp, shuffle_direction=3, cross_boundaries=True)

    def test_307_validator_rejects_valid_interval(self):
        # https://github.com/biocommons/hgvs/issues/307/
        # before fix, raises error; after fix, should pass
        self.hv.validate(self.hp.parse_hgvs_variant("NM_002858.3:c.1903-573_*1108del"))

    def test_314_parsing_identity_variant(self):
        v = self.hp.parse_hgvs_variant("NM_206933.2:c.6317=")
        self.assertEqual(str(v), "NM_206933.2:c.6317=")

        v = self.hp.parse_hgvs_variant("NM_206933.2:c.6317C=")
        self.assertEqual(str(v), "NM_206933.2:c.6317=")

    def test_316_provide_ttog_and_gtot_methods(self):
        h_c_g, h_c_t = "NC_000023.10:g.152864509_152864510insA", "NM_152274.3:c.21_22insT"
        h_nc_g, h_nc_t = "NC_000001.10:g.12776161G>A", "NR_111984.1:n.44G>A"

        v_c_g = self.hp.parse_hgvs_variant(h_c_g)
        v_c_t = self.hp.parse_hgvs_variant(h_c_t)
        v_nc_g = self.hp.parse_hgvs_variant(h_nc_g)
        v_nc_t = self.hp.parse_hgvs_variant(h_nc_t)

        self.assertEqual(h_c_t, str(self.am37.g_to_t(v_c_g, tx_ac=v_c_t.ac)))
        self.assertEqual(h_c_g, str(self.am37.t_to_g(v_c_t)))
        self.assertEqual(h_nc_t, str(self.am37.g_to_t(v_nc_g, tx_ac=v_nc_t.ac)))
        self.assertEqual(h_nc_g, str(self.am37.t_to_g(v_nc_t)))

    def test_317_improve_p_repr_of_syn_variants(self):
        # from original issue:
        v317 = self.hp.parse_hgvs_variant("NM_000059.3:c.7791A>G")
        self.assertEqual(str("NP_000050.2:p.(Lys2597=)"), str(self.am37.c_to_p(v317)))

    #def test_370_handle_multicodon_syn_variants(self):
    #    # Verify behavior of syn MNVs
    #    # NM_000059.3:
    #    #   :   3   6   9  12  15  18  21  24  27  30  33  36  39  42
    #    # 1 : ATG CCT ATT GGA TCC AAA GAG AGG CCA ACA TTT TTT GAA ATT
    #    #       M   P   I   G   W   F   E   R   P   T
    #    #                               ^      v21 (snv)
    #    #                                 ^    v22 (snv)
    #    #                                 ^^^  v2224 (mnv, 1 codon)
    #    #                               ^ ^    v2122 (mnv, 2 codons)
    #    # Glu (E): GA[AG]
    #    # Arg (R): AG[AG],  CG[UCAG]
    #    #
    #    # All of the following are synonymous changes:
    #    v21 = self.hp.parse_hgvs_variant("NM_000059.3:c.21G>A")  # GAG -> GAA (E)
    #    v22 = self.hp.parse_hgvs_variant("NM_000059.3:c.22A>C")  # AGG -> CGA (R)
    #    v2224 = self.hp.parse_hgvs_variant("NM_000059.3:c.22_24delinsCGA") # AGG -> CGA (R)
    #    v2122 = self.hp.parse_hgvs_variant("NM_000059.3:c.21_22delinsAC")  # GAGAGG -> GAACAG (ER)
    #
    #    self.assertEqual(b"NP_000050.2:p.(Glu7=)", str(self.am37.c_to_p(v21)))
    #    self.assertEqual(b"NP_000050.2:p.(Arg8=)", str(self.am37.c_to_p(v22)))
    #    self.assertEqual(b"NP_000050.2:p.(Arg8=)", str(self.am37.c_to_p(v2224)))
    #    self.assertEqual(b"NP_000050.2:p.(Glu7_Arg8=)", str(self.am37.c_to_p(v2122)))

    def test_322_raise_exception_when_mapping_bogus_variant(self):
        v = self.hp.parse_hgvs_variant("chrX:g.71684476delTGGAGinsAC")
        with self.assertRaises(HGVSInvalidVariantError):
            self.am37.g_to_c(v, "NM_018486.2")

    def test_324_error_normalizing_simple_inversion(self):
        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")    # no change

        vg = self.am37.c_to_g(v)
        v2 = self.am37.g_to_c(vg, tx_ac=v.ac)
        self.assertEqual(str(v), str(v2))    # no change after roundtrip

    def test_325_3p_UTR_treated_as_intronic(self):
        v = self.hp.parse_hgvs_variant("NM_000023.2:c.*6C>T")
        self.assertTrue(self.hv.validate(v))

    def test_326_handle_variants_in_par(self):
        "For a variant in the PAR, allow user user choose preference"
        hgvs_c = "NM_000451.3:c.584G>A"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)

        self.am37.in_par_assume = None
        with self.assertRaises(HGVSError):
            var_g = self.am37.c_to_g(var_c)

        self.am37.in_par_assume = 'X'
        self.assertEqual(self.am37.c_to_g(var_c).ac, "NC_000023.10")

        self.am37.in_par_assume = 'Y'
        self.assertEqual(self.am37.c_to_g(var_c).ac, "NC_000024.9")

        self.am37.in_par_assume = None

    def test_330_incorrect_end_datum_post_ter(self):
        # https://github.com/biocommons/hgvs/issues/330/
        # In a variant like NM_004006.2:c.*87_91del, the interval
        # start and end are parsed independently. The * binds to
        # start but not end, causing end to have an incorrect datum.
        v = self.hp.parse_hgvs_variant("NM_004006.2:c.*87_91del")
        self.assertEqual(v.posedit.pos.start.datum, Datum.CDS_END)
        self.assertEqual(v.posedit.pos.end.datum, Datum.CDS_END)

    def test_334_delins_normalization(self):
        # also tests 335 re: inv including sequence (e.g., NC_000009.11:g.36233991_36233992invCA)
        v = self.hp.parse_hgvs_variant("NC_000009.11:g.36233991_36233992delCAinsAC")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NC_000009.11:g.36233991_36233992delinsAC")

        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674delCCinsGG")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")

        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1_3delATGinsCAT")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1_3inv")

    def test_346_reject_partial_alignments(self):
        # hgvs-346: verify that alignment data covers full-length transcript
        with self.assertRaises(HGVSDataNotAvailableError):
            hgvs.alignmentmapper.AlignmentMapper(
                self.hdp, tx_ac="NM_001290223.1", alt_ac="NC_000010.10", alt_aln_method="splign")

    def test_349_incorrect_normalization_insGG_to_dup2(self):
        # With dupN support, NM_005877.4:c.1104_1105insGG would
        # normalize to NM_005877.4:c.1104_1105dup2 (which is
        # different).
        original_var = "NM_005877.4:c.1104_1105insGG"
        self.assertEqual(original_var, str(
            self.hn.normalize(self.hp.parse_c_variant(original_var))))

    def test_379_move_replace_reference_to_variantmapper(self):
        # replace_reference code was in am, not vm. That meant that using vm directly
        # resulted in variants that were not reference corrected.
        g_var = self.hp.parse_hgvs_variant("NC_000006.11:g.44275011T=")
        c_var = self.hp.parse_hgvs_variant(
            "NM_020745.3:c.1015G>A")    # correct projection with ref replacement
        self.assertEqual(c_var, self.am37.g_to_c(g_var, "NM_020745.3"))    # previously okay
        self.assertEqual(c_var, self.vm_rr.g_to_c(g_var, "NM_020745.3"))    # previously wrong

    def test_381_c_to_p_error_with_del_variants(self):
        hgvs_c = "NM_000302.3:c.1594_1596del"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        self.am37.c_to_p(var_c)    # raises exception before fixing

    def test_393_posedit_for_unknown_p_effect(self):
        hgvs_c = "NM_001330368.1:c.641-3353C>A"
        hgvs_p = "NP_001317297.1:p.?"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        translated_var_p = self.am37.c_to_p(var_c)
        parsed_var_p = self.hp.parse_hgvs_variant(hgvs_p)
        self.assertIsNone(parsed_var_p.posedit)
        self.assertIsNone(translated_var_p.posedit)
        self.assertEqual(hgvs_p, str(translated_var_p))
