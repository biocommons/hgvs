# -*- coding: utf-8 -*-

"""Tests uta postgresql client"""


from __future__ import absolute_import, division, print_function, unicode_literals

import os
import re
import unittest

from hgvs.exceptions import HGVSDataNotAvailableError
import hgvs.dataproviders.uta
import hgvs.edit
import hgvs.location
import hgvs.posedit
import hgvs.variantmapper
import hgvs.variant


class UTA_Base(object):
    def test_get_acs_for_protein_seq(self):
        exp = ['NP_001005405.1', 'MD5_8fc09b1d9a38a8c55176a0fa922df227']
        s = """
        mgccgcsggc gsgcggcgsg sggcgsgcgg cgssccvpic cckpvcccvp acscsscgsc
        ggskggcgsc gsskggcgsc gcsqsncckp ccsssgcgsf ccqsscskpc ccqssccqss
        cckpcccqss ccqsscfkpc ccqssccvpv ccqcki
        """

        s = re.sub('\s+', '', s.upper())
        self.assertEqual(sorted(self.hdp.get_acs_for_protein_seq(s)), sorted(exp))

        exp = ['NP_071928.2', 'MD5_ffb0d4adbd5e0b5d71678228b3696984']
        s = """
        masetektha llqtcstesl isslglgafc lvadrllqfs tiqqndwlra lsdnavhcvi
        gmwswavvtg ikkktdfgei ilagflasvi dvdhfflags mslkaaltlp rrpflhcstv
        ipvvvltlkf tmhlfklkds wcflpwmlfi swtshhirdg irhglwicpf gktsplpfwl
        yviitsslph icsfvmyltg trqmmsskhg vridv
        """

        s = re.sub('\s+', '', s.upper())
        self.assertEqual(sorted(self.hdp.get_acs_for_protein_seq(s)), sorted(exp))

    def test_get_gene_info(self):
        gene_info = self.hdp.get_gene_info('VHL')
        self.assertEqual('VHL', gene_info['hgnc'])
        self.assertEqual('3p25.3', gene_info['maploc'])
        self.assertEqual(6, len(gene_info))

    def test_get_tx_exons(self):
        tx_exons = self.hdp.get_tx_exons('NM_000551.3', 'NC_000003.11', 'splign')
        self.assertEqual(3, len(tx_exons))

    def test_get_tx_exons_invalid_tx_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons('NM_999999.9', 'NC_000003.11', 'splign')

    def test_get_tx_exons_invalid_alt_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons('NM_000551.3', 'NC_000999.9', 'splign')

    def test_get_tx_exons_invalid_alt_aln_method(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons('NM_000551.3', 'NC_000999.9', 'best')

    def test_get_tx_for_gene(self):
        tig = self.hdp.get_tx_for_gene('VHL')
        self.assertEqual(12, len(tig))

    def test_get_tx_for_gene_invalid_gene(self):
        tig = self.hdp.get_tx_for_gene('GENE')
        self.assertEqual(0, len(tig))

    def test_get_tx_info(self):
        tx_info = self.hdp.get_tx_info('NM_000051.3', 'AC_000143.1', 'splign')
        self.assertEqual(385, tx_info['cds_start_i'])
        self.assertEqual(9556, tx_info['cds_end_i'])
        self.assertEqual('AC_000143.1', tx_info['alt_ac'])

    def test_get_tx_info_invalid_tx_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_info('NM_999999.9', 'AC_000143.1', 'splign')

    def test_get_tx_mapping_options(self):
        tx_mapping_options = self.hdp.get_tx_mapping_options('NM_000551.3')
        self.assertIn(['NM_000551.3', 'NC_000003.11', 'splign'], tx_mapping_options)
        self.assertIn(['NM_000551.3', 'NC_000003.11', 'blat'], tx_mapping_options)

    def test_get_tx_mapping_options_invalid(self):
        tx_info_options = self.hdp.get_tx_mapping_options('NM_999999.9')
        self.assertEqual(tx_info_options, [])

    def test_get_tx_seq(self):
        tsi = self.hdp._get_tx_seq('NM_000551.3')
        seq = tsi
        self.assertEqual(4560, len(seq))
        self.assertTrue(seq.startswith('cctcgcctccgttacaacgg'.upper()))
        self.assertTrue(seq.endswith('aaaggaaatag'.upper()))

    def test_get_tx_seq_invalid_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp._get_tx_seq('NM_999999.9')


class Test_hgvs_dataproviders_uta_UTA_default(unittest.TestCase, UTA_Base):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect()


class Test_hgvs_dataproviders_uta_UTA_default_with_pooling(unittest.TestCase, UTA_Base):
    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(pooling=True)


class TestUTACache(Test_hgvs_dataproviders_uta_UTA_default):
    def _create_cdna_variant(self):
        start = hgvs.location.SimplePosition(118898437)
        end = hgvs.location.SimplePosition(118898437)
        iv = hgvs.location.Interval(start=start, end=end)
        edit = hgvs.edit.NARefAlt(ref='G', alt='T')
        posedit = hgvs.posedit.PosEdit(pos=iv, edit=edit)
        genomic_variant = hgvs.variant.SequenceVariant(ac='NC_000011.9', type='g', posedit=posedit, )
        variantmapper = hgvs.variantmapper.VariantMapper(self.hdp)
        return variantmapper.g_to_c(genomic_variant, 'NM_001164277.1')

    def test_deterministic_cache_results(self):
        """
        Check that identical request to the UTA yields the same results.
        """
        var1 = self._create_cdna_variant()
        var2 = self._create_cdna_variant()
        self.assertEqual(str(var1), str(var2))


if __name__ == '__main__':
    unittest.main()

    ## <LICENSE>
    ## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
