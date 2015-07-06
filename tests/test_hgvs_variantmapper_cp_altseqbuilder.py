# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import unittest

from Bio.Seq import Seq
import recordtype

import hgvs.parser
import hgvs.utils.altseqbuilder as altseqbuilder

import framework.mock_input_source as mock_input_data_source


class TestAltSeqBuilder(unittest.TestCase):

    # root sequence = ""
    fn = os.path.join(os.path.dirname(__file__), 'data', 'sanity_cp.tsv')
    _datasource = mock_input_data_source.MockInputSource(fn)
    _parser = hgvs.parser.Parser()

    def test_substitution_start(self):
        hgvsc = "NM_999999.1:c.1A>T"
        expected_sequence = "AAAATCAAATTGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_substitution_middle(self):
        hgvsc = "NM_999999.1:c.6A>T"
        expected_sequence = "AAAATCAAAATGAATGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_substitution_end(self):
        hgvsc = "NM_999999.1:c.30G>C"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATACGGG"
        self._run_comparison(hgvsc, expected_sequence)

    # TODO - build in support when system can handle variants in 5'utr region
    # def test_insertion_before_start(self):
    #     hgvsc = "NM_999999.1:c.-1_1insGGG"
    #     expected_sequence = "AAAATCAAAGGGATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
    #     self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_start(self):
        hgvsc = "NM_999999.1:c.1_2insAAA"
        expected_sequence = "AAAATCAAAAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_middle(self):
        hgvsc = "NM_999999.1:c.22_23insT"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGTCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_insertion_end(self):
        hgvsc = "NM_999999.1:c.29_30insGG"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    # TODO - build in support when system can handle variants in 3'utr region
    # def test_insertion_after_end(self):
    #     hgvsc = "NM_999999.1:c.30_*1insAA"
    #     expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGAAGGGN"
    #     self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_start(self):
        hgvsc = "NM_999999.1:c.1del"
        expected_sequence = "AAAATCAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_middle(self):
        hgvsc = "NM_999999.1:c.2_7del"
        expected_sequence = "AAAATCAAAACGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_deletion_end(self):
        hgvsc = "NM_999999.1:c.30del"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_start(self):
        hgvsc = "NM_999999.1:c.1delinsTTTT"
        expected_sequence = "AAAATCAAATTTTTGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_middle(self):
        hgvsc = "NM_999999.1:c.2_3delinsAA"
        expected_sequence = "AAAATCAAAAAAAAAGCGAAAGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delins_end(self):
        hgvsc = "NM_999999.1:c.30delinsCCCC"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATACCCCGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_dup(self):
        hgvsc = "NM_999999.1:c.16_24dup"
        expected_sequence = "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGTTTCGCGCGAAATAGGGG"
        self._run_comparison(hgvsc, expected_sequence)

    def test_delete_gene(self):
        hgvsc = "NM_999999.1:c.-3_*1del"
        expected_sequence = ""
        self._run_comparison(hgvsc, expected_sequence)

    # def test_2_substitutions(self):
    #     pass
    #
    # def test_2_indel_no_net_frameshift(self):
    #     pass
    #
    # def test_2_indel_net_frameshift(self):
    #     pass

    def _run_comparison(self, hgvsc, expected_sequence):

        # test replicates the internal class of p_to_c
        class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                                      ['transcript_sequence', 'aa_sequence', 'cds_start', 'cds_stop',
                                                       'protein_accession'])):
            @classmethod
            def setup_transcript_data(cls, ac, ac_p, db, ref='GRCh37.p10'):
                """helper for generating RefTranscriptData from for c_to_p"""
                tx_info = db.get_tx_info(ac)
                tx_seq = db.get_tx_seq(ac)

                if tx_info is None or tx_seq is None:
                    raise hgvs.exceptions.HGVSError("Missing transcript data for accession: {}".format(ac))

                # use 1-based hgvs coords
                cds_start = tx_info['cds_start_i'] + 1
                cds_stop = tx_info['cds_end_i']

                # padding list so biopython won't complain during the conversion
                tx_seq_to_translate = tx_seq[cds_start - 1:cds_stop]
                if len(tx_seq_to_translate) % 3 != 0:
                    ''.join(list(tx_seq_to_translate).extend(['N'] * ((3 - len(tx_seq_to_translate) % 3) % 3)))

                tx_seq_cds = Seq(tx_seq_to_translate)
                protein_seq = str(tx_seq_cds.translate())

                transcript_data = RefTranscriptData(tx_seq, protein_seq, cds_start, cds_stop, ac_p)

                return transcript_data

        ac_p = "DUMMY"
        var = self._parser.parse_hgvs_variant(hgvsc)
        transcript_data = RefTranscriptData.setup_transcript_data(var.ac, ac_p, self._datasource)

        builder = altseqbuilder.AltSeqBuilder(var, transcript_data)
        insert_result = builder.build_altseq()
        actual_sequence = insert_result[0].transcript_sequence
        msg = "expected: {}\nactual  : {}".format(expected_sequence, actual_sequence)
        self.assertEqual(expected_sequence, actual_sequence, msg)


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
