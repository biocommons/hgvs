# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import unittest

from nose.plugins.attrib import attr

import hgvs.dataproviders.uta

import hgvs.location
import hgvs.parser
from hgvs.exceptions import HGVSError
from hgvs.transcriptmapper import TranscriptMapper


@attr(tags=["quick"])
class Test_transcriptmapper(unittest.TestCase):
    ref = 'GRCh37.p10'

    @classmethod
    def setUp(cls):
        cls.hdp = hgvs.dataproviders.uta.connect()

    def test_transcriptmapper_failures(self):
        self.assertRaises(HGVSError, TranscriptMapper, self.hdp,
                          tx_ac='bogus',
                          alt_ac='NM_033089.6',
                          alt_aln_method='splign')
        self.assertRaises(HGVSError, TranscriptMapper, self.hdp,
                          tx_ac='NM_033089.6',
                          alt_ac='bogus',
                          alt_aln_method='splign')
        self.assertRaises(HGVSError, TranscriptMapper, self.hdp,
                          tx_ac='NM_000051.3',
                          alt_ac='NC_000011.9',
                          alt_aln_method='bogus')

    def test_transcriptmapper_TranscriptMapper_LCE3C_uncertain(self):
        """Use NM_178434.2 tests to test mapping with uncertain positions"""
        tx_ac = 'NM_178434.2'
        alt_ac = 'NC_000001.10'
        tm = TranscriptMapper(self.hdp, tx_ac, alt_ac, alt_aln_method='splign')
        parser = hgvs.parser.Parser()
        test_cases = [
            {
                'g': parser.parse_g_interval('(152573138)'),
                'r': parser.parse_r_interval('(1)'),
                'c': parser.parse_c_interval('(-70)')
            },
            {
                'g': parser.parse_g_interval('(152573138_152573139)'),
                'r': parser.parse_r_interval('(1_2)'),
                'c': parser.parse_c_interval('(-70_-69)')
            },
            # ? is not yet supported
            # {'g': parser.parse_g_interval('(?_152573139)'), 'r': parser.parse_r_interval('(?_2)'), 'c': parser.parse_c_interval('(?_-69)')},
            # {'g': parser.parse_g_interval('(152573138_?)'), 'r': parser.parse_r_interval('(1_?)'), 'c': parser.parse_c_interval('(-70_?)')},
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_LCE3C(self):
        """NM_178434.2: LCE3C single exon, strand = +1, all coordinate input/output are in HGVS"""
        tx_ac = 'NM_178434.2'
        alt_ac = 'NC_000001.10'
        tm = TranscriptMapper(self.hdp, tx_ac, alt_ac, alt_aln_method='splign')
        parser = hgvs.parser.Parser()
        test_cases = [
            # 5'
            {
                'g': parser.parse_g_interval('152573138'),
                'r': parser.parse_r_interval('1'),
                'c': parser.parse_c_interval('-70')
            },
            {
                'g': parser.parse_g_interval('152573140'),
                'r': parser.parse_r_interval('3'),
                'c': parser.parse_c_interval('-68')
            },
            # cds
            {
                'g': parser.parse_g_interval('152573207'),
                'r': parser.parse_r_interval('70'),
                'c': parser.parse_c_interval('-1')
            },
            {
                'g': parser.parse_g_interval('152573208'),
                'r': parser.parse_r_interval('71'),
                'c': parser.parse_c_interval('1')
            },
            # 3'
            {
                'g': parser.parse_g_interval('152573492'),
                'r': parser.parse_r_interval('355'),
                'c': parser.parse_c_interval('285')
            },
            {
                'g': parser.parse_g_interval('152573493'),
                'r': parser.parse_r_interval('356'),
                'c': parser.parse_c_interval('*1')
            },
            {
                'g': parser.parse_g_interval('152573560'),
                'r': parser.parse_r_interval('423'),
                'c': parser.parse_c_interval('*68')
            },
            {
                'g': parser.parse_g_interval('152573562'),
                'r': parser.parse_r_interval('425'),
                'c': parser.parse_c_interval('*70')
            },
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_HIST3H2A(self):
        """NM_033445.2: LCE3C single exon, strand = -1, all coordinate input/output are in HGVS"""
        tx_ac = 'NM_033445.2'
        alt_ac = 'NC_000001.10'
        tm = TranscriptMapper(self.hdp, tx_ac, alt_ac, alt_aln_method='splign')
        parser = hgvs.parser.Parser()
        test_cases = [
            # 3'
            {
                'g': parser.parse_g_interval('228645560'),
                'r': parser.parse_r_interval('1'),
                'c': parser.parse_c_interval('-42')
            },
            {
                'g': parser.parse_g_interval('228645558'),
                'r': parser.parse_r_interval('3'),
                'c': parser.parse_c_interval('-40')
            },
            # cds
            {
                'g': parser.parse_g_interval('228645519'),
                'r': parser.parse_r_interval('42'),
                'c': parser.parse_c_interval('-1')
            },
            {
                'g': parser.parse_g_interval('228645518'),
                'r': parser.parse_r_interval('43'),
                'c': parser.parse_c_interval('1')
            },
            # 5'
            {
                'g': parser.parse_g_interval('228645126'),
                'r': parser.parse_r_interval('435'),
                'c': parser.parse_c_interval('393')
            },
            {
                'g': parser.parse_g_interval('228645125'),
                'r': parser.parse_r_interval('436'),
                'c': parser.parse_c_interval('*1')
            },
            {
                'g': parser.parse_g_interval('228645124'),
                'r': parser.parse_r_interval('437'),
                'c': parser.parse_c_interval('*2')
            },
            {
                'g': parser.parse_g_interval('228645065'),
                'r': parser.parse_r_interval('496'),
                'c': parser.parse_c_interval('*61')
            },
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_LCE2B(self):
        """NM_014357.4: LCE2B, two exons, strand = +1, all coordinate input/output are in HGVS"""
        tx_ac = 'NM_014357.4'
        alt_ac = 'NC_000001.10'
        tm = TranscriptMapper(self.hdp, tx_ac, alt_ac, alt_aln_method='splign')
        parser = hgvs.parser.Parser()
        test_cases = [
            # 5'
            {
                'g': parser.parse_g_interval('152658599'),
                'r': parser.parse_r_interval('1'),
                'c': parser.parse_c_interval('-54')
            },
            {
                'g': parser.parse_g_interval('152658601'),
                'r': parser.parse_r_interval('3'),
                'c': parser.parse_c_interval('-52')
            },
            # cds
            {
                'g': parser.parse_g_interval('152659319'),
                'r': parser.parse_r_interval('54'),
                'c': parser.parse_c_interval('-1')
            },
            {
                'g': parser.parse_g_interval('152659320'),
                'r': parser.parse_r_interval('55'),
                'c': parser.parse_c_interval('1')
            },
            # around end of exon 1
            {
                'g': parser.parse_g_interval('152658632'),
                'r': parser.parse_r_interval('34'),
                'c': parser.parse_c_interval('-21')
            },
            {
                'g': parser.parse_g_interval('152658633'),
                'r': parser.parse_r_interval('34+1'),
                'c': parser.parse_c_interval('-21+1')
            },
            # span
            {
                'g': parser.parse_g_interval('152658633_152659299'),
                'r': parser.parse_r_interval('34+1_35-1'),
                'c': parser.parse_c_interval('-21+1_-20-1')
            },
            # around beginning of exon 2
            {
                'g': parser.parse_g_interval('152659300'),
                'r': parser.parse_r_interval('35'),
                'c': parser.parse_c_interval('-20')
            },
            {
                'g': parser.parse_g_interval('152659299'),
                'r': parser.parse_r_interval('35-1'),
                'c': parser.parse_c_interval('-20-1')
            },
            # around end of exon 2
            {
                'g': parser.parse_g_interval('152659652'),
                'r': parser.parse_r_interval('387'),
                'c': parser.parse_c_interval('333')
            },
            {
                'g': parser.parse_g_interval('152659653'),
                'r': parser.parse_r_interval('388'),
                'c': parser.parse_c_interval('*1')
            },
            # span
            {
                'g': parser.parse_g_interval('152659651_152659654'),
                'r': parser.parse_r_interval('386_389'),
                'c': parser.parse_c_interval('332_*2')
            },
            # 3'
            {
                'g': parser.parse_g_interval('152659877'),
                'r': parser.parse_r_interval('612'),
                'c': parser.parse_c_interval('*225')
            },
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_PTH2(self):
        """NM_178449.3: PTH2, two exons, strand = -1, all coordinate input/output are in HGVS"""
        tx_ac = 'NM_178449.3'
        alt_ac = 'NC_000019.9'
        tm = TranscriptMapper(self.hdp, tx_ac, alt_ac, alt_aln_method='splign')
        parser = hgvs.parser.Parser()
        test_cases = [
            # 3'
            {
                'g': parser.parse_g_interval('49926698'),
                'r': parser.parse_r_interval('1'),
                'c': parser.parse_c_interval('-102')
            },
            # cds
            {
                'g': parser.parse_g_interval('49926597'),
                'r': parser.parse_r_interval('102'),
                'c': parser.parse_c_interval('-1')
            },
            {
                'g': parser.parse_g_interval('49926596'),
                'r': parser.parse_r_interval('103'),
                'c': parser.parse_c_interval('1')
            },
            # around end of exon 1
            {
                'g': parser.parse_g_interval('49926469'),
                'r': parser.parse_r_interval('230'),
                'c': parser.parse_c_interval('128')
            },
            {
                'g': parser.parse_g_interval('49926468'),
                'r': parser.parse_r_interval('230+1'),
                'c': parser.parse_c_interval('128+1')
            },
            # span
            {
                'g': parser.parse_g_interval('49925901_49926467'),
                'r': parser.parse_r_interval('230+2_231-2'),
                'c': parser.parse_c_interval('128+2_129-2')
            },
            # around beginning of exon 2
            {
                'g': parser.parse_g_interval('49925900'),
                'r': parser.parse_r_interval('231-1'),
                'c': parser.parse_c_interval('129-1')
            },
            {
                'g': parser.parse_g_interval('49925899'),
                'r': parser.parse_r_interval('231'),
                'c': parser.parse_c_interval('129')
            },
            # around end of exon 2
            {
                'g': parser.parse_g_interval('49925725'),
                'r': parser.parse_r_interval('405'),
                'c': parser.parse_c_interval('303')
            },
            {
                'g': parser.parse_g_interval('49925724'),
                'r': parser.parse_r_interval('406'),
                'c': parser.parse_c_interval('*1')
            },
            {
                'g': parser.parse_g_interval('49925671'),
                'r': parser.parse_r_interval('459'),
                'c': parser.parse_c_interval('*54')
            },
        ]
        self.run_cases(tm, test_cases)

    def run_cases(self, tm, test_cases):
        for test_case in test_cases:
            self.assertEquals(tm.g_to_n(test_case['g']), test_case['r'])
            self.assertEquals(tm.n_to_g(test_case['r']), test_case['g'])
            self.assertEquals(tm.n_to_c(test_case['r']), test_case['c'])
            self.assertEquals(tm.c_to_n(test_case['c']), test_case['r'])
            self.assertEquals(tm.g_to_c(test_case['g']), test_case['c'])
            self.assertEquals(tm.c_to_g(test_case['c']), test_case['g'])


if __name__ == '__main__':
    unittest.main()

    # TODO: Reintegrate older tests, especially those with indels
    ### harder tests ###
    #def test_transcriptmapper_TranscriptMapper_1_ZCCHC3(self):
    #    """
    #    reece=> select * from uta.tx_info where ac='NM_033089.6';
    #       gene  | strand |     ac      | cds_start_i | cds_end_i |                 descr                 | summary
    #    --------+--------+-------------+-------------+-----------+---------------------------------------+---------
    #      ZCCHC3 |      1 | NM_033089.6 |          24 |      1236 | zinc finger, CCHC domain containing 3 |
    #
    #    reece=> select * from uta.tx_exons where ac='NM_033089.6';
    #          ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i |    cigar    |
    #    -------------+-----+------+-----------+---------+------------+-----------+---------+-------------+------------------------
    #      NM_033089.6 |   1 | 1    |         0 |    2759 | GRCh37.p10 |    278203 |  280965 | 484M3D2275M | GGAGGATGCTGGGAAGGAGGTAA
    #    """
    #    # http://tinyurl.com/mattx8u
    #    #
    #    # Around the deletion
    #    # http://tinyurl.com/jwt3txg
    #    # 687       690
    #    # C | C G G | C
    #    #   \___ ___/
    #    #     C | C
    #    #      484
    #
    #    ### Add one to g., r., and c. because we are returning hgvs coordinates ###
    #    ac = 'NM_033089.6'
    #    tm = TranscriptMapper(self.hdp, ac, self.ref)
    #    cds = 24 + 1 # hgvs
    #    # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
    #    test_cases = [
    #        {'gs': 278204, 'ge': 278204, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
    #        {'gs': 278214, 'ge': 278214, 'rs': 11, 're': 11, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 11-cds, 'ce': 11-cds},
    #        {'gs': 278204, 'ge': 278214, 'rs': 1, 're': 11, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 11-cds},
    #
    #        # around cds (cds can't be zero)
    #        {'gs': 278227, 'ge': 278227, 'rs': 24, 're': 24, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 24-cds, 'ce': 24-cds},
    #
    #        # beyond cds add 1 due to hgvs
    #        {'gs': 278228, 'ge': 278228, 'rs': 25, 're': 25, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 25-cds+1, 'ce': 25-cds+1},
    #        {'gs': 278229, 'ge': 278229, 'rs': 26, 're': 26, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 26-cds+1, 'ce': 26-cds+1},
    #        {'gs': 280966, 'ge': 280966, 'rs': 2760, 're': 2760, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 2760-cds+1, 'ce': 2760-cds+1},
    #        {'gs': 278687, 'ge': 278687, 'rs': 484, 're': 484, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 484-cds+1, 'ce': 484-cds+1},
    #        {'gs': 278687, 'ge': 278688, 'rs': 484, 're': 485, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 484-cds+1, 'ce': 485-cds+1},
    #        {'gs': 278688, 'ge':278691, 'rs': 485, 're': 485, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 485-cds+1, 'ce': 485-cds+1},
    #
    #        # around cds_start (24) and cds_end (1236), mindful of *coding* del (3D)
    #        {'gs': 278204+24, 'ge': 278204+1236, 'rs': 25, 're': 1237-3, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 25-cds+1, 'ce': 1237-cds-3+1},
    #        {'gs': 280956, 'ge': 280966, 'rs': 2750, 're': 2760, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 2750-cds+1, 'ce': 2760-cds+1},
    #    ]
    #    self.run_cases(tm, test_cases)
    #
    #def test_transcriptmapper_TranscriptMapper_2_MCL1(self):
    #    """
    #    reece=> select * from uta.tx_info where ac='NM_182763.2';
    #      gene | strand |     ac      | cds_start_i | cds_end_i |                      descr                      |
    #    ------+--------+-------------+-------------+-----------+-------------------------------------------------+----------------
    #      MCL1 |     -1 | NM_182763.2 |         208 |      1024 | myeloid cell leukemia sequence 1 (BCL2-related) | This gene encod
    #
    #    reece=> select * from uta.tx_exons where ac='NM_182763.2';
    #          ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i |  g_end_i  |    cigar     |
    #    -------------+-----+------+-----------+---------+------------+-----------+-----------+--------------+---------------------
    #      NM_182763.2 |   1 | 1b   |         0 |     896 | GRCh37.p10 | 150551318 | 150552214 | 896M         |
    #      NM_182763.2 |   2 | 3    |       896 |    3841 | GRCh37.p10 | 150547026 | 150549967 | 1077M4I1864M | GATGGGTTTGTGGAGTTCTT
    #    """
    #
    #    ### Add one to g., r., and c. because we are returning hgvs coordinates ###
    #
    #    ac = 'NM_182763.2'
    #    tm = TranscriptMapper(self.hdp, ac, self.ref)
    #    cds = 208 + 1 # hgvs
    #    test_cases = [
    #        {'gs': 150552215, 'ge': 150552215, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 1-cds, 'ce': 1-cds},
    #        {'gs': 150552214, 'ge': 150552214, 'rs': 2, 're': 2, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 2-cds, 'ce': 2-cds},
    #
    #        # beyond cds add 1 due to hgvs
    #        {'gs': 150552007, 'ge': 150552007, 'rs': 209, 're': 209, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 209-cds+1, 'ce': 209-cds+1},
    #        {'gs': 150547027, 'ge': 150547027, 'rs': 3842, 're': 3842, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 3842-cds+1, 'ce': 3842-cds+1},
    #
    #        #{'gs': 150549968, 'ge': 150549968, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #        {'gs': 150551318, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 1, 'eo': 1, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #        {'gs': 150551318, 'ge': 150551319, 'rs': 897, 're': 897, 'so': 1, 'eo': 0, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #        {'gs': 150551317, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 2, 'eo': 1, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #        {'gs': 150549968, 'ge': 150549969, 'rs': 897, 're': 897, 'so': 0, 'eo': -1, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #        {'gs': 150549969, 'ge': 150549970, 'rs': 897, 're': 897, 'so': -1, 'eo': -2, 'd': hgvs.location.SEQ_START , 'cs': 897-cds+1, 'ce': 897-cds+1},
    #
    #        # exon 2, 4nt insertion ~ r.2760
    #        # See http://tinyurl.com/mwegybw
    #        # The coords of this indel via NW alignment differ from those at NCBI, but are the same canonicalized
    #        # variant.  Nothing to do about that short of running Splign ourselves.  Test a few examples.
    #        {'gs': 150548892, 'ge': 150548892, 'rs': 1973, 're': 1973, 'so': 0, 'eo':0, 'd': hgvs.location.SEQ_START , 'cs': 1973-cds+1, 'ce': 1973-cds+1},
    #        #? {'gs': 150548891, 'ge': 150548892, 'rs': 1972, 're': 1973, 'so': 0, 'eo':0, 'd': hgvs.location.SEQ_START , 'cs': 1972-cds+1, 'ce': 1973-cds+1},
    #        {'gs': 150548890, 'ge': 150548892, 'rs': 1973, 're': 1979, 'so': 0, 'eo':0, 'd': hgvs.location.SEQ_START , 'cs': 1973-cds+1, 'ce': 1979-cds+1},
    #    ]
    #    self.run_cases(tm, test_cases)
    #
    #    ## exon 2, 4nt insertion ~ r.2760
    #    ## See http://tinyurl.com/mwegybw
    #    ## The coords of this indel via NW alignment differ from those at
    #    ## NCBI, but are the same canonicalized variant.  Nothing to do
    #    ## about that short of running Splign ourselves.
    #    #self.assertEqual(tm.n_to_g(1972, 1972), (150548891, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1973), (150548890, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1974), (150548890, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1975), (150548890, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1976), (150548890, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1977), (150548890, 150548891))
    #    #self.assertEqual(tm.n_to_g(1972, 1978), (150548889, 150548891))
    #    #
    #    #self.assertEqual(tm.g_to_n(150548891, 150548891), (1972, 1972, 0, 0))
    #    #self.assertEqual(tm.g_to_n(150548890, 150548891), (1972, 1973, 0, 0))
    #    #self.assertEqual(tm.g_to_n(150548889, 150548891), (1972, 1978, 0, 0))
    #    #
    #    ## around cds_start (208) and cds_end (1024), mindful of *non-coding* ins (4I)
    #    ## i.e., we *don't* need to account for the 4nt insertion here
    #    #self.assertEquals(tm.n_to_c(208, 1024), (0, 1024 - 208, 0, 0))
    #    #self.assertEquals(tm.c_to_n(0, 1024 - 208), (208, 1024, 0, 0))
    #    #self.assertEquals(tm.g_to_c(150552214 - 208, 150552214 - 208), (0, 0, 0, 0))
    #    #self.assertEquals(tm.c_to_g(0, 0), (150552214 - 208, 150552214 - 208))
    #    ## cds_end is in 2nd exon
    #    #self.assertEquals(tm.g_to_c(150549967 - (1024 - 896), 150549967 - (1024 - 896)), (1024 - 208, 1024 - 208, 0, 0))
    #    #self.assertEquals(tm.c_to_g(1024 - 208, 1024 - 208), (150549967 - (1024 - 896), 150549967 - (1024 - 896)))
    #
    #
    #def test_transcriptmapper_TranscriptMapper_3_IFI27L1(self):
    #    """
    #    #reece=> select * from uta.tx_info where ac='NM_145249.2';
    #    #  gene   | chr | strand |     ac      | cds_start_i | cds_end_i |                     descr                     | summary
    #    #---------+-----+--------+-------------+-------------+-----------+-----------------------------------------------+---------
    #    # IFI27L1 | 14  |      1 | NM_145249.2 |         254 |       569 | interferon, alpha-inducible protein 27-like 1 |
    #    #(1 row)
    #    # reece=>select * from uta.tx_exons where ac = 'NM_145249.2';
    #    #
    #    #      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i  | g_cigar | g_seq_a | t_seq_a
    #    # -------------+-----+------+-----------+---------+------------+-----------+----------+---------+---------+---------
    #    #  NM_145249.2 |   1 | 1    |         0 |     157 | GRCh37.p10 |  94547638 | 94547795 | 157M    |         |
    #    #  NM_145249.2 |   2 | 2a   |       157 |     282 | GRCh37.p10 |  94563186 | 94563311 | 125M    |         |
    #    #  NM_145249.2 |   3 | 3    |       282 |     315 | GRCh37.p10 |  94567084 | 94567117 | 33M     |         |
    #    #  NM_145249.2 |   4 | 4    |       315 |     477 | GRCh37.p10 |  94568159 | 94568321 | 162M    |         |
    #    #  NM_145249.2 |   5 | 5    |       477 |     715 | GRCh37.p10 |  94568822 | 94569060 | 238M    |         |
    #    """
    #
    #    ### Add one to g., r., and c. because we are returning hgvs coordinates ###
    #
    #    ac = 'NM_145249.2'
    #    tm = TranscriptMapper(self.hdp, ac, self.ref)
    #    cds = 254 + 1 # hgvs
    #    test_cases = [
    #        #{'gs': 94547639, 'ge': 94547639, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
    #        #{'gs': 94547796, 'ge': 94547796, 'rs': 158, 're': 158, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 158-cds, 'ce': 158-cds},
    #        #{'gs': 94563185, 'ge': 94563185, 'rs': 159, 're': 159, 'so': -2, 'eo': -2, 'd': hgvs.location.SEQ_START, 'cs': 159-cds, 'ce': 159-cds},
    #
    #        # beyond cds add 1 due to hgvs
    #        #{'gs': 94567118, 'ge': 94567120, 'rs': 316, 're': 316, 'so': 0, 'eo': 2, 'd': hgvs.location.SEQ_START, 'cs': 316-cds+1, 'ce': 316-cds+1},
    #        {'gs': 94567115, 'ge': 94567118, 'rs': 313, 're': 316, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 313-cds+1, 'ce': 316-cds+1},
    #
    #        # intron in the middle between exon 1 and 2
    #        #{'gs': 94555500, 'ge': 94555501, 'rs': 157, 're': 158, 'so': 7686, 'eo': -7685, 'd': hgvs.location.SEQ_START, 'cs': 157-cds+1, 'ce': 158-cds+1},
    #        #{'gs': 94555481, 'ge': 94555501, 'rs': 157, 're': 158, 'so': 7686, 'eo': -7685, 'd': hgvs.location.SEQ_START, 'cs': 157-cds+1, 'ce': 158-cds+1},
    #    ]
    #    self.run_cases(tm, test_cases)

    ### ANOTHER POSSIBLE TEST CASE ###
    # reece=> select * from uta.tx_info where ac = 'NM_145171.3';
    #  gene  | strand |     ac      | cds_start_i | cds_end_i |            descr            | summary
    # -------+--------+-------------+-------------+-----------+-----------------------------+-----------------------------------
    #  GPHB5 |     -1 | NM_145171.3 |          57 |       450 | glycoprotein hormone beta 5 | GPHB5 is a cystine knot-forming...
    #
    # reece=> select * from uta.tx_exons where ac = 'NM_145171.3' order by g_start_i;
    #      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i  |   cigar   | g_seq_a
    # -------------+-----+------+-----------+---------+------------+-----------+----------+-----------+-------------------------
    #  NM_145171.3 |   3 | 3    |       261 |     543 | GRCh37.p10 |  63779548 | 63779830 | 282M      |
    #  NM_145171.3 |   2 | 2    |        56 |     261 | GRCh37.p10 |  63784360 | 63784564 | 156M1I48M | CATGAAGCTGGCATTCCTCTT...
    #  NM_145171.3 |   1 | 1    |         0 |      56 | GRCh37.p10 |  63785537 | 63785593 | 56M       |
    # def test_transcriptmapper_TranscriptMapper_GPHB5(self):
    #     ac = 'NM_145171.3'
    #     tm = TranscriptMapper(self.hdp,ac,self.ref)
    #     pass

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
