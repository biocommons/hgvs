import unittest
import hgvs.location
from uta.db.transcriptdb import TranscriptDB
from hgvs.exceptions import *
from hgvs.transcriptmapper import TranscriptMapper

class Test_transcriptmapper(unittest.TestCase):
    ref = 'GRCh37.p10'

    def setUp(self):
        self.db = TranscriptDB()

    def test_transcriptmapper_failures(self):
        self.assertRaises(HGVSError, TranscriptMapper, self.db, ref=self.ref, ac='bogus')
        self.assertRaises(HGVSError, TranscriptMapper, self.db, ref='bogus', ac='NM_033089.6')

    def test_transcriptmapper_TranscriptMapper_LCE3C(self):
        """NM_178434.2: LCE3C single exon, strand = +1, all coordinate input/output are in HGVS"""
        ac = 'NM_178434.2'
        tm = TranscriptMapper(self.db, ac, self.ref)
        cds = 70 + 1 # hgvs
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'gs': 152573138, 'ge': 152573138, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
            {'gs': 152573138+2, 'ge': 152573138+2, 'rs': 3, 're': 3, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 3-cds, 'ce': 3-cds},
            # cds
            {'gs': 152573138+70, 'ge': 152573138+70, 'rs': 71, 're': 71, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 71-cds+1, 'ce': 71-cds+1},
            # beyond cds add 1 due to hgvs
            {'gs': 152573562, 'ge': 152573562, 'rs': 425, 're': 425, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 425-cds+1, 'ce': 425-cds+1},
            {'gs': 152573562-2, 'ge': 152573562-2, 'rs': 423, 're': 423, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 423-cds+1, 'ce': 423-cds+1},
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_HIST3H2A(self):
        """NM_033445.2: LCE3C single exon, strand = -1, all coordinate input/output are in HGVS"""
        ac = 'NM_033445.2'
        tm = TranscriptMapper(self.db, ac, self.ref)
        cds = 42 + 1 # hgvs
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'gs': 228645560, 'ge': 228645560, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
            {'gs': 228645560-2, 'ge': 228645560-2, 'rs': 3, 're': 3, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 3-cds, 'ce': 3-cds},
            # beyond cds add 1 due to hgvs
            {'gs': 228645065, 'ge': 228645065, 'rs': 496, 're': 496, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 496-cds+1, 'ce': 496-cds+1},
            {'gs': 228645065+2, 'ge': 228645065+2, 'rs': 494, 're': 494, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 494-cds+1, 'ce': 494-cds+1},
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_LCE2B(self):
        """NM_014357.4: LCE2B, two exons, strand = +1, all coordinate input/output are in HGVS"""
        ac = 'NM_014357.4'
        tm = TranscriptMapper(self.db, ac, self.ref)
        cds = 54 + 1 # hgvs
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'gs': 152658599, 'ge': 152658599, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
            {'gs': 152658599+2, 'ge': 152658599+2, 'rs': 3, 're': 3, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 3-cds, 'ce': 3-cds},
            # around end of exon 1
            {'gs': 152658632, 'ge': 152658632, 'rs': 34, 're': 34, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 34-cds, 'ce': 34-cds},
            {'gs': 152658632, 'ge': 152658632+2, 'rs': 34, 're': 34, 'so': 0, 'eo': 2, 'd': hgvs.location.SEQ_START, 'cs': 34-cds, 'ce': 34-cds},
            {'gs': 152658632+1, 'ge': 152658632+1, 'rs': 34, 're': 34, 'so': 1, 'eo': 1, 'd': hgvs.location.SEQ_START, 'cs': 34-cds, 'ce': 34-cds},
            # around beginning of exon 2
            {'gs': 152659300, 'ge': 152659300, 'rs': 35, 're': 35, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 35-cds, 'ce': 35-cds},
            {'gs': 152659300-2, 'ge': 152659300-2, 'rs': 35, 're': 35, 'so': -2, 'eo': -2, 'd': hgvs.location.SEQ_START, 'cs': 35-cds, 'ce': 35-cds},
            {'gs': 152659300-2, 'ge': 152659301, 'rs': 35, 're': 36, 'so': -2, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 35-cds, 'ce': 36-cds},
            # beyond cds add 1 due to hgvs
            {'gs': 152659877, 'ge': 152659877, 'rs': 612, 're': 612, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 612-cds+1, 'ce': 612-cds+1},
            {'gs': 152659877-2, 'ge': 152659877-2, 'rs': 610, 're': 610, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 610-cds+1, 'ce': 610-cds+1},
        ]
        self.run_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_PTH2(self):
        """NM_178449.3: PTH2, two exons, strand = -1, all coordinate input/output are in HGVS"""
        ac = 'NM_178449.3'
        tm = TranscriptMapper(self.db, ac, self.ref)
        cds = 102 + 1 # hgvs
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'gs': 49926698, 'ge': 49926698, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 1-cds, 'ce': 1-cds},
            {'gs': 49926698-2, 'ge': 49926698-2, 'rs': 3, 're': 3, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 3-cds, 'ce': 3-cds},
            # around end of exon 1
            {'gs': 49926469, 'ge': 49926469, 'rs': 230, 're': 230, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 230-cds+1, 'ce': 230-cds+1},
            {'gs': 49926469-2, 'ge': 49926469-2, 'rs': 230, 're': 230, 'so': 2, 'eo': 2, 'd': hgvs.location.SEQ_START, 'cs': 230-cds+1, 'ce': 230-cds+1},
            {'gs': 49926470, 'ge': 49926470, 'rs': 229, 're': 229, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 229-cds+1, 'ce': 229-cds+1},
            {'gs': 49926469-2, 'ge': 49926470, 'rs': 229, 're': 230, 'so': 0, 'eo': 2, 'd': hgvs.location.SEQ_START, 'cs': 229-cds+1, 'ce': 230-cds+1},
            # around beginning of exon 2
            {'gs': 49925899, 'ge': 49925899, 'rs': 231, 're': 231, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 231-cds+1, 'ce': 231-cds+1},
            {'gs': 49925899+1, 'ge': 49925899+1, 'rs': 231, 're': 231, 'so': -1, 'eo': -1, 'd': hgvs.location.SEQ_START, 'cs': 231-cds+1, 'ce': 231-cds+1},
            {'gs': 49925899-1, 'ge': 49925899+2, 'rs': 231, 're': 232, 'so': -2, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 231-cds+1, 'ce': 232-cds+1},
            # beyond cds add 1 due to hgvs
            {'gs': 49925671, 'ge': 49925671, 'rs': 459, 're': 459, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 459-cds+1, 'ce': 459-cds+1},
            {'gs': 49925671+2, 'ge': 49925671+2, 'rs': 457, 're': 457, 'so': 0, 'eo': 0, 'd': hgvs.location.SEQ_START, 'cs': 457-cds+1, 'ce': 457-cds+1},
        ]
        self.run_cases(tm, test_cases)


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
    #    tm = TranscriptMapper(self.db, ac, self.ref)
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
    #    tm = TranscriptMapper(self.db, ac, self.ref)
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
    #    #self.assertEqual(tm.r_to_g(1972, 1972), (150548891, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1973), (150548890, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1974), (150548890, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1975), (150548890, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1976), (150548890, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1977), (150548890, 150548891))
    #    #self.assertEqual(tm.r_to_g(1972, 1978), (150548889, 150548891))
    #    #
    #    #self.assertEqual(tm.g_to_r(150548891, 150548891), (1972, 1972, 0, 0))
    #    #self.assertEqual(tm.g_to_r(150548890, 150548891), (1972, 1973, 0, 0))
    #    #self.assertEqual(tm.g_to_r(150548889, 150548891), (1972, 1978, 0, 0))
    #    #
    #    ## around cds_start (208) and cds_end (1024), mindful of *non-coding* ins (4I)
    #    ## i.e., we *don't* need to account for the 4nt insertion here
    #    #self.assertEquals(tm.r_to_c(208, 1024), (0, 1024 - 208, 0, 0))
    #    #self.assertEquals(tm.c_to_r(0, 1024 - 208), (208, 1024, 0, 0))
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
    #    tm = TranscriptMapper(self.db, ac, self.ref)
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

    def run_cases(self, tm, test_cases):
        for test_case in test_cases:
            g = hgvs.location.Interval(start=test_case['gs'], end=test_case['ge'])
            r = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['rs'], offset=test_case['so'], datum=test_case['d']),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['re'], offset=test_case['eo'], datum=test_case['d']))
            c = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['cs'], offset=test_case['so'], datum=test_case['d'] + 1),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['ce'], offset=test_case['eo'], datum=test_case['d'] + 1))
            self.assertEquals(tm.hgvsg_to_hgvsr(g), r)
            self.assertEquals(tm.hgvsr_to_hgvsg(r), g)
            self.assertEquals(tm.hgvsr_to_hgvsc(r), c)
            self.assertEquals(tm.hgvsc_to_hgvsr(c), r)
            self.assertEquals(tm.hgvsg_to_hgvsc(g), c)
            self.assertEquals(tm.hgvsc_to_hgvsg(c), g)


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
#     tm = TranscriptMapper(self.db,ac,self.ref)
#     pass

if __name__ == '__main__':
    unittest.main()

