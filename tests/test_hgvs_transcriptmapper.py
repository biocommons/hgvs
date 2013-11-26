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

    def test_transcriptmapper_TranscriptMapper_1_ZCCHC3(self):
        """
        reece=> select * from uta.tx_info where ac='NM_033089.6';
           gene  | strand |     ac      | cds_start_i | cds_end_i |                 descr                 | summary
        --------+--------+-------------+-------------+-----------+---------------------------------------+---------
          ZCCHC3 |      1 | NM_033089.6 |          24 |      1236 | zinc finger, CCHC domain containing 3 |

        reece=> select * from uta.tx_exons where ac='NM_033089.6';
              ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i |    cigar    |
        -------------+-----+------+-----------+---------+------------+-----------+---------+-------------+------------------------
          NM_033089.6 |   1 | 1    |         0 |    2759 | GRCh37.p10 |    278203 |  280965 | 484M3D2275M | GGAGGATGCTGGGAAGGAGGTAA
        """
        # http://tinyurl.com/mattx8u
        #
        # Around the deletion
        # http://tinyurl.com/jwt3txg
        # 687       690
        # C | C G G | C
        #   \___ ___/
        #     C | C
        #      484

        ### Add one to r. and c. because we are returning hgvs coordinates ###

        ac = 'NM_033089.6'
        tm = TranscriptMapper(self.db, ac, self.ref)
        # test cases: type indicates the type of conversation the test will work.  Some are only for g -> r or r -> g
        # due to intronic regions
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'type': 'grc', 'gs': 278203, 'ge': 278203, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': 0, 'cs': 1-24, 'ce': 1-24},
            {'type': 'grc', 'gs': 278213, 'ge': 278213, 'rs': 11, 're': 11, 'so': 0, 'eo': 0, 'd': 0, 'cs': 11-24, 'ce': 11-24},
            {'type': 'grc', 'gs': 280965, 'ge': 280965, 'rs': 2760, 're': 2760, 'so': 0, 'eo': 0, 'd': 0, 'cs': 2760-24, 'ce': 2760-24},
            {'type': 'grc', 'gs': 278203, 'ge': 278213, 'rs': 1, 're': 11, 'so': 0, 'eo': 0, 'd': 0, 'cs': 1-24, 'ce': 11-24},
            {'type': 'grc', 'gs': 278227, 'ge': 278227, 'rs': 25, 're': 25, 'so': 0, 'eo': 0, 'd': 0, 'cs': 25-24, 'ce': 25-24},
            {'type': 'grc', 'gs': 278686, 'ge': 278686, 'rs': 484, 're': 484, 'so': 0, 'eo': 0, 'd': 0, 'cs': 484-24, 'ce': 484-24},
            {'type': 'grc', 'gs': 278686, 'ge': 278687, 'rs': 484, 're': 485, 'so': 0, 'eo': 0, 'd': 0, 'cs': 484-24, 'ce': 485-24},
            {'type': 'grc', 'gs': 278687, 'ge': 278690, 'rs': 485, 're': 485, 'so': 0, 'eo': 0, 'd': 0, 'cs': 485-24, 'ce': 485-24},

            # around cds_start (24) and cds_end (1236), mindful of *coding* del (3D)
            {'type': 'grc', 'gs': 278203+24, 'ge': 278203+1236, 'rs': 25, 're': 1237-3, 'so': 0, 'eo': 0, 'd': 0, 'cs': 25-24, 'ce': 1237-24-3},

            {'type': 'grc', 'gs': 280955, 'ge': 280965, 'rs': 2750, 're': 2760, 'so': 0, 'eo': 0, 'd': 0, 'cs': 2750-24, 'ce': 2760-24},
        ]
        self.run_test_cases(tm, test_cases)

    def test_transcriptmapper_TranscriptMapper_2_MCL1(self):
        """
        reece=> select * from uta.tx_info where ac='NM_182763.2';
          gene | strand |     ac      | cds_start_i | cds_end_i |                      descr                      |
        ------+--------+-------------+-------------+-----------+-------------------------------------------------+----------------
          MCL1 |     -1 | NM_182763.2 |         208 |      1024 | myeloid cell leukemia sequence 1 (BCL2-related) | This gene encod

        reece=> select * from uta.tx_exons where ac='NM_182763.2';
              ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i |  g_end_i  |    cigar     |
        -------------+-----+------+-----------+---------+------------+-----------+-----------+--------------+---------------------
          NM_182763.2 |   1 | 1b   |         0 |     896 | GRCh37.p10 | 150551318 | 150552214 | 896M         |
          NM_182763.2 |   2 | 3    |       896 |    3841 | GRCh37.p10 | 150547026 | 150549967 | 1077M4I1864M | GATGGGTTTGTGGAGTTCTT
        """

        ### Add one to r. and c. because we are returning hgvs coordinates ###

        ac = 'NM_182763.2'
        tm = TranscriptMapper(self.db, ac, self.ref)
        # test cases: type indicates the type of conversation the test will work.  Some are only for g -> r or r -> g
        # due to intronic regions
        test_cases = [
            {'type': 'grc', 'gs': 150552214, 'ge': 150552214, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': 0, 'cs': 1-208, 'ce': 1-208},
            {'type': 'grc', 'gs': 150552213, 'ge': 150552213, 'rs': 2, 're': 2, 'so': 0, 'eo': 0, 'd': 0, 'cs': 2-208, 'ce': 2-208},
            {'type': 'g', 'gs': 150551318, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'g', 'gs': 150549967, 'ge': 150549967, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549967, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551317, 'ge': 150551317, 'rs': 897, 're': 897, 'so': 1, 'eo': 1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551317, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 1, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551316, 'ge': 150551317, 'rs': 897, 're': 897, 'so': 2, 'eo': 1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549967, 'ge': 150549968, 'rs': 897, 're': 897, 'so': 0, 'eo': -1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549968, 'ge': 150549969, 'rs': 897, 're': 897, 'so': -1, 'eo': -2, 'd': 0, 'cs': 897-208, 'ce': 897-208},

            # exon 2, 4nt insertion ~ r.2760
            # See http://tinyurl.com/mwegybw
            # The coords of this indel via NW alignment differ from those at NCBI, but are the same canonicalized
            # variant.  Nothing to do about that short of running Splign ourselves.  Test a few examples.
            {'type': 'grc', 'gs': 150548891, 'ge': 150548891, 'rs': 1973, 're': 1973, 'so': 0, 'eo':0, 'd': 0, 'cs': 1973-208, 'ce': 1973-208},
            #? {'type': 'grc', 'gs': 150548890, 'ge': 150548891, 'rs': 1972, 're': 1973, 'so': 0, 'eo':0, 'd': 0, 'cs': 1972-208, 'ce': 1973-208},
            {'type': 'grc', 'gs': 150548889, 'ge': 150548891, 'rs': 1973, 're': 1979, 'so': 0, 'eo':0, 'd': 0, 'cs': 1973-208, 'ce': 1979-208},
        ]
        self.run_test_cases(tm, test_cases)

        ## exon 2, 4nt insertion ~ r.2760
        ## See http://tinyurl.com/mwegybw
        ## The coords of this indel via NW alignment differ from those at
        ## NCBI, but are the same canonicalized variant.  Nothing to do
        ## about that short of running Splign ourselves.
        #self.assertEqual(tm.r_to_g(1972, 1972), (150548891, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1973), (150548890, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1974), (150548890, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1975), (150548890, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1976), (150548890, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1977), (150548890, 150548891))
        #self.assertEqual(tm.r_to_g(1972, 1978), (150548889, 150548891))
        #
        #self.assertEqual(tm.g_to_r(150548891, 150548891), (1972, 1972, 0, 0))
        #self.assertEqual(tm.g_to_r(150548890, 150548891), (1972, 1973, 0, 0))
        #self.assertEqual(tm.g_to_r(150548889, 150548891), (1972, 1978, 0, 0))
        #
        ## around cds_start (208) and cds_end (1024), mindful of *non-coding* ins (4I)
        ## i.e., we *don't* need to account for the 4nt insertion here
        #self.assertEquals(tm.r_to_c(208, 1024), (0, 1024 - 208, 0, 0))
        #self.assertEquals(tm.c_to_r(0, 1024 - 208), (208, 1024, 0, 0))
        #self.assertEquals(tm.g_to_c(150552214 - 208, 150552214 - 208), (0, 0, 0, 0))
        #self.assertEquals(tm.c_to_g(0, 0), (150552214 - 208, 150552214 - 208))
        ## cds_end is in 2nd exon
        #self.assertEquals(tm.g_to_c(150549967 - (1024 - 896), 150549967 - (1024 - 896)), (1024 - 208, 1024 - 208, 0, 0))
        #self.assertEquals(tm.c_to_g(1024 - 208, 1024 - 208), (150549967 - (1024 - 896), 150549967 - (1024 - 896)))


    def test_transcriptmapper_TranscriptMapper_3_IFI27L1(self):
        """
        #reece=> select * from uta.tx_info where ac='NM_145249.2';
        #  gene   | chr | strand |     ac      | cds_start_i | cds_end_i |                     descr                     | summary
        #---------+-----+--------+-------------+-------------+-----------+-----------------------------------------------+---------
        # IFI27L1 | 14  |      1 | NM_145249.2 |         254 |       569 | interferon, alpha-inducible protein 27-like 1 |
        #(1 row)
        # reece=>select * from uta.tx_exons where ac = 'NM_145249.2';
        #
        #      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i  | g_cigar | g_seq_a | t_seq_a
        # -------------+-----+------+-----------+---------+------------+-----------+----------+---------+---------+---------
        #  NM_145249.2 |   1 | 1    |         0 |     157 | GRCh37.p10 |  94547638 | 94547795 | 157M    |         |
        #  NM_145249.2 |   2 | 2a   |       157 |     282 | GRCh37.p10 |  94563186 | 94563311 | 125M    |         |
        #  NM_145249.2 |   3 | 3    |       282 |     315 | GRCh37.p10 |  94567084 | 94567117 | 33M     |         |
        #  NM_145249.2 |   4 | 4    |       315 |     477 | GRCh37.p10 |  94568159 | 94568321 | 162M    |         |
        #  NM_145249.2 |   5 | 5    |       477 |     715 | GRCh37.p10 |  94568822 | 94569060 | 238M    |         |
        """

        ### Add one to r. and c. because we are returning hgvs coordinates ###

        ac = 'NM_145249.2'
        tm = TranscriptMapper(self.db, ac, self.ref)

        # test cases: type indicates the type of conversation the test will work.  Some are only for g -> r or r -> g
        # due to intronic regions
        test_cases = [
            {'type': 'grc', 'gs': 94547638, 'ge': 94547638, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': 0, 'cs': 1-254, 'ce': 1-254},
            {'type': 'g', 'gs': 94547795, 'ge': 94547795, 'rs': 158, 're': 158, 'so': 0, 'eo': 0, 'd': 0, 'cs': 158-254, 'ce': 158-254},
            {'type': 'grc', 'gs': 94563184, 'ge': 94563184, 'rs': 158, 're': 158, 'so': -2, 'eo': -2, 'd': 0, 'cs': 158-254, 'ce': 158-254},
            {'type': 'grc', 'gs': 94567117, 'ge': 94567119, 'rs': 316, 're': 316, 'so': 0, 'eo': 2, 'd': 0, 'cs': 316-254, 'ce': 316-254},

            # intron in the middle between exon 1 and 2
            {'type': 'grc', 'gs': 94555480, 'ge': 94555500, 'rs': 158, 're': 158, 'so': 7685, 'eo': -7686, 'd': 0, 'cs': 158-254, 'ce': 158-254},
        ]
        self.run_test_cases(tm, test_cases)

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

    def run_test_cases(self, tm, test_cases):
        for test_case in test_cases:
            g = hgvs.location.Interval(start=test_case['gs'], end=test_case['ge'])
            r = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['rs'], offset=test_case['so'], datum=test_case['d']),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['re'], offset=test_case['eo'], datum=test_case['d']))
            c = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['cs'], offset=test_case['so'], datum=test_case['d'] + 1),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['ce'], offset=test_case['eo'], datum=test_case['d'] + 1))
            if test_case['type'] == 'grc' or test_case['type'] == 'g':
                self.assertEquals(tm.hgvsg_to_hgvsr(g), r)
            if test_case['type'] == 'grc':
                self.assertEquals(tm.hgvsr_to_hgvsg(r), g)
                self.assertEquals(tm.hgvsr_to_hgvsc(r), c)
                self.assertEquals(tm.hgvsc_to_hgvsr(c), r)
                self.assertEquals(tm.hgvsg_to_hgvsc(g), c)
                self.assertEquals(tm.hgvsc_to_hgvsg(c), g)


if __name__ == '__main__':
    unittest.main()
