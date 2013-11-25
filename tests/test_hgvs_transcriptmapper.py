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


    # reece=> select * from uta.tx_info where ac='NM_033089.6';
    #   gene  | strand |     ac      | cds_start_i | cds_end_i |                 descr                 | summary
    # --------+--------+-------------+-------------+-----------+---------------------------------------+---------
    #  ZCCHC3 |      1 | NM_033089.6 |          24 |      1236 | zinc finger, CCHC domain containing 3 |
    #
    # reece=> select * from uta.tx_exons where ac='NM_033089.6';
    #      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i |    cigar    |
    # -------------+-----+------+-----------+---------+------------+-----------+---------+-------------+------------------------
    #  NM_033089.6 |   1 | 1    |         0 |    2759 | GRCh37.p10 |    278203 |  280965 | 484M3D2275M | GGAGGATGCTGGGAAGGAGGTAA

    def test_transcriptmapper_TranscriptMapper_1_ZCCHC3(self):
        ###### Add one because we are returning hgvs coordinates ######
        ac = 'NM_033089.6'
        tm = TranscriptMapper(self.db, ac, self.ref)
        # test cases: type indicates the type of conversation the test will work.  Some are only for g -> r or r -> g
        # due to intronic regions
        # gs, ge = genomic start/end; rs,re = rna start/end; cs, ce = cdna start/end; so, eo = start offset/end offset
        test_cases = [
            {'type': 'grc', 'gs': 278203, 'ge': 278203, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': 0, 'cs': -23, 'ce': -23},
            {'type': 'grc', 'gs': 278213, 'ge': 278213, 'rs': 11, 're': 11, 'so': 0, 'eo': 0, 'd': 0, 'cs': -13, 'ce': -13},
            {'type': 'grc', 'gs': 280965, 'ge': 280965, 'rs': 2760, 're': 2760, 'so': 0, 'eo': 0, 'd': 0, 'cs': 2736, 'ce': 2736},
            {'type': 'grc', 'gs': 278203, 'ge': 278213, 'rs': 1, 're': 11, 'so': 0, 'eo': 0, 'd': 0, 'cs': -23, 'ce': -13},
            {'type': 'grc', 'gs': 278227, 'ge': 278227, 'rs': 25, 're': 25, 'so': 0, 'eo': 0, 'd': 0, 'cs': 1, 'ce': 1},
        ]
        self.run_test_cases(tm, test_cases)

    #
    #    # http://tinyurl.com/mattx8u
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278203), (1, 1, 0 , 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278213), (10, 10, 0 , 0))
    #
    #
    #    ##self.assertEquals(tm.g_to_r(278203, 278203), (0, 0, 0, 0))
    #    ##self.assertEquals(tm.g_to_r(278203, 278213), (0, 10, 0, 0))
    #    #self.assertEquals(tm.r_to_g(0, 0), (278203, 278203))
    #    #self.assertEquals(tm.r_to_g(0, 10), (278203, 278213))
    #    #
    #    ## http://tinyurl.com/kvrtkan
    #    self.assertEquals(tm.hgvsg_to_hgvsr(280955), (2749, 2749, 0, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(280965), (2759, 2759, 0, 0))
    #    ## self.assertEquals( tm.g_to_r(2759,2759),(280965,280965)  )  # 0-width not supported
    #    ##self.assertEquals(tm.g_to_r(280955, 280965), (2749, 2759, 0, 0))
    #    #self.assertEquals(tm.r_to_g(2759, 2759), (280965, 280965))
    #    #self.assertEquals(tm.r_to_g(2749, 2759), (280955, 280965))
    #    #
    #    ## Around the deletion
    #    ## http://tinyurl.com/jwt3txg
    #    ## 687       690
    #    ## C | C G G | C
    #    ##   \___ ___/
    #    ##     C | C
    #    ##      484
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278686), (483, 483, 0, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278687), (484, 484, 0, 0))
    #    ##self.assertEquals(tm.g_to_r(278686, 278687), (483, 484, 0, 0))
    #    #self.assertEquals(tm.r_to_g(483, 484), (278686, 278687))
    #    ## self.assertEquals( tm.r_to_g(484,484), (278687, 278690) ) # 0-width not supported yet
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278690), (484, 484, 0, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(278691), (485, 485, 0, 0))
    #    ##self.assertEquals(tm.g_to_r(278690, 278691), (484, 485, 0, 0))
    #    #self.assertEquals(tm.r_to_g(484, 485), (278690, 278691))
    #    #
    #    ## around cds_start (24) and cds_end (1236), mindful of *coding* del (3D)
    #    #self.assertEquals(tm.r_to_c(24, 1236), (0, 1236 - 24, 0, 0))
    #    #self.assertEquals(tm.c_to_r(0, 1236 - 24), (24, 1236, 0, 0))
    #    #self.assertEquals(tm.g_to_c(278203 + 24, 278203 + 1236), (0, 1236 - 24 - 3, 0, 0))
    #    #self.assertEquals(tm.c_to_g(0, 1236 - 24 - 3), (278203 + 24, 278203 + 1236))


    # reece=> select * from uta.tx_info where ac='NM_182763.2';
    #  gene | strand |     ac      | cds_start_i | cds_end_i |                      descr                      |
    # ------+--------+-------------+-------------+-----------+-------------------------------------------------+----------------
    #  MCL1 |     -1 | NM_182763.2 |         208 |      1024 | myeloid cell leukemia sequence 1 (BCL2-related) | This gene encod
    #
    # reece=> select * from uta.tx_exons where ac='NM_182763.2';
    #      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i |  g_end_i  |    cigar     |
    # -------------+-----+------+-----------+---------+------------+-----------+-----------+--------------+---------------------
    #  NM_182763.2 |   1 | 1b   |         0 |     896 | GRCh37.p10 | 150551318 | 150552214 | 896M         |
    #  NM_182763.2 |   2 | 3    |       896 |    3841 | GRCh37.p10 | 150547026 | 150549967 | 1077M4I1864M | GATGGGTTTGTGGAGTTCTT

    def test_transcriptmapper_TranscriptMapper_2_MCL1(self):
        ###### Add one because we are returning hgvs coordinates ######
        ac = 'NM_182763.2'
        tm = TranscriptMapper(self.db, ac, self.ref)
        # test cases: type indicates the type of conversation the test will work.  Some are only for g -> r or r -> g
        # due to intronic regions
        test_cases = [
            {'type': 'grc', 'gs': 150552214, 'ge': 150552214, 'rs': 1, 're': 1, 'so': 0, 'eo': 0, 'd': 0, 'cs': -207, 'ce': -207},
            {'type': 'grc', 'gs': 150552213, 'ge': 150552213, 'rs': 2, 're': 2, 'so': 0, 'eo': 0, 'd': 0, 'cs': -206, 'ce': -206},
            {'type': 'g', 'gs': 150551318, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'g', 'gs': 150549967, 'ge': 150549967, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549967, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 0, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551317, 'ge': 150551317, 'rs': 897, 're': 897, 'so': 1, 'eo': 1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551317, 'ge': 150551318, 'rs': 897, 're': 897, 'so': 1, 'eo': 0, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150551316, 'ge': 150551317, 'rs': 897, 're': 897, 'so': 2, 'eo': 1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549967, 'ge': 150549968, 'rs': 897, 're': 897, 'so': 0, 'eo': -1, 'd': 0, 'cs': 897-208, 'ce': 897-208},
            {'type': 'grc', 'gs': 150549968, 'ge': 150549969, 'rs': 897, 're': 897, 'so': -1, 'eo': -2, 'd': 0, 'cs': 897-208, 'ce': 897-208},
        ]
        self.run_test_cases(tm, test_cases)

        #self.assertEquals(tm.hgvsg_to_hgvsr(150552214), (1, 1, 0, 0))
        #self.assertEquals(tm.hgvsg_to_hgvsr(150547026), (3842, 3842, 0, 0))
        #self.assertEqual(tm.g_to_r(150552214, 150552214), (0, 0, 0, 0))
        #self.assertEqual(tm.r_to_g(0, 0, 0, 0), ((150552214, 150552214)))
        #
        ## ends of exon 1 and 2 -> map to same r coord.
        #self.assertEqual(tm.g_to_r(150549967, 150549967), (896, 896, 0, 0))
        #self.assertEqual(tm.g_to_r(150551318, 150551318), (896, 896, 0, 0))
        #self.assertEqual(tm.r_to_g(896, 896, 0, 0), (150549967, 150551318))
        #
        #self.assertEqual(tm.g_to_r(150547026, 150547026), (3841, 3841, 0, 0))
        #self.assertEqual(tm.r_to_g(3841, 3841, 0, 0), (150547026, 150547026))
        #
        ## check that the default 0-based offsets are working
        #self.assertEqual(tm.r_to_g(3841, 3841), (150547026, 150547026))
        #
        ## exon 1, around the start
#        self.assertEquals(tm.hgvsg_to_hgvsr(150552213), (896, 896, 0, 0))
#        self.assertEquals(tm.hgvsg_to_hgvsr(150552214), (897, 897, 0, 0))
#        self.assertEquals(tm.hgvsg_to_hgvsr(150552215), (896, 896, -1, 0))
#        self.assertEquals(tm.hgvsg_to_hgvsr(150552216), (896, 896, -2, 0))

        ##self.assertEqual(tm.g_to_r(150552213, 150552214), (0, 1, 0, 0))
        #self.assertEqual(tm.r_to_g(0, 1, 0, 0), (150552213, 150552214))
        #
#        self.assertEquals(tm.hgvsg_to_hgvsr(150551319), (896, 896, -1, 0))
        ##self.assertEqual(tm.g_to_r(150551318, 150551319), (895, 896, 0, 0))
        #self.assertEqual(tm.r_to_g(895, 896, 0, 0), (150551318, 150551319))
        #
        ## intron 1nt before exon 1
#        self.assertEquals(tm.hgvsg_to_hgvsr(150551317), (896, 896, 1, 0))
#        self.assertEquals(tm.hgvsg_to_hgvsr(150551318), (896, 896, 0, 0))
        #self.assertEqual(tm.g_to_r(150551317, 150551318), (896, 896, 0, 1))
        #self.assertEqual(tm.r_to_g(896, 896, 0, 1), (150551317, 150551318))
        #
        ## intron 3nts before exon 1
        #self.assertEqual(tm.g_to_r(150551315, 150551317), (896, 896, 1, 3))
        #self.assertEqual(tm.r_to_g(896, 896, 1, 3), (150551315, 150551317))
        #
        ## intron 3nts after exon 2
        #self.assertEqual(tm.g_to_r(150549968, 150549970), (896, 896, -3, -1))
        #self.assertEqual(tm.r_to_g(896, 896, -3, -1), (150549968, 150549970))
        #
        ## intron 1nt after exon 2
        #self.assertEqual(tm.g_to_r(150549967, 150549968), (896, 896, -1, 0))
        #self.assertEqual(tm.r_to_g(896, 896, -1, 0), (150549967, 150549968))
        #
        ## exon 2, 1nt at either end
        #self.assertEqual(tm.g_to_r(150549966, 150549967), (896, 897, 0, 0))
        #self.assertEqual(tm.r_to_g(896, 897, 0, 0), (150549966, 150549967))
        #
        #self.assertEqual(tm.g_to_r(150547026, 150547027), (3840, 3841, 0, 0))
        #self.assertEqual(tm.r_to_g(3840, 3841, 0, 0), (150547026, 150547027))
        #
        #
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
    #def test_transcriptmapper_TranscriptMapper_3_IFI27L1(self):
    #    ac = 'NM_145249.2'
    #    tm = TranscriptMapper(self.db, ac, self.ref)
    #
    #    # some base cases
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94547795), (157, 157, 0, 0))
    #
    #    #self.assertEqual(tm.r_to_g(0, 157, 0, 0), (94547638, 94547795))
    #    #self.assertEqual(tm.g_to_r(94547638, 94547795), (0, 157, 0, 0))
    #    #
    #    #self.assertEqual(tm.r_to_g(157, 282, 0, 0), (94563186, 94563311))
    #    #self.assertEqual(tm.g_to_r(94563186, 94563311), (157, 282, 0, 0))
    #    #
    #    #self.assertEqual(tm.r_to_g(477, 715), (94568822, 94569060))
    #    #self.assertEqual(tm.g_to_r(94568822, 94569060), (477, 715, 0, 0))
    #    #
    #    # intron around end of exon 1
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94547795 + 1), (157, 157, 1, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94547795 + 2), (157, 157, 2, 0))
    #    #self.assertEqual(tm.r_to_g(157, 157, 0, 1), (94547795, 94547796))
    #    #self.assertEqual(tm.g_to_r(94547795, 94547796), (157, 157, 0, 1))
    #    #
    #    ## around beginning of exon 2
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94563186 - 2), (157, 157, -2, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94563186 - 1), (157, 157, -1, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94563186 - 0), (157, 157, 0, 0))
    #    #self.assertEqual(tm.r_to_g(157, 157, -1, 0), (94563185, 94563186))
    #    #self.assertEqual(tm.g_to_r(94563185, 94563186), (157, 157, -1, 0))
    #    #
    #    ## intron in the middle between exon 1 and 2
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94555480), (157, 157, 7685, 0))
    #    self.assertEquals(tm.hgvsg_to_hgvsr(94555500), (157, 157, -7686, 0))
    #    #self.assertEqual(tm.g_to_r(94555480, 94555500), (157, 157, 7685, -7686))
    #    #self.assertEqual(tm.r_to_g(157, 157, 7685, -7686), (94555480, 94555500))
    #
    #
    #
    #
    #
    #
    #
    #
    ### reece=> select * from uta.tx_info where ac = 'NM_145171.3';
    ###  gene  | strand |     ac      | cds_start_i | cds_end_i |            descr            | summary
    ### -------+--------+-------------+-------------+-----------+-----------------------------+-----------------------------------
    ###  GPHB5 |     -1 | NM_145171.3 |          57 |       450 | glycoprotein hormone beta 5 | GPHB5 is a cystine knot-forming...
    ###
    ### reece=> select * from uta.tx_exons where ac = 'NM_145171.3' order by g_start_i;
    ###      ac      | ord | name | t_start_i | t_end_i |    ref     | g_start_i | g_end_i  |   cigar   | g_seq_a
    ### -------------+-----+------+-----------+---------+------------+-----------+----------+-----------+-------------------------
    ###  NM_145171.3 |   3 | 3    |       261 |     543 | GRCh37.p10 |  63779548 | 63779830 | 282M      |
    ###  NM_145171.3 |   2 | 2    |        56 |     261 | GRCh37.p10 |  63784360 | 63784564 | 156M1I48M | CATGAAGCTGGCATTCCTCTT...
    ###  NM_145171.3 |   1 | 1    |         0 |      56 | GRCh37.p10 |  63785537 | 63785593 | 56M       |
    ### def test_transcriptmapper_TranscriptMapper_GPHB5(self):
    ###     ac = 'NM_145171.3'
    ###     tm = TranscriptMapper(self.db,ac,self.ref)
    ###     pass

    def run_test_cases(self, tm, test_cases):
        for test_case in test_cases:
            g = hgvs.location.Interval(start=test_case['gs'], end=test_case['ge'])
            r = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['rs'], offset=test_case['so'], datum=test_case['d']),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['re'], offset=test_case['eo'], datum=test_case['d']))
            c = hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=test_case['cs'], offset=test_case['so'], datum=test_case['d'] + 1),
                    end=hgvs.location.BaseOffsetPosition(base=test_case['ce'], offset=test_case['eo'], datum=test_case['d'] + 1))
            try:
                if test_case['type'] == 'g' or test_case['type'] == 'grc':
                    self.assertEquals(tm.hgvsg_to_hgvsr(g), r)
                if test_case['type'] == 'grc':
                    self.assertEquals(tm.hgvsr_to_hgvsg(r), g)
                    self.assertEquals(tm.hgvsr_to_hgvsc(r), c)
                    self.assertEquals(tm.hgvsc_to_hgvsr(c), r)
                    self.assertEquals(tm.hgvsg_to_hgvsc(g), c)
                    self.assertEquals(tm.hgvsc_to_hgvsg(c), g)
            except Exception, msg:
                print 'FAIL: test case: %s' % test_case
                print msg


if __name__ == '__main__':
    unittest.main()
