#
# Translate hgvsc tag to hgvsp tag
#
from Bio.Seq import Seq
import recordtype

import hgvs.edit
#import hgvs.parser
import hgvs.stopgap
import hgvs.utils
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder


class HgvsCToP(object):

    def __init__(self,db=None,cache_transcripts=False):
        self.db = db
        self.cache_transcripts = cache_transcripts
        #self.hgvs_parser = hgvs.parser.Parser()
        self.tm_cache = {}

    def hgvsc_to_hgvsp(self, var_c):
        """Convert hgvsc tag to hgvsp tag

        :param var_c: hgvsc tag
        :type SequenceVariant
        :return hgvsp tag
        :type SequenceVariant
        """

        class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                                      ['transcript_sequence', 'aa_sequence',
                                                       'cds_start', 'cds_stop', 'protein_accession'])):

            @classmethod
            def setup_transcript_data(cls, ac, db, ref='GRCh37.p10'):
                """helper for generating RefTranscriptData from for hgvsc_to_hgvsp"""
                ts_exons = db.get_tx_exons(ac)
                ts_info = db.get_tx_info(ac)

                # use 1-based hgvs coords
                cds_start = ts_info['cds_start_i'] + 1
                cds_stop = ts_info['cds_stop_i']

                # concatenate sequences by exon
                ord_seq = {x['ord']: x['t_seq_a'] for x in ts_exons}
                transcript_sequence = ''.join(ord_seq[x] for x in sorted(ord_seq.keys()))

                transcript_dna_cds = Seq(transcript_sequence[cds_start - 1:cds_stop])
                protein_seq = str(transcript_dna_cds.translate())
                protein_acc = hgvs.stopgap.pseq_to_ac(protein_seq)

                transcript_data = RefTranscriptData(transcript_sequence, protein_seq, cds_start,
                                                    cds_stop, protein_acc)

                return transcript_data


        reference_data = RefTranscriptData.setup_transcript_data(var_c.seqref, self.db)
        builder = altseqbuilder.AltSeqBuilder(var_c, reference_data)

        # TODO - handle case where you get 2+ alt sequences back; currently get list of 1 element
        # loop structure implemented to handle this, but doesn't really do anything currently.
        alt_data = builder.build_altseq()

        var_ps = []
        for alt in alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data.aa_sequence,
                                                    alt.aa_sequence,
                                                    reference_data.protein_accession,
                                                    alt.frameshift_start
                                                    )
            var_p = builder.build_hgvsp()
            var_ps.append(var_p)

        var_p = var_ps[0]

        return var_p

