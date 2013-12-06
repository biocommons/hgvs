#
# Translate hgvsc tag to hgvsp tag
#
from Bio.Seq import Seq
import recordtype

import hgvs.edit
import hgvs.parser
import hgvs.stopgap
import hgvs.utils
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder


class HgvsCToP(object):

    def __init__(self,db=None,cache_transcripts=False):
        self.db = db
        self.cache_transcripts = cache_transcripts
        self.hgvs_parser = hgvs.parser.Parser()
        self.tm_cache = {}

    def hgvsc_to_hgvsp(self, hgvsc):
        """Convert hgvsc tag to hgvsp tag

        :param hgvsc: hgvsc tag
        :type str
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
                ordered_keys = ord_seq.keys()
                ordered_keys.sort()
                transcript_sequence = ''.join(ord_seq[x] for x in ordered_keys)

                transcript_dna_cds = Seq(transcript_sequence[cds_start - 1:cds_stop])
                protein_seq = str(transcript_dna_cds.translate())
                protein_acc = hgvs.stopgap.pseq_to_ac(protein_seq)

                transcript_data = RefTranscriptData(transcript_sequence, protein_seq, cds_start,
                                                    cds_stop, protein_acc)

                return transcript_data


        var = self.hgvs_parser.parse_hgvs_variant(hgvsc)
        reference_data = RefTranscriptData.setup_transcript_data(var.seqref, self.db)
        builder = altseqbuilder.AltSeqBuilder(var, reference_data)

        # TODO - handle case where you get 2+ alt sequences back; currently get list of 1 element
        # loop structure implemented to handle this, but doesn't really do anything currently.
        alt_data = builder.build_altseq()

        hgvsps = []
        for alt in alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data.aa_sequence,
                                                    alt.aa_sequence,
                                                    reference_data.protein_accession,
                                                    alt.frameshift_start
                                                    )
            hgvsp = builder.build_hgvsp()
            hgvsps.append(hgvsp)

        hgvsp = hgvsps[0]

        return str(hgvsp)

