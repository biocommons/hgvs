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


class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                              ['transcript_sequence', 'aa_sequence',
                                               'cds_start', 'cds_stop', 'protein_accession'])):
    pass


class HgvsCToP(object):

    def __init__(self, datasource):
        """Constructor

        :param datasource
        :type base_input_source
        """
        self.datasource = datasource
        self.parser = hgvs.parser.Parser()

    def convert(self, hgvsc):
        """Convert hgvsc tag to hgvsp tag

        :param hgvsc: hgvsc tag
        :type str
        :return hgvsp tag
        :type SequenceVariant
        """

        # translate hgvsc
        var = self.parser.parse_hgvs_variant(hgvsc)

        # get transcript from datasource & convert to AA
        reference_data = self._setup_transcript_data(var.seqref)

        # incorporate DNA change into sequence(s) and translate
        builder = altseqbuilder.AltSeqBuilder(var, reference_data)
        alt_data = builder.build_altseq()

        # perform comparison to get hgvs tag
        hgvsps = []
        for alt in alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data.aa_sequence,
                                                    alt.aa_sequence,
                                                    reference_data.protein_accession,
                                                    alt.frameshift_start
                                                    )
            hgvsp = builder.build_hgvsp()
            hgvsps.append(hgvsp)

        hgvsp = hgvsps[0]    # TODO - handle case of more than 1 variant

        return hgvsp

    #
    # internal methods
    #

    def _setup_transcript_data(self, ac):
        """wrapper to convert transcript data from external source into common form

        :param ac accession #
        :type str
        :return transcript info
        :type recordtype
        """
        ts_exons = self.datasource.fetch_transcript_exons(ac, 'GRCh37.p10')
        ts_info = self.datasource.fetch_transcript_info(ac)

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

        transcript_data = RefTranscriptData(transcript_sequence,
                                            protein_seq,
                                            cds_start,
                                            cds_stop,
                                            protein_acc)

        return transcript_data