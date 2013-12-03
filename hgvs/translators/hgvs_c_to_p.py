#
# Translate hgvsc tag to hgvsp tag
#
import pprint

from Bio.Seq import Seq
import recordtype

import hgvs.edit
import hgvs.parser
import hgvs.utils
import hgvs.translators.tools.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.translators.tools.altseqbuilder as altseqbuilder

DBG = True  # DEBUG flag

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
        """

        # translate hgvsc
        var = self.parser.parse_hgvs_variant(hgvsc)
        if DBG:
            print "START DEBUG"
            print hgvsc

        # get transcript from datasource & convert to AA
        reference_data = self._setup_transcript_data(var.seqref)

        # incorporate DNA change into sequence(s) and translate
        builder = altseqbuilder.AltSeqBuilder(var, reference_data)
        alt_data = builder.build_altseq()

        if DBG:
            pprint.pprint(reference_data)
            pprint.pprint(alt_data)

        # perform comparison to get hgvs tag
        hgvsps = []
        for alt in alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data.aa_sequence,
                                                    alt.aa_sequence, alt.frameshift_start)
            hgvsp_no_acc = builder.build_hgvsp()

            hgvsps.append(self._add_accession(hgvsp_no_acc, alt.protein_accession))

        hgvsp = hgvsps[0]    # TODO - handle case of more than 1 variant

        if DBG:
            print "END DEBUG"

        return hgvsp

    #
    # internal methods
    #

    def _setup_transcript_data(self, ac):
        """wrapper to convert transcript data from external source into common form
        :param ac accession #
        :type str
        :return transcript info
        :type dict
        """
        # TODO - take whatever form you get the input from & convert to a uniform format
        # TODO - handle what UTA returns & convert test inputs to mimic that
        td = self.datasource.fetch_transcript_exons(ac)

        transcript_dna_cds = Seq(td['transcript_sequence'][td['cds_start'] - 1:td['cds_stop']])

        transcript_data = RefTranscriptData(td['transcript_sequence'],
                                            str(transcript_dna_cds.translate()),
                                            td['cds_start'],
                                            td['cds_stop'],
                                            td['protein_accession'])

        return transcript_data

    def _add_accession(self, hgvsps, acc):
        """Helper to add protein accession to hgvsp
        """
        for hgvsp in hgvsps:
            hgvsp.ac = acc

        if len(hgvsps) > 1:
            raise NotImplementedError("Can only handle single variants")
        else:
            hgvsp = hgvsps[0]

        return hgvsp



def main():
    pass


if __name__ == "__main__":
    main()