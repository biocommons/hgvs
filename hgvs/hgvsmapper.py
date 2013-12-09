import os,re,sys

from Bio.Seq import Seq
import recordtype

import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.stopgap
import hgvs.transcriptmapper
from hgvs.utils import reverse_complement
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder
import hgvs.variant

class HGVSMapper(object):
    """
    Maps HGVS variants to and from g., r., c., and p. representations.
    All methods require and return objects of type hgvs.variant.Variant.
    """

    def __init__(self,db=None,cache_transcripts=False):
        self.db = db
        self.cache_transcripts = cache_transcripts
        self.__tm_cache = {}

    def hgvsg_to_hgvsc(self,var_g,ac,ref='GRCh37.p10'):
        """Given a genomic (g.) HGVS variant, return a transcript (c.) variant on the specified transcript.
        hgvs must be an HGVS-formatted variant or variant position.
        """

        if not (var_g.type == 'g'):
            raise hgvs.exception.InvalidHGVSVariantError('Expected a genomic (g.); got '+var_g)

        tm = self._fetch_TranscriptMapper(ac=ac,ref=ref)
        
        pos_c = tm.hgvsg_to_hgvsc( var_g.posedit.pos )
        if isinstance(var_g.posedit.edit,hgvs.edit.NARefAlt):
            if tm.strand == 1:
                edit_c = var_g.posedit.edit
            else:
                edit_g = var_g.posedit.edit
                edit_c = hgvs.edit.NARefAlt(
                    ref = reverse_complement(edit_g.ref),
                    alt = reverse_complement(edit_g.alt),
                    )
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        var_c = hgvs.variant.SequenceVariant(ac=ac,
                                             type='c',
                                             posedit=hgvs.posedit.PosEdit( pos_c, edit_c ) )
        return var_c

    def hgvsc_to_hgvsp(self, var_c):
        """Convert hgvsc tag to hgvsp tag

        :param var_c: hgvsc tag
        :type SequenceVariant
        :return hgvsp tag
        :type SequenceVariant
        """

        # TODO - error handling - invalid hgvs tag
        # TODO - error handling - transcript not found
        # TODO - error handling - unsupported hgvs tag transformation

        class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                                      ['transcript_sequence', 'aa_sequence',
                                                       'cds_start', 'cds_stop', 'protein_accession'])):

            @classmethod
            def setup_transcript_data(cls, ac, db, ref='GRCh37.p10'):
                """helper for generating RefTranscriptData from for hgvsc_to_hgvsp"""
                tx_info = db.get_tx_info(ac)
                tx_seq = db.get_tx_seq(ac)

                # use 1-based hgvs coords
                cds_start = tx_info['cds_start_i'] + 1
                cds_stop = tx_info['cds_stop_i']

                tx_seq_cds = Seq(tx_seq[cds_start - 1:cds_stop])
                protein_seq = str(tx_seq_cds.translate())
                protein_acc = hgvs.stopgap.pseq_to_ac(protein_seq)

                transcript_data = RefTranscriptData(tx_seq, protein_seq, cds_start,
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

    ############################################################################
    ## Internal methods

    def _fetch_TranscriptMapper(self,ac,ref='GRCh37.p10'):
        """
        Get a new TranscriptMapper for the given transcript accession (ac),
        possibly caching the result.
        """
        try:
            tm = self.__tm_cache[ac]
        except KeyError:
            tm = hgvs.transcriptmapper.TranscriptMapper(self.db, ref = ref, ac = ac)
            if self.cache_transcripts:
                self.__tm_cache[ac] = tm
        return tm
 

if __name__ == '__main__':
    import uta.db.transcriptdb
    import hgvs.parser

    hp = hgvs.parser.Parser()
    uta_conn = uta.db.transcriptdb.TranscriptDB()
    hm = HGVSMapper(uta_conn, cache_transcripts=True)

    # From garcia.tsv:
    # AOAH    NM_001177507.1:c.1486G>A      
    hgvs_g = 'NC_000007.13:g.36561662C>T'
    hgvs_c = 'NM_001637.3:c.1582G>A'
    
    var_g = hp.parse_hgvs_variant(hgvs_g)
    var_c = hm.hgvsg_to_hgvsc( var_g, 'NM_001637.3' )

    print( str(var_g) + ' --> ' + str(var_c) )
