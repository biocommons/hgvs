import os,re,sys

import hgvs.variant
import hgvs.posedit
import hgvs.location
import hgvs.transcriptmapper
from hgvs.utils import reverse_complement

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
            raise InvalidHGVSVariantError('Expected a genomic (g.); got '+var_g)

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
