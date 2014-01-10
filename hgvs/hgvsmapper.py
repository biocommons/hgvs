"""
hgvs.hgvsmapper
"""

import os,re,sys

from Bio.Seq import Seq
import recordtype

import hgvs.exceptions
import hgvs.variant
import hgvs.posedit
import hgvs.location
import hgvs.transcriptmapper
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder
from hgvs.utils import reverse_complement
from hgvs.utils import chr_to_nc

class HGVSMapper(object):
    """
    Maps HGVS variants to and from g., r., c., and p. representations.
    All methods require and return objects of type hgvs.variant.Variant.
    """

    def __init__(self,bdi,cache_transcripts=False):
        self.bdi = bdi
        self.cache_transcripts = cache_transcripts
        self.__tm_cache = {}

    def hgvsg_to_hgvsc(self,var_g, ac, ref='GRCh37.p10'):
        """Given a genomic (g.) HGVS variant, return a transcript (c.) variant on the specified transcript.
        hgvs must be an HGVS-formatted variant or variant position.
        """

        if not (var_g.type == 'g'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a genomic (g.); got '+ str(var_g))

        tm = self._fetch_TranscriptMapper(ac=ac,ref=ref)
        
        pos_c = tm.hgvsg_to_hgvsc( var_g.posedit.pos )
        if isinstance(var_g.posedit.edit, hgvs.edit.NARefAlt):
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

    def hgvsg_to_hgvsr(self,var_g, ac, ref='GRCh37.p10'):
        """Given a genomic (g.) HGVS variant, return a transcript (r.) variant on the specified transcript.
        hgvs must be an HGVS-formatted variant or variant position.
        """

        if not (var_g.type == 'g'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a genomic (g.); got '+ str(var_g))

        tm = self._fetch_TranscriptMapper(ac=ac,ref=ref)

        pos_r = tm.hgvsg_to_hgvsr( var_g.posedit.pos )
        if isinstance(var_g.posedit.edit, hgvs.edit.NARefAlt):
            if tm.strand == 1:
                edit_r = var_g.posedit.edit
            else:
                edit_g = var_g.posedit.edit
                edit_r = hgvs.edit.NARefAlt(
                    ref = reverse_complement(edit_g.ref),
                    alt = reverse_complement(edit_g.alt),
                    )
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        var_r = hgvs.variant.SequenceVariant(ac=ac,
                                             type='r',
                                             posedit=hgvs.posedit.PosEdit( pos_r, edit_r ) )
        return var_r

    def hgvsr_to_hgvsg(self,var_r, ref='GRCh37.p10'):
        """Given an RNA (r.) HGVS variant, return a genomic (g.) variant from the inferred transcript.
        hgvs must be an HGVS-formatted variant or variant position.
        """

        if not (var_r.type == 'r'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a RNA (r.); got '+ str(var_r))

        tm = self._fetch_TranscriptMapper(ac=var_r.ac,ref=ref)

        pos_g = tm.hgvsr_to_hgvsg( var_r.posedit.pos )
        if isinstance(var_r.posedit.edit, hgvs.edit.NARefAlt):
            if tm.strand == 1:
                edit_r = var_r.posedit.edit
            else:
                edit_r = var_r.posedit.edit
                edit_g = hgvs.edit.NARefAlt(
                    ref = reverse_complement(edit_r.ref),
                    alt = reverse_complement(edit_r.alt),
                    )
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        # get NC accession for g.
        var_g = hgvs.variant.SequenceVariant(ac=chr_to_nc(tm.tx_info['chr']),
                                             type='g',
                                             posedit=hgvs.posedit.PosEdit( pos_g, edit_g ) )
        return var_g
    
    def hgvsc_to_hgvsg(self, var_c, ref='GRCh37.p10'):
        """Given a cDNA (c.) HGVS variant and an inferred transcript, return a genomic (g.) variant.
        hgvs must be an HGVS-formatted variant or variant position.
        """

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a cDNA (c.); got ' + str(var_c))

        tm = self._fetch_TranscriptMapper(ac=var_c.ac, ref=ref)

        pos_g = tm.hgvsc_to_hgvsg(var_c.posedit.pos)
        if isinstance(var_c.posedit.edit, hgvs.edit.NARefAlt):
            if tm.strand == 1:
                edit_g = var_c.posedit.edit
            else:
                edit_c = var_c.posedit.edit
                edit_g = hgvs.edit.NARefAlt(
                    ref = reverse_complement(edit_c.ref),
                    alt = reverse_complement(edit_c.alt),
                )
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        # get NC accession for g.
        var_g = hgvs.variant.SequenceVariant(ac=chr_to_nc(tm.tx_info['chr']),
                                             type='g',
                                             posedit=hgvs.posedit.PosEdit(pos_g, edit_g))
        return var_g

    def hgvsc_to_hgvsr(self, var_c, ref='GRCh37.p10'):
        """Given a cDNA (c.) HGVS variant and an inferred transcript, return an RNA (r.) variant.
        hgvs must be an HGVS-formatted variant or variant position
        """

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a cDNA (c.); got ' + str(var_c))

        tm = self._fetch_TranscriptMapper(ac=var_c.ac, ref=ref)
        pos_r = tm.hgvsc_to_hgvsr(var_c.posedit.pos)

        # not necessary to check strand
        if isinstance(var_c.posedit.edit, hgvs.edit.NARefAlt):
            edit_r = var_c.posedit.edit
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        var_r = hgvs.variant.SequenceVariant(ac=var_c.ac,
                                             type='r',
                                             posedit=hgvs.posedit.PosEdit( pos_r, edit_r ) )
        return var_r

    def hgvsr_to_hgvsc(self, var_r, ref='CRCh37.p10'):
        """Given a RNA (r.) HGVS variant and an inferred transcript, return an cDNA 
        hgvs must be an HGVS-formatted variant or variant position
        """

        if not (var_r.type == 'r'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a cDNA (c.); got ' + str(var_r))

        tm = self._fetch_TranscriptMapper(ac=var_r.ac, ref=ref)
        pos_c = tm.hgvsr_to_hgvsc(var_r.posedit.pos)

        # not necessary to check strand
        if isinstance(var_r.posedit.edit, hgvs.edit.NARefAlt):
            edit_c = var_r.posedit.edit
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        var_c = hgvs.variant.SequenceVariant(ac=var_r.ac,
                                             type='c',
                                             posedit=hgvs.posedit.PosEdit( pos_c, edit_c ) )
        return var_c

    def hgvsc_to_hgvsp(self, var_c, ac_p):
        """
        Converts a c. SequenceVariant to a p. SequenceVariant on the specified protein accession

        :param var_c: hgvsc tag
        :type var_c: SequenceVariant
        :param ac_p: protein accession
        :type ac_p: str
        :rtype: SequenceVariant

        """

        class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                                      ['transcript_sequence', 'aa_sequence',
                                                       'cds_start', 'cds_stop', 'protein_accession'])):

            @classmethod
            def setup_transcript_data(cls, ac, ac_p, bdi, ref='GRCh37.p10'):
                """helper for generating RefTranscriptData from for hgvsc_to_hgvsp"""
                tx_info = bdi.get_tx_info(ac)
                tx_seq = bdi.get_tx_seq(ac)

                if tx_info is None or tx_seq is None:
                    raise hgvs.exceptions.HGVSError("Missing transcript data for accession: {}".format(ac))

                # use 1-based hgvs coords
                cds_start = tx_info['cds_start_i'] + 1
                cds_stop = tx_info['cds_end_i']

                # padding list so biopython won't complain during the conversion
                tx_seq_to_translate = tx_seq[cds_start - 1:cds_stop]
                if len(tx_seq_to_translate) % 3 != 0:
                    ''.join(list(tx_seq_to_translate).extend(['N']*((3-len(tx_seq_to_translate) % 3) % 3)))

                tx_seq_cds = Seq(tx_seq_to_translate)
                protein_seq = str(tx_seq_cds.translate())
                
                if ac_p is None:
                    # get_acs... will always return at least the MD5_ accession
                    ac_p = bdi.get_acs_for_protein_seq(protein_seq)[0]

                transcript_data = RefTranscriptData(tx_seq, protein_seq, cds_start,
                                                    cds_stop, ac_p)

                return transcript_data

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.InvalidHGVSVariantError('Expected a cDNA (c.); got ' + str(var_c))

        reference_data = RefTranscriptData.setup_transcript_data(var_c.ac, ac_p, self.bdi)
        builder = altseqbuilder.AltSeqBuilder(var_c, reference_data)

        # TODO - handle case where you get 2+ alt sequences back; currently get list of 1 element
        # loop structure implemented to handle this, but doesn't really do anything currently.
        all_alt_data = builder.build_altseq()

        var_ps = []
        for alt_data in all_alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data, alt_data)
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
            tm = hgvs.transcriptmapper.TranscriptMapper(self.bdi, ref = ref, ac = ac)
            if self.cache_transcripts:
                self.__tm_cache[ac] = tm
        return tm
 

if __name__ == '__main__':
    import hgvs.parser
    import hgvs.hgvsmapper
    from bdi.sources.uta0_sqlite import UTA0

    bdi = UTA0('/tmp/uta-0.0.4.db')
    hp = hgvs.parser.Parser()
    hm = hgvs.hgvsmapper.HGVSMapper(bdi, cache_transcripts=True)

    ref = 'GRCh37.p10'

    # From garcia.tsv:
    # AOAH    NM_001177507.1:c.1486G>A      
    hgvs_g = 'NC_000007.13:g.36561662C>T'
    hgvs_c = 'NM_001637.3:c.1582G>A'
    
    var_g = hp.parse_hgvs_variant(hgvs_g)
    var_c = hm.hgvsg_to_hgvsc( var_g, 'NM_001637.3' )

    print( str(var_g) + ' --> ' + str(var_c) )

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
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
