import copy

from Bio.Seq import Seq
import recordtype

import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.transcriptmapper
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder
import hgvs.variant

from hgvs.utils import reverse_complement
from hgvs.utils.deprecated import deprecated



class VariantMapper(object):
    """
    Maps HGVS variants to and from g., r., c., and p. representations.
    All methods require and return objects of type :class:`hgvs.variant.SequenceVariant`.
    """

    def __init__(self,hdp,cache_transcripts=False):
        self.hdp = hdp
        self.cache_transcripts = cache_transcripts
        self.__tm_cache = {}


    def g_to_c(self, var_g, tx_ac, alt_aln_method='splign'):
        """Given a genomic (g.) parsed HGVS variant, return a transcript (c.) variant on the specified transcript using
        the specified alignment method (default is 'splign' from NCBI).

        :param hgvs.variant.SequenceVariant var_g: a variant object
        :param str tx_ac: a transcript accession (e.g., NM_012345.6 or ENST012345678)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`) using CDS coordinates
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_g is not of type 'g'

        """

        if not (var_g.type == 'g'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a genomic (g.) variant; got '+ str(var_g))

        tm = self._fetch_TranscriptMapper(tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method=alt_aln_method)
        
        pos_c = tm.g_to_c( var_g.posedit.pos )
        edit_c = self._convert_edit_check_strand(tm.strand, var_g.posedit.edit)
        var_c = hgvs.variant.SequenceVariant(ac=tx_ac,
                                             type='c',
                                             posedit=hgvs.posedit.PosEdit( pos_c, edit_c ) )
        return var_c


    def g_to_r(self, var_g, tx_ac, alt_aln_method='splign'):
        """Given a genomic (g.) parsed HGVS variant, return a transcript (r.) variant on the specified transcript using
        the specified alignment method (default is 'splign' from NCBI).

        :param hgvs.variant.SequenceVariant var_g: a variant object
        :param str tx_ac: a transcript accession (e.g., NM_012345.6 or ENST012345678)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`) using transcript (r.) coordinates
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_g is not of type 'g'

        """

        if not (var_g.type == 'g'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a genomic (g.); got '+ str(var_g))

        tm = self._fetch_TranscriptMapper(tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method=alt_aln_method)

        pos_r = tm.g_to_r( var_g.posedit.pos )
        edit_r = self._convert_edit_check_strand(tm.strand, var_g.posedit.edit)
        var_r = hgvs.variant.SequenceVariant(ac=tx_ac,
                                             type='r',
                                             posedit=hgvs.posedit.PosEdit( pos_r, edit_r ) )
        return var_r


    def r_to_g(self, var_r, alt_ac, alt_aln_method='splign'):
        """Given an RNA (r.) parsed HGVS variant, return a genomic (g.) variant on the specified transcript using
        the specified alignment method (default is 'splign' from NCBI).

        :param hgvs.variant.SequenceVariant var_r: a variant object
        :param str alt_ac: a reference sequence accession (e.g., NC_000001.11)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`)
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_r is not of type 'r'

        """

        if not (var_r.type == 'r'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a RNA (r.); got '+ str(var_r))

        tm = self._fetch_TranscriptMapper(tx_ac=var_r.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method)

        pos_g = tm.r_to_g( var_r.posedit.pos )
        edit_g = self._convert_edit_check_strand(tm.strand, var_r.posedit.edit)

        var_g = hgvs.variant.SequenceVariant(ac=alt_ac,
                                             type='g',
                                             posedit=hgvs.posedit.PosEdit( pos_g, edit_g ) )
        return var_g

    
    def c_to_g(self, var_c, alt_ac, alt_aln_method='splign'):
        """Given a cDNA (c.) parsed HGVS variant, return a genomic (g.) variant on the specified transcript using
        the specified alignment method (default is 'splign' from NCBI).

        :param hgvs.variant.SequenceVariant var_c: a variant object
        :param str alt_ac: a reference sequence accession (e.g., NC_000001.11)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`)
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_c is not of type 'c'

        """

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a cDNA (c.); got ' + str(var_c))

        tm = self._fetch_TranscriptMapper(tx_ac=var_c.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method)

        pos_g = tm.c_to_g(var_c.posedit.pos)
        edit_g = self._convert_edit_check_strand(tm.strand, var_c.posedit.edit)

        var_g = hgvs.variant.SequenceVariant(ac=alt_ac,
                                             type='g',
                                             posedit=hgvs.posedit.PosEdit(pos_g, edit_g))
        return var_g


    def c_to_r(self, var_c):
        """Given a cDNA (c.) parsed HGVS variant, return a RNA (r.) variant on the specified transcript using
        the specified alignment method (default is 'transcript' indicating a self alignment).

        :param hgvs.variant.SequenceVariant var_c: a variant object
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`)
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_c is not of type 'c'

        """

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a cDNA (c.); got ' + str(var_c))

        tm = self._fetch_TranscriptMapper(tx_ac=var_c.ac, alt_ac=var_c.ac, alt_aln_method='transcript')
        pos_r = tm.c_to_r(var_c.posedit.pos)

        # not necessary to check strand
        if isinstance(var_c.posedit.edit, hgvs.edit.NARefAlt) or isinstance(var_c.posedit.edit, hgvs.edit.Dup):
            edit_r = var_c.posedit.edit
        else:
            raise NotImplemented('Only NARefAlt/Dup types are currently implemented')

        var_r = hgvs.variant.SequenceVariant(ac=var_c.ac,
                                             type='r',
                                             posedit=hgvs.posedit.PosEdit( pos_r, edit_r ) )
        return var_r


    def r_to_c(self, var_r):
        """Given an RNA (r.) parsed HGVS variant, return a cDNA (c.) variant on the specified transcript using
        the specified alignment method (default is 'transcript' indicating a self alignment).

        :param hgvs.variant.SequenceVariant var_r: a variant object
        :returns: variant object (:class:`hgvs.variant.SequenceVariant`)
        :raises hgvs.exceptions.HGVSInvalidVariantError: if var_r is not of type 'r'

        """

        if not (var_r.type == 'r'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected RNA (r.); got ' + str(var_r))

        tm = self._fetch_TranscriptMapper(tx_ac=var_r.ac, alt_ac=var_r.ac, alt_aln_method='transcript')
        pos_c = tm.r_to_c(var_r.posedit.pos)

        # not necessary to check strand
        if isinstance(var_r.posedit.edit, hgvs.edit.NARefAlt) or isinstance(var_r.posedit.edit, hgvs.edit.Dup):
            edit_c = var_r.posedit.edit
        else:
            raise NotImplemented('Only NARefAlt types are currently implemented')

        var_c = hgvs.variant.SequenceVariant(ac=var_r.ac,
                                             type='c',
                                             posedit=hgvs.posedit.PosEdit( pos_c, edit_c ) )
        return var_c


    # TODO: c_to_p needs refactoring
    # TODO: data prep belongs in the data interface
    def c_to_p(self, var_c, pro_ac=None):
        """
        Converts a c. SequenceVariant to a p. SequenceVariant on the specified protein accession
        Author: Rudy Rico

        :param SequenceVariant var_c: hgvsc tag
        :param str pro_ac: protein accession
        :rtype: hgvs.variant.SequenceVariant

        """

        class RefTranscriptData(recordtype.recordtype('RefTranscriptData',
                                                      ['transcript_sequence', 'aa_sequence',
                                                       'cds_start', 'cds_stop', 'protein_accession'])):

            @classmethod
            def setup_transcript_data(cls, hdp, tx_ac, pro_ac):
                """helper for generating RefTranscriptData from for c_to_p"""
                tx_info = hdp.get_tx_identity_info(var_c.ac)
                tx_seq = hdp.get_tx_seq(tx_ac)

                if tx_info is None or tx_seq is None:
                    raise hgvs.exceptions.HGVSError("Missing transcript data for accession: {}".format(tx_ac))

                # use 1-based hgvs coords
                cds_start = tx_info['cds_start_i'] + 1
                cds_stop = tx_info['cds_end_i']

                # padding list so biopython won't complain during the conversion
                tx_seq_to_translate = tx_seq[cds_start - 1:cds_stop]
                if len(tx_seq_to_translate) % 3 != 0:
                    ''.join(list(tx_seq_to_translate).extend(['N']*((3-len(tx_seq_to_translate) % 3) % 3)))

                tx_seq_cds = Seq(tx_seq_to_translate)
                protein_seq = str(tx_seq_cds.translate())

                if pro_ac is None:
                    # get_acs... will always return at least the MD5_ accession
                    pro_ac = hdp.get_acs_for_protein_seq(protein_seq)[0]

                transcript_data = RefTranscriptData(tx_seq, protein_seq, cds_start,
                                                    cds_stop, pro_ac)

                return transcript_data

        if not (var_c.type == 'c'):
            raise hgvs.exceptions.HGVSInvalidVariantError('Expected a cDNA (c.); got ' + str(var_c))

        reference_data = RefTranscriptData.setup_transcript_data(self.hdp, var_c.ac, pro_ac)
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
    ## DEPRECATED METHODS
    @deprecated(use_instead='g_to_c(...)')
    def hgvsg_to_hgvsc(self,*args,**kwargs):
        return self.g_to_c(*args,**kwargs)
    @deprecated(use_instead='g_to_r(...)')
    def hgvsg_to_hgvsr(self,*args,**kwargs):
        return self.g_to_r(*args,**kwargs)
    @deprecated(use_instead='r_to_g(...)')
    def hgvsr_to_hgvsg(self,*args,**kwargs):
        return self.r_to_g(*args,**kwargs)
    @deprecated(use_instead='c_to_g(...)')
    def hgvsc_to_hgvsg(self,*args,**kwargs):
        return self.c_to_g(*args,**kwargs)
    @deprecated(use_instead='c_to_r(...)')
    def hgvsc_to_hgvsr(self,*args,**kwargs):
        return self.c_to_r(*args,**kwargs)
    @deprecated(use_instead='r_to_c(...)')
    def hgvsr_to_hgvsc(self,*args,**kwargs):
        return self.r_to_c(*args,**kwargs)
    @deprecated(use_instead='c_to_p(...)')
    def hgvsc_to_hgvsp(self,*args,**kwargs):
        return self.c_to_p(*args,**kwargs)


    ############################################################################
    ## Internal methods

    def _fetch_TranscriptMapper(self, tx_ac, alt_ac, alt_aln_method):
        """
        Get a new TranscriptMapper for the given transcript accession (ac),
        possibly caching the result.
        """
        key = '_'.join([tx_ac, alt_ac, alt_aln_method])
        try:
            tm = self.__tm_cache[key]
        except KeyError:
            tm = hgvs.transcriptmapper.TranscriptMapper(self.hdp, tx_ac=tx_ac, alt_ac=alt_ac,
                                                        alt_aln_method=alt_aln_method)
            if self.cache_transcripts:
                self.__tm_cache[key] = tm
        return tm


    @staticmethod
    def _convert_edit_check_strand(strand, edit_in):
        """
        Convert an edit from one type to another, based on the stand and type
        """
        if isinstance(edit_in, hgvs.edit.NARefAlt):
            if strand == 1:
                edit_out = copy.deepcopy(edit_in)
            else:
                try:
                    # if smells like an int, do nothing
                    # TODO: should use ref_n, right?
                    int(edit_in.ref)
                    ref = edit_in.ref
                except (ValueError, TypeError):
                    ref = reverse_complement(edit_in.ref)
                edit_out = hgvs.edit.NARefAlt(
                    ref = ref,
                    alt = reverse_complement(edit_in.alt),
                    )
        elif isinstance(edit_in, hgvs.edit.Dup):
            if strand == 1:
                edit_out = copy.deepcopy(edit_in)
            else:
                edit_out = hgvs.edit.Dup(
                    seq = reverse_complement(edit_in.seq)
                )
        else:
            raise NotImplemented('Only NARefAlt/Dup types are currently implemented')
        return edit_out


# TODO: move assembly data to UTA and query through HDPI
primary_assembly_accessions = {
    'GRCh37': {
        'NC_000001.10', 'NC_000002.11', 'NC_000003.11',
        'NC_000004.11', 'NC_000005.9' , 'NC_000006.11',
        'NC_000007.13', 'NC_000008.10', 'NC_000009.11',
        'NC_000010.10', 'NC_000011.9' , 'NC_000012.11',
        'NC_000013.10', 'NC_000014.10', 'NC_000015.9',
        'NC_000016.9' , 'NC_000017.10', 'NC_000018.9',
        'NC_000019.9' , 'NC_000020.10', 'NC_000021.8',
        'NC_000022.10', 'NC_000023.10', 'NC_000024.9',
        },
    }

class EasyVariantMapper(VariantMapper):
    """Provides simplified variant mapping for a single assembly and
    transcript-reference alignment method.
    
    EasyVariantMapper is instantiated with a primary_assembly and
    alt_aln_method. These enable the following conveniences over
    VariantMapper:

    * The primary assembly and alignment method are used to
      automatically select an appropriate chromosomal reference
      sequence when mapping from a transcript to a genome (i.e.,
      c_to_g(...) and r_to_g(...)).

    * A new method, relevant_trancripts(g_variant), returns a list of
      transcript accessions available for the specified variant. These
      accessions are candidates mapping from genomic to trancript
      coordinates (i.e., g_to_c(...) and g_to_r(...)).

    [tests occur in module doc (rather than in method doc) to use a
    single db connection]

    IMPORTANT: Callers should be prepared to catch HGVSError
    exceptions. These will be thrown whenever a transcript maps
    ambiguously to a chromosome, such as for pseudoautosomal region
    transcripts.
    """

    def __init__(self,hdp,primary_assembly='GRCh37',alt_aln_method='splign',cache_transcripts=False):
        super(EasyVariantMapper,self).__init__(hdp=hdp,cache_transcripts=cache_transcripts)
        self.primary_assembly = primary_assembly
        self.alt_aln_method = alt_aln_method
        self.primary_assembly_accessions = set(primary_assembly_accessions[primary_assembly])

    def g_to_c(self, var_g, tx_ac): 
        return super(EasyVariantMapper,self).g_to_c(var_g, tx_ac, alt_aln_method=self.alt_aln_method)
    def g_to_r(self, var_g, tx_ac): 
        return super(EasyVariantMapper,self).g_to_r(var_g, tx_ac, alt_aln_method=self.alt_aln_method)

    def c_to_g(self, var_c): 
        alt_ac = self._alt_ac_for_tx_ac(var_c.ac)
        return super(EasyVariantMapper,self).c_to_g(var_c, alt_ac, alt_aln_method=self.alt_aln_method)
    def r_to_g(self, var_r): 
        alt_ac = self._alt_ac_for_tx_ac(var_r.ac)
        return super(EasyVariantMapper,self).r_to_g(var_r, alt_ac, alt_aln_method=self.alt_aln_method)

    def c_to_r(self, var_c): 
        return super(EasyVariantMapper,self).c_to_r(var_c)
    def r_to_c(self, var_r): 
        return super(EasyVariantMapper,self).r_to_c(var_r)

    def c_to_p(self, var_c): 
        return super(EasyVariantMapper,self).c_to_p(var_c)

    def relevant_transcripts(self,var_g):
        """return list of transcripts accessions (strings) for given variant,
        selected by genomic overlap"""
        tx = self.hdp.get_tx_for_region(var_g.ac,
                                        self.alt_aln_method,
                                        var_g.posedit.pos.start.base,
                                        var_g.posedit.pos.end.base)
        return [ e['tx_ac'] for e in tx ]

    def _alt_ac_for_tx_ac(self,tx_ac):
        """return chromosomal accession for given transcript accession (and
        the primary_assembly and aln_method setting used to
        instantiate this EasyVariantMapper)

        """
        alt_acs = [e['alt_ac'] 
                   for e in  self.hdp.get_tx_mapping_options(tx_ac)
                   if e['alt_aln_method'] == self.alt_aln_method
                   and e['alt_ac'] in self.primary_assembly_accessions]
        if len(alt_acs) > 1:
            raise hgvs.exceptions.HGVSError("Multiple chromosomal alignments for {tx_ac} in {pa} using {am} (likely paralog or pseudoautosomal region)".format(
                tx_ac=tx_ac, pa=self.primary_assembly, am=self.alt_aln_method))
        if len(alt_acs) == 0:
            raise hgvs.exceptions.HGVSError("No alignments for {tx_ac} in {pa} using {am}".format(
                tx_ac=tx_ac, pa=self.primary_assembly, am=self.alt_aln_method))
        return alt_acs[0]       # exactly one remains



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
