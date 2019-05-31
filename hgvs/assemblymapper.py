# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

import logging

import hgvs
import hgvs.normalizer

from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSInvalidVariantError, HGVSUnsupportedOperationError
from hgvs.variantmapper import VariantMapper

_logger = logging.getLogger(__name__)


class AssemblyMapper(VariantMapper):
    """Provides simplified variant mapping for a single assembly and
    transcript-reference alignment method.

    AssemblyMapper is instantiated with an assembly name and
    alt_aln_method. These enable the following conveniences over
    VariantMapper:

    * The assembly and alignment method are used to
      automatically select an appropriate chromosomal reference
      sequence when mapping from a transcript to a genome (i.e.,
      c_to_g(...) and n_to_g(...)).

    * A new method, relevant_trancripts(g_variant), returns a list of
      transcript accessions available for the specified variant. These
      accessions are candidates mapping from genomic to trancript
      coordinates (i.e., g_to_c(...) and g_to_n(...)).

    Note: AssemblyMapper supports only chromosomal references (e.g.,
    NC_000006.11). It does not support contigs or other genomic
    sequences (e.g., NT_167249.1).

    """

    def __init__(self,
                 hdp,
                 assembly_name=hgvs.global_config.mapping.assembly,
                 alt_aln_method=hgvs.global_config.mapping.alt_aln_method,
                 normalize=hgvs.global_config.mapping.normalize,
                 prevalidation_level=hgvs.global_config.mapping.prevalidation_level,
                 in_par_assume=hgvs.global_config.mapping.in_par_assume,
                 replace_reference=hgvs.global_config.mapping.replace_reference,
                 add_gene_symbol=hgvs.global_config.mapping.add_gene_symbol,
                 *args,
                 **kwargs):
        """
        :param object hdp: instance of hgvs.dataprovider subclass
        :param bool replace_reference: replace reference (entails additional network access)

        :param str assembly_name: name of assembly ("GRCh38.p5")
        :param str alt_aln_method: genome-transcript alignment method ("splign", "blat", "genewise")
        :param bool normalize: normalize variants
        :param str prevalidation_level: None or Intrinsic or Extrinsic validation before mapping
        :param str in_par_assume: during x_to_g, assume this chromosome name if alignment is ambiguous

        :raises HGVSError subclasses: for a variety of mapping and data lookup failures
        """

        super(AssemblyMapper, self).__init__(
            hdp=hdp,
            replace_reference=replace_reference,
            prevalidation_level=prevalidation_level,
            add_gene_symbol=add_gene_symbol,
            *args,
            **kwargs)
        self.assembly_name = assembly_name
        self.alt_aln_method = alt_aln_method
        self.normalize = normalize
        self.in_par_assume = in_par_assume
        self._norm = None
        if self.normalize:
            self._norm = hgvs.normalizer.Normalizer(
                hdp, alt_aln_method=alt_aln_method, validate=False)
        self._assembly_map = {
            k: v
            for k, v in hdp.get_assembly_map(self.assembly_name).items() if k.startswith("NC_")
        }
        self._assembly_accessions = set(self._assembly_map.keys())

    def __repr__(self):
        return ("{self.__module__}.{t.__name__}(alt_aln_method={self.alt_aln_method}, "
                "assembly_name={self.assembly_name}, normalize={self.normalize}, "
                "prevalidation_level={self.prevalidation_level}, "
                "replace_reference={self.replace_reference})".format(self=self, t=type(self)))

    def g_to_c(self, var_g, tx_ac):
        var_out = super(AssemblyMapper, self).g_to_c(
            var_g, tx_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def g_to_n(self, var_g, tx_ac):
        var_out = super(AssemblyMapper, self).g_to_n(
            var_g, tx_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def g_to_t(self, var_g, tx_ac):
        var_out = super(AssemblyMapper, self).g_to_t(
            var_g, tx_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def c_to_g(self, var_c):
        alt_ac = self._alt_ac_for_tx_ac(var_c.ac)
        var_out = super(AssemblyMapper, self).c_to_g(
            var_c, alt_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def n_to_g(self, var_n):
        alt_ac = self._alt_ac_for_tx_ac(var_n.ac)
        var_out = super(AssemblyMapper, self).n_to_g(
            var_n, alt_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def t_to_g(self, var_t):
        alt_ac = self._alt_ac_for_tx_ac(var_t.ac)
        var_out = super(AssemblyMapper, self).t_to_g(
            var_t, alt_ac, alt_aln_method=self.alt_aln_method)
        return self._maybe_normalize(var_out)

    def t_to_p(self, var_t):
        """Return a protein variant, or "non-coding" for non-coding variant types

        CAUTION: Unlike other x_to_y methods that always return
        SequenceVariant instances, this method returns a string when
        the variant type is ``n``.  This is intended as a convenience,
        particularly when looping over ``relevant_transcripts``,
        projecting with ``g_to_t``, then desiring a protein
        representation for coding transcripts.

        """
        if var_t.type == "n":
            return "non-coding"
        if var_t.type == "c":
            return self.c_to_p(var_t)
        raise HGVSInvalidVariantError("Expected a coding (c.) or non-coding (n.) variant; got " +
                                      str(var_t))

    def c_to_n(self, var_c):
        var_out = super(AssemblyMapper, self).c_to_n(var_c)
        return self._maybe_normalize(var_out)

    def n_to_c(self, var_n):
        var_out = super(AssemblyMapper, self).n_to_c(var_n)
        return self._maybe_normalize(var_out)

    def c_to_p(self, var_c):
        var_out = super(AssemblyMapper, self).c_to_p(var_c)
        return self._maybe_normalize(var_out)

    def relevant_transcripts(self, var_g):
        """return list of transcripts accessions (strings) for given variant,
        selected by genomic overlap"""
        tx = self.hdp.get_tx_for_region(var_g.ac, self.alt_aln_method, var_g.posedit.pos.start.base,
                                        var_g.posedit.pos.end.base)
        return [e["tx_ac"] for e in tx]

    def _alt_ac_for_tx_ac(self, tx_ac):
        """return chromosomal accession for given transcript accession (and
        the_assembly and aln_method setting used to instantiate this
        AssemblyMapper)

        """
        alt_acs = [
            e["alt_ac"] for e in self.hdp.get_tx_mapping_options(tx_ac) if
            e["alt_aln_method"] == self.alt_aln_method and e["alt_ac"] in self._assembly_accessions
        ]

        if not alt_acs:
            raise HGVSDataNotAvailableError("No alignments for {tx_ac} in {an} using {am}".format(
                tx_ac=tx_ac, an=self.assembly_name, am=self.alt_aln_method))

        # TODO: conditional is unnecessary; remove
        if len(alt_acs) > 1:
            names = set(self._assembly_map[ac] for ac in alt_acs)
            if names != set("XY"):
                alts = ", ".join(
                    ["{ac} ({n})".format(ac=ac, n=self._assembly_map[ac]) for ac in alt_acs])
                raise HGVSError("Multiple chromosomal alignments for {tx_ac} in {an}"
                                " using {am} (non-pseudoautosomal region) [{alts}]".format(
                                    tx_ac=tx_ac,
                                    an=self.assembly_name,
                                    am=self.alt_aln_method,
                                    alts=alts))

            # assume PAR
            if self.in_par_assume is None:
                raise HGVSError("Multiple chromosomal alignments for {tx_ac} in {an}"
                                " using {am} (likely pseudoautosomal region)".format(
                                    tx_ac=tx_ac, an=self.assembly_name, am=self.alt_aln_method))

            alt_acs = [ac for ac in alt_acs if self._assembly_map[ac] == self.in_par_assume]
            if len(alt_acs) != 1:
                raise HGVSError("Multiple chromosomal alignments for {tx_ac} in {an}"
                                " using {am}; in_par_assume={ipa} selected {n} of them".format(
                                    tx_ac=tx_ac,
                                    an=self.assembly_name,
                                    am=self.alt_aln_method,
                                    ipa=self.in_par_assume,
                                    n=len(alt_acs)))

        assert len(alt_acs) == 1, "Should have exactly one alignment at this point"
        return alt_acs[0]

    def _fetch_AlignmentMapper(self, tx_ac, alt_ac=None, alt_aln_method=None):
        """convenience version of VariantMapper._fetch_AlignmentMapper that
        derives alt_ac from transcript, assembly, and alt_aln_method
        used to instantiate the AssemblyMapper instance

        """

        if alt_ac is None:
            alt_ac = self._alt_ac_for_tx_ac(tx_ac)
        if alt_aln_method is None:
            alt_aln_method = self.alt_aln_method
        return super(AssemblyMapper, self)._fetch_AlignmentMapper(tx_ac, alt_ac, alt_aln_method)

    def _maybe_normalize(self, var):
        """normalize variant if requested, and ignore HGVSUnsupportedOperationError
        This is better than checking whether the variant is intronic because
        future UTAs will support LRG, which will enable checking intronic variants.
        """
        if self.normalize:
            try:
                return self._norm.normalize(var)
            except HGVSUnsupportedOperationError as e:
                _logger.warning(str(e) + "; returning unnormalized variant")
                # fall through to return unnormalized variant
        return var


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
