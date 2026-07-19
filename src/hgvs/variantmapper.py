"""Projects variants between sequences using AlignmentMapper."""

import copy
import logging

from bioutils.sequences import TranslationTable, reverse_complement

import hgvs
import hgvs.alignmentmapper
import hgvs.dataproviders
import hgvs.dataproviders.interface
import hgvs.edit
import hgvs.location
import hgvs.normalizer
import hgvs.posedit
import hgvs.sequencevariant
import hgvs.validator
from hgvs.decorators.lru_cache import lru_cache
from hgvs.enums import PrevalidationLevel
from hgvs.exceptions import HGVSInvalidVariantError, HGVSUnsupportedOperationError
from hgvs.utils import altseq_to_hgvsp, altseqbuilder
from hgvs.utils.position import get_start_end, get_start_end_interbase
from hgvs.utils.reftranscriptdata import RefTranscriptData

_logger = logging.getLogger(__name__)


class VariantMapper:
    r"""Maps SequenceVariant objects between g., n., r., c., and p. representations.

    g⟷{c,n,r} projections are similar in that c, n, and r variants
    may use intronic coordinates. There are two essential differences
    that distinguish the three types:

    * Sequence start: In n and r variants, position 1 is the sequence
      start; in c variants, 1 is the transcription start site.
    * Alphabet: In n and c variants, sequences are DNA; in
      r. variants, sequences are RNA.

    This differences are summarized in this diagram::

      g ----acgtatgcac--gtctagacgt----      ----acgtatgcac--gtctagacgt----      ----acgtatgcac--gtctagacgt----
            \         \/         /              \         \/         /              \         \/         /
      c      acgtATGCACGTCTAGacgt         n      acgtatgcacgtctagacgt         r      acguaugcacgucuagacgu
                 1                               1                                   1
      p          MetHisValTer

    The g excerpt and exon structures are identical. The g⟷n
    transformation, which is the most basic, accounts for the offset
    of the aligned sequences (shown with "1") and the exon structure.
    The g⟷c transformation is akin to g⟷n transformation, but
    requires an addition offset to account for the translation start
    site (c.1).  The CDS in uppercase. The g⟷r transformation is
    akin to g⟷n transformation with a change of alphabet.

    Therefore, this this code uses g⟷n as the core transformation
    between genomic and c, n, and r variants: All c⟷g and r⟷g
    transformations use n⟷g after accounting for the above
    differences. For example, c_to_g accounts for the transcription
    start site offset, then calls n_to_g.

    All methods require and return objects of type
    :class:`hgvs.sequencevariant.SequenceVariant`.

    """

    def __init__(
        self,
        hdp: hgvs.dataproviders.interface.Interface,
        replace_reference=hgvs.global_config.mapping.replace_reference,
        prevalidation_level=hgvs.global_config.mapping.prevalidation_level,
        add_gene_symbol=hgvs.global_config.mapping.add_gene_symbol,
    ):
        """
        :param bool replace_reference: replace reference (entails additional network access)
        :param str prevalidation_level: None or Intrinsic or Extrinsic validation before mapping

        """
        self.hdp = hdp
        self.replace_reference = replace_reference
        self.add_gene_symbol = add_gene_symbol
        if prevalidation_level is None:
            self.prevalidation_level = PrevalidationLevel.NONE
        else:
            self.prevalidation_level = PrevalidationLevel[prevalidation_level.upper()]
        if self.prevalidation_level == PrevalidationLevel.NONE:
            self._validator = None
        elif self.prevalidation_level == PrevalidationLevel.INTRINSIC:
            self._validator = hgvs.validator.IntrinsicValidator(strict=False)
        else:
            self._validator = hgvs.validator.Validator(self.hdp, strict=False)
        self.left_normalizer = hgvs.normalizer.Normalizer(
            hdp, shuffle_direction=5, variantmapper=self
        )

    # ############################################################################
    # g⟷t

    def g_to_t(self, var_g, tx_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        if var_g.type not in "gm":
            raise HGVSInvalidVariantError("Expected a g. or m. variant; got " + str(var_g))

        if self._validator:
            self._validator.validate(var_g)
        var_g.fill_ref(self.hdp)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method=alt_aln_method
        )
        if mapper.is_coding_transcript:
            var_out = VariantMapper.g_to_c(
                self, var_g=var_g, tx_ac=tx_ac, alt_aln_method=alt_aln_method
            )
        else:
            var_out = VariantMapper.g_to_n(
                self, var_g=var_g, tx_ac=tx_ac, alt_aln_method=alt_aln_method
            )
        return var_out

    def t_to_g(self, var_t, alt_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        if var_t.type not in "cn":
            raise HGVSInvalidVariantError("Expected a c. or n. variant; got " + str(var_t))
        if self._validator:
            self._validator.validate(var_t)
        var_t.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        if var_t.type == "c":
            var_out = VariantMapper.c_to_g(
                self, var_c=var_t, alt_ac=alt_ac, alt_aln_method=alt_aln_method
            )
        else:
            var_out = VariantMapper.n_to_g(
                self, var_n=var_t, alt_ac=alt_ac, alt_aln_method=alt_aln_method
            )
        return var_out

    # ############################################################################
    # g⟷n
    def g_to_n(self, var_g, tx_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        """Given a parsed g. variant, return a n. variant on the specified
        transcript using the specified alignment method (default is
        "splign" from NCBI).

        :param hgvs.sequencevariant.SequenceVariant var_g: a variant object
        :param str tx_ac: a transcript accession (e.g., NM_012345.6 or ENST012345678)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`) using transcript (n.) coordinates
        :raises HGVSInvalidVariantError: if var_g is not of type "g"

        """

        if var_g.type not in "gm":
            raise HGVSInvalidVariantError("Expected a g. or m. variant; got " + str(var_g))
        if self._validator:
            self._validator.validate(var_g)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method=alt_aln_method
        )

        if (
            mapper.strand == -1
            and not hgvs.global_config.mapping.strict_bounds
            and not mapper.g_interval_is_inbounds(var_g.posedit.pos)
        ):
            _logger.info("Renormalizing out-of-bounds minus strand variant on genomic sequence")
            var_g = self.left_normalizer.normalize(var_g)

        var_g.fill_ref(self.hdp)
        pos_n = mapper.g_to_n(var_g.posedit.pos)
        if not pos_n.uncertain:
            edit_n = self._convert_edit_check_strand(mapper.strand, var_g.posedit.edit)
            if (
                edit_n.type == "ins"
                and pos_n.start.offset == 0
                and pos_n.end.offset == 0
                and pos_n.end - pos_n.start > 1
            ):
                pos_n.start.base += 1
                pos_n.end.base -= 1
                edit_n.ref = ""
            elif self._variant_has_internal_gap(mapper, var_g):
                edit_n, pos_n = self._get_altered_tx_sequence(
                    mapper.strand, mapper, var_g.posedit.pos, var_g, pos_n
                )
            elif self._variant_spans_i_segment(mapper, var_g):
                expanded_pos_g = self._expand_pos_g_for_adjacent_gap(mapper, var_g)
                edit_n, pos_n = self._get_altered_tx_sequence(
                    mapper.strand, mapper, expanded_pos_g, var_g, pos_n
                )
        else:
            # variant at alignment gap
            pos_g = mapper.n_to_g(pos_n)
            if self._variant_spans_i_segment(mapper, var_g):
                edit_n, pos_n = self._get_altered_tx_sequence(mapper.strand, mapper, pos_g, var_g, pos_n)
            else:
                edit_n = hgvs.edit.NARefAlt(
                    ref="", alt=self._get_altered_sequence(mapper.strand, pos_g, var_g)
                )
        pos_n.uncertain = var_g.posedit.pos.uncertain
        var_n = hgvs.sequencevariant.SequenceVariant(
            ac=tx_ac, type="n", posedit=hgvs.posedit.PosEdit(pos_n, edit_n)
        )

        s, e = get_start_end(var_n)

        if (
            self.replace_reference
            and s is not None
            and s.base >= 0
            and e is not None
            and e.base < mapper.tgt_len
        ):
            self._replace_reference(var_n)
        if self.add_gene_symbol:
            self._update_gene_symbol(var_n, var_g.gene)
        return var_n

    def n_to_g(self, var_n, alt_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        """Given a parsed n. variant, return a g. variant on the specified
        transcript using the specified alignment method (default is
        "splign" from NCBI).

        :param hgvs.sequencevariant.SequenceVariant var_n: a variant object
        :param str alt_ac: a reference sequence accession (e.g., NC_000001.11)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`)
        :raises HGVSInvalidVariantError: if var_n is not of type "n"

        """

        if not (var_n.type == "n"):
            raise HGVSInvalidVariantError("Expected a n. variant; got " + str(var_n))
        if self._validator:
            self._validator.validate(var_n)
        var_n.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=var_n.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method
        )
        pos_g = mapper.n_to_g(var_n.posedit.pos)
        if not pos_g.uncertain:
            edit_g = self._convert_edit_check_strand(mapper.strand, var_n.posedit.edit)
            if edit_g.type == "ins" and pos_g.end - pos_g.start > 1:
                pos_g.start.base += 1
                pos_g.end.base -= 1
                edit_g.ref = ""
        else:
            # variant at alignment gap
            pos_n = mapper.g_to_n(pos_g)
            edit_g = hgvs.edit.NARefAlt(
                ref="", alt=self._get_altered_sequence(mapper.strand, pos_n, var_n)
            )
        pos_g.uncertain = var_n.posedit.pos.uncertain
        var_g = hgvs.sequencevariant.SequenceVariant(
            ac=alt_ac, type="g", posedit=hgvs.posedit.PosEdit(pos_g, edit_g)
        )
        if self.replace_reference:
            self._replace_reference(var_g, alt_ac, alt_aln_method)
        # No gene symbol for g. variants (actually, *should* for NG, but no way to distinguish)
        return var_g

    # ############################################################################
    # g⟷c
    def g_to_c(self, var_g, tx_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        """Given a parsed g. variant, return a c. variant on the specified
        transcript using the specified alignment method (default is
        "splign" from NCBI).

        :param hgvs.sequencevariant.SequenceVariant var_g: a variant object
        :param str tx_ac: a transcript accession (e.g., NM_012345.6 or ENST012345678)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`) using CDS coordinates
        :raises HGVSInvalidVariantError: if var_g is not of type "g"

        """

        if var_g.type not in "gm":
            raise HGVSInvalidVariantError("Expected a g. or m. variant; got " + str(var_g))
        if self._validator:
            self._validator.validate(var_g)

        var_g.fill_ref(self.hdp)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method=alt_aln_method
        )
        pos_c = mapper.g_to_c(var_g.posedit.pos)
        if not pos_c.uncertain:
            edit_c = self._convert_edit_check_strand(mapper.strand, var_g.posedit.edit)
            if (
                edit_c.type == "ins"
                and pos_c.start.offset == 0
                and pos_c.end.offset == 0
                and pos_c.start.datum == pos_c.end.datum
                and pos_c.end - pos_c.start > 1
            ):
                pos_c.start.base += 1
                pos_c.end.base -= 1
                edit_c.ref = ""
            elif self._variant_has_internal_gap(mapper, var_g):
                # Gap segment lies inside the variant interval; both endpoints are in = regions
                # so pos_c.uncertain is False, but the transcript length differs from genomic length.
                # Recalculate the transcript edit from sequences to absorb the gap.
                edit_c, pos_c = self._get_altered_tx_sequence(
                    mapper.strand, mapper, var_g.posedit.pos, var_g, pos_c
                )
            elif self._variant_spans_i_segment(mapper, var_g):
                # I-segment is immediately adjacent to the variant interval (not internal, not uncertain).
                # The alt allele may cancel with the adjacent I-segment yielding a simpler tx edit.
                expanded_pos_g = self._expand_pos_g_for_adjacent_gap(mapper, var_g)
                edit_c, pos_c = self._get_altered_tx_sequence(
                    mapper.strand, mapper, expanded_pos_g, var_g, pos_c
                )
        else:
            # variant at alignment gap
            pos_g = mapper.c_to_g(pos_c)
            if self._variant_spans_i_segment(mapper, var_g):
                edit_c, pos_c = self._get_altered_tx_sequence(mapper.strand, mapper, pos_g, var_g, pos_c)
            else:
                edit_c = hgvs.edit.NARefAlt(
                    ref="", alt=self._get_altered_sequence(mapper.strand, pos_g, var_g)
                )
        pos_c.uncertain = var_g.posedit.pos.uncertain
        var_c = hgvs.sequencevariant.SequenceVariant(
            ac=tx_ac, type="c", posedit=hgvs.posedit.PosEdit(pos_c, edit_c)
        )
        if self.replace_reference:
            self._replace_reference(var_c)
        if self.add_gene_symbol:
            self._update_gene_symbol(var_c, var_g.gene)
        return var_c

    def c_to_g(self, var_c, alt_ac, alt_aln_method=hgvs.global_config.mapping.alt_aln_method):
        """Given a parsed c. variant, return a g. variant on the specified
        transcript using the specified alignment method (default is
        "splign" from NCBI).

        :param hgvs.sequencevariant.SequenceVariant var_c: a variant object
        :param str alt_ac: a reference sequence accession (e.g., NC_000001.11)
        :param str alt_aln_method: the alignment method; valid values depend on data source
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`)
        :raises HGVSInvalidVariantError: if var_c is not of type "c"

        """

        if not (var_c.type == "c"):
            raise HGVSInvalidVariantError("Expected a cDNA (c.); got " + str(var_c))
        if self._validator:
            self._validator.validate(var_c)
        var_c.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=var_c.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method
        )
        pos_g = mapper.c_to_g(var_c.posedit.pos)

        if not pos_g.uncertain:
            edit_g = self._convert_edit_check_strand(mapper.strand, var_c.posedit.edit)
            if edit_g.type == "ins" and pos_g.end - pos_g.start > 1:
                pos_g.start.base += 1
                pos_g.end.base -= 1
                edit_g.ref = ""
        else:
            # variant at alignment gap
            var_n = copy.deepcopy(var_c)
            var_n.posedit.pos = mapper.c_to_n(var_c.posedit.pos)
            var_n.type = "n"
            pos_n = mapper.g_to_n(pos_g)
            edit_g = hgvs.edit.NARefAlt(
                ref="", alt=self._get_altered_sequence(mapper.strand, pos_n, var_n)
            )
        pos_g.uncertain = var_c.posedit.pos.uncertain
        var_g = hgvs.sequencevariant.SequenceVariant(
            ac=alt_ac, type="g", posedit=hgvs.posedit.PosEdit(pos_g, edit_g)
        )
        if self.replace_reference:
            self._replace_reference(var_g, alt_ac, alt_aln_method)
        return var_g

    # ############################################################################
    # c⟷n
    def c_to_n(
        self,
        var_c,
        alt_ac=None,
        alt_aln_method=hgvs.global_config.mapping.alt_aln_method,
    ):
        """Given a parsed c. variant, return a n. variant on the specified
        transcript using the specified alignment method (default is
        "transcript" indicating a self alignment).

        :param hgvs.sequencevariant.SequenceVariant var_c: a variant object
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`)
        :raises HGVSInvalidVariantError: if var_c is not of type "c"

        """

        if not (var_c.type == "c"):
            raise HGVSInvalidVariantError("Expected a cDNA (c.); got " + str(var_c))
        if self._validator:
            self._validator.validate(var_c)
        var_c.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=var_c.ac, alt_ac=var_c.ac, alt_aln_method="transcript"
        )
        pos_n = mapper.c_to_n(var_c.posedit.pos)
        if isinstance(var_c.posedit.edit, hgvs.edit.NARefAlt | hgvs.edit.Dup | hgvs.edit.Inv):
            edit_n = copy.deepcopy(var_c.posedit.edit)
        else:
            msg = "Only NARefAlt/Dup/Inv types are currently implemented"
            raise HGVSUnsupportedOperationError(msg)
        var_n = hgvs.sequencevariant.SequenceVariant(
            ac=var_c.ac, type="n", posedit=hgvs.posedit.PosEdit(pos_n, edit_n)
        )
        if self.replace_reference:
            self._replace_reference(var_n, alt_ac, alt_aln_method)
        if self.add_gene_symbol:
            self._update_gene_symbol(var_n, var_c.gene)
        return var_n

    def n_to_c(
        self,
        var_n,
        alt_ac=None,
        alt_aln_method=hgvs.global_config.mapping.alt_aln_method,
    ):
        """Given a parsed n. variant, return a c. variant on the specified
        transcript using the specified alignment method (default is
        "transcript" indicating a self alignment).

        :param hgvs.sequencevariant.SequenceVariant var_n: a variant object
        :returns: variant object (:class:`hgvs.sequencevariant.SequenceVariant`)
        :raises HGVSInvalidVariantError: if var_n is not of type "n"

        """

        if not (var_n.type == "n"):
            raise HGVSInvalidVariantError("Expected n. variant; got " + str(var_n))
        if self._validator:
            self._validator.validate(var_n)
        var_n.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        mapper = self._fetch_AlignmentMapper(
            tx_ac=var_n.ac, alt_ac=var_n.ac, alt_aln_method="transcript"
        )
        pos_c = mapper.n_to_c(var_n.posedit.pos)
        if isinstance(var_n.posedit.edit, hgvs.edit.NARefAlt | hgvs.edit.Dup | hgvs.edit.Inv):
            edit_c = copy.deepcopy(var_n.posedit.edit)
        else:
            msg = "Only NARefAlt/Dup/Inv types are currently implemented"
            raise HGVSUnsupportedOperationError(msg)
        var_c = hgvs.sequencevariant.SequenceVariant(
            ac=var_n.ac, type="c", posedit=hgvs.posedit.PosEdit(pos_c, edit_c)
        )
        if self.replace_reference:
            self._replace_reference(var_c, alt_ac, alt_aln_method)
        if self.add_gene_symbol:
            self._update_gene_symbol(var_c, var_n.gene)
        return var_c

    # ############################################################################
    # c ⟶ p
    def c_to_p(
        self,
        var_c,
        pro_ac=None,
        alt_ac=None,
        alt_aln_method=hgvs.global_config.mapping.alt_aln_method,
        translation_table=TranslationTable.standard,
    ):
        """
        Converts a c. SequenceVariant to a p. SequenceVariant on the specified protein accession
        Author: Rudy Rico

        :param SequenceVariant var_c: hgvsc tag
        :param str pro_ac: protein accession
        :rtype: hgvs.sequencevariant.SequenceVariant

        """

        if not (var_c.type == "c"):
            raise HGVSInvalidVariantError("Expected a cDNA (c.) variant; got " + str(var_c))
        if self._validator:
            self._validator.validate(var_c)
        var_c.fill_ref(self.hdp, alt_ac=alt_ac, alt_aln_method=alt_aln_method)
        reference_data = RefTranscriptData(
            self.hdp, var_c.ac, pro_ac, translation_table=translation_table
        )
        # Use the actual translation table from reference_data (may have been auto-detected)
        builder = altseqbuilder.AltSeqBuilder(
            var_c, reference_data, translation_table=reference_data.translation_table
        )

        # TODO: handle case where you get 2+ alt sequences back;
        # currently get list of 1 element loop structure implemented
        # to handle this, but doesn't really do anything currently.
        all_alt_data = builder.build_altseq()

        var_ps = []
        for alt_data in all_alt_data:
            builder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data, alt_data)
            var_p = builder.build_hgvsp()
            var_ps.append(var_p)

        var_p = var_ps[0]

        if self.add_gene_symbol:
            self._update_gene_symbol(var_p, var_c.gene)

        return var_p

    ############################################################################
    # Internal methods

    def _replace_reference(  # noqa: PLR0912
        self, var, alt_ac=None, alt_aln_method=hgvs.global_config.mapping.alt_aln_method
    ):
        """fetch reference sequence for variant and update (in-place) if necessary"""

        if var.type not in "cgmnr":
            msg = "Can only update references for type c, g, m, n, r"
            raise HGVSUnsupportedOperationError(msg)

        if var.posedit.edit.type in ("ins", "con"):
            # these types have no reference sequence (zero-width), so return as-is
            return var

        pos = var.posedit.pos
        mapper = None
        if (isinstance(pos.start, hgvs.location.BaseOffsetPosition) and pos.start.offset != 0) or (
            isinstance(pos.end, hgvs.location.BaseOffsetPosition) and pos.end.offset != 0
        ):
            if var.type == "r":
                _logger.info("Can't update reference sequence for intronic variant %s", var)
                return var

            if alt_ac is None:
                _logger.info(
                    "Can't update reference sequence for intronic variant %s without alt_ac",
                    var,
                )
                return var

            # For c. variants, we need coords on underlying sequences
            if var.type == "c":
                mapper = self._fetch_AlignmentMapper(
                    tx_ac=var.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method
                )
                pos = mapper.c_to_g(var.posedit.pos)
                ac = alt_ac
                _type = "g"
            if var.type == "n":
                mapper = self._fetch_AlignmentMapper(
                    tx_ac=var.ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method
                )
                pos = mapper.n_to_g(var.posedit.pos)
                ac = alt_ac
                _type = "g"
        elif var.type == "c":
            # For c. variants, we need coords on underlying sequences
            mapper = self._fetch_AlignmentMapper(
                tx_ac=var.ac, alt_ac=var.ac, alt_aln_method="transcript"
            )
            pos = mapper.c_to_n(var.posedit.pos)
            ac = var.ac
            _type = var.type
        else:
            pos = var.posedit.pos
            ac = var.ac
            _type = var.type

        seq_start, seq_end = get_start_end_interbase(pos, outer_confidence=True)
        # When strict_bounds is False and an error occurs, return
        # variant as-is

        if seq_start is None or seq_end is None or seq_start < 0:
            # this is an out-of-bounds variant
            return var

        seq = self.hdp.get_seq(ac, seq_start, seq_end)

        if len(seq) != seq_end - seq_start:
            # tried to read beyond seq end; this is an out-of-bounds variant
            return var

        if _type == "g" and mapper and mapper.strand == -1:
            seq = reverse_complement(seq)

        edit = var.posedit.edit
        if edit.ref != seq:
            _logger.debug("Replaced reference sequence in %s with %s", var, seq)
            edit.ref = seq

        return var

    @lru_cache(maxsize=hgvs.global_config.lru_cache.maxsize)
    def _fetch_AlignmentMapper(self, tx_ac, alt_ac, alt_aln_method):
        """
        Get a new AlignmentMapper for the given transcript accession (ac),
        possibly caching the result.
        """
        return hgvs.alignmentmapper.AlignmentMapper(
            self.hdp, tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method
        )

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
                    ref=ref,
                    alt=reverse_complement(edit_in.alt),
                )
        elif isinstance(edit_in, hgvs.edit.Dup):
            if strand == 1:
                edit_out = copy.deepcopy(edit_in)
            else:
                edit_out = hgvs.edit.Dup(ref=reverse_complement(edit_in.ref))
        elif isinstance(edit_in, hgvs.edit.Inv):
            if strand == 1:
                edit_out = copy.deepcopy(edit_in)
            else:
                try:
                    int(edit_in.ref)
                    ref = edit_in.ref
                except (ValueError, TypeError):
                    ref = reverse_complement(edit_in.ref)
                edit_out = hgvs.edit.Inv(ref=ref)
        else:
            msg = "Only NARefAlt/Dup/Inv types are currently implemented"
            raise NotImplementedError(msg)
        return edit_out

    def _get_altered_sequence(self, strand, interval, var):
        seq = list(self.hdp.get_seq(var.ac, interval.start.base - 1, interval.end.base))
        # positions are 0-based and half-open
        pos_start = var.posedit.pos.start.base - interval.start.base
        pos_end = var.posedit.pos.end.base - interval.start.base + 1
        edit = var.posedit.edit

        if edit.type == "sub":
            seq[pos_start] = edit.alt
        elif edit.type == "del":
            del seq[pos_start:pos_end]
        elif edit.type == "ins":
            seq.insert(pos_start + 1, edit.alt)
        elif edit.type == "delins":
            del seq[pos_start:pos_end]
            seq.insert(pos_start, edit.alt)
        elif edit.type == "dup":
            seq.insert(pos_end, "".join(seq[pos_start:pos_end]))
        elif edit.type == "inv":
            seq[pos_start:pos_end] = list(reverse_complement("".join(seq[pos_start:pos_end])))
        elif edit.type == "identity":
            pass
        else:
            msg = f"Getting altered sequence for {edit.type} is unsupported"
            raise HGVSUnsupportedOperationError(msg)

        seq = "".join(seq)
        if strand == -1:
            seq = reverse_complement(seq)
        return seq

    def _variant_spans_i_segment(self, mapper, var_g):
        """Return True if var_g straddles a normal-to-I-segment boundary in the CIGAR.

        I-segments are genomic positions with no transcript equivalent. The fix is needed
        only when the variant spans *across* an I-segment boundary (one endpoint in a matching
        "=" region, the other landing in an I-segment). When the entire variant lies within an
        I-segment, the existing uncertain-path logic handles it correctly.
        """
        gc_offset = mapper.gc_offset
        ref_pos = mapper.cigarmapper.ref_pos
        cigar_op = mapper.cigarmapper.cigar_op
        # Guard: partially uncertain variants may have Interval positions without .base
        try:
            start_offset = var_g.posedit.pos.start.base - 1 - gc_offset
            end_offset = var_g.posedit.pos.end.base - gc_offset
        except AttributeError:
            return False

        # Determine the CIGAR op at each endpoint
        start_op = end_op = None
        for i, op in enumerate(cigar_op):
            if start_op is None and start_offset < ref_pos[i + 1]:
                start_op = op
            if end_op is None and end_offset < ref_pos[i + 1]:
                end_op = op
            if start_op is not None and end_op is not None:
                break

        # If either endpoint falls outside the alignment, we can't determine the op
        if start_op is None or end_op is None:
            return False

        # Apply fix only when at least one endpoint is in "I" and at least one is not
        ops = {start_op, end_op}
        return "I" in ops and ops != {"I"}

    def _expand_pos_g_for_adjacent_gap(self, mapper, var_g):
        """Return a copy of var_g.posedit.pos expanded to include any adjacent I-segment.

        When the variant's interbase end sits at the start of an I-segment (or its
        interbase start sits at the end of an I-segment), extend pos_g so that
        _get_altered_tx_sequence sees the full I-segment in seq and i_offsets.
        """
        gc_offset = mapper.gc_offset
        ref_pos = mapper.cigarmapper.ref_pos
        cigar_op = mapper.cigarmapper.cigar_op

        pos_g = copy.deepcopy(var_g.posedit.pos)

        # Check right side: interbase end = var_g.posedit.pos.end.base - gc_offset
        end_offset = var_g.posedit.pos.end.base - gc_offset
        for i, op in enumerate(cigar_op):
            if op == "I" and ref_pos[i] == end_offset:
                # I-segment starts exactly at the interbase end of the variant
                pos_g.end.base = ref_pos[i + 1] + gc_offset
                break

        # Check left side: interbase start = var_g.posedit.pos.start.base - 1 - gc_offset
        start_offset = var_g.posedit.pos.start.base - 1 - gc_offset
        for i, op in enumerate(cigar_op):
            if op == "I" and ref_pos[i + 1] == start_offset:
                # I-segment ends exactly at the interbase start of the variant
                pos_g.start.base = ref_pos[i] + gc_offset + 1
                break

        return pos_g

    def _gap_segments_within_pos_g(self, mapper, pos_g):
        """Return gap segments (I and D CIGAR ops) that fall within genomic interval pos_g.

        I-segments: genomic bases with no transcript equivalent.
        D-segments: transcript-only positions at an interbase genomic location.

        Returns:
            dict with keys:
              "I": set of 0-based offsets (relative to pos_g.start.base) for I-segment bases
              "D": list of (interbase_offset, tgt_start, tgt_end) for each D-segment
                   where interbase_offset is 0-based relative to pos_g.start.base,
                   tgt_start/tgt_end are 0-based transcript positions of the D-segment bases
        """
        i_offsets = set()
        d_segments = []
        cigar_ops = mapper.cigarmapper.cigar_op
        ref_pos = mapper.cigarmapper.ref_pos
        tgt_pos = mapper.cigarmapper.tgt_pos
        gc_offset = mapper.gc_offset

        for i, op in enumerate(cigar_ops):
            if op == "I":
                # I advances genome/ref only — absolute genomic positions (1-based, inclusive)
                abs_start = ref_pos[i] + gc_offset + 1
                abs_end = ref_pos[i + 1] + gc_offset  # inclusive
                overlap_start = max(abs_start, pos_g.start.base)
                overlap_end = min(abs_end, pos_g.end.base)
                for p in range(overlap_start, overlap_end + 1):
                    i_offsets.add(p - pos_g.start.base)
            elif op == "D":
                # D advances transcript/tgt only — genomic interbase position is ref_pos[i] (same as ref_pos[i+1])
                # The D-segment sits at interbase position gc_offset + ref_pos[i]
                # It is "inside" pos_g if its genomic interbase position falls strictly between start and end
                genomic_interbase = ref_pos[i] + gc_offset  # 0-based interbase genomic position
                # pos_g uses 1-based inclusive; interbase is between pos_g.start.base-1 and pos_g.start.base
                # A D-segment is inside pos_g if: pos_g.start.base <= genomic_interbase+1 <= pos_g.end.base
                # i.e., it sits between two bases both within pos_g
                if pos_g.start.base <= genomic_interbase < pos_g.end.base:
                    interbase_offset = genomic_interbase - (pos_g.start.base - 1)  # offset relative to pos_g start
                    d_segments.append((interbase_offset, tgt_pos[i], tgt_pos[i + 1]))

        return {"I": i_offsets, "D": d_segments}

    def _variant_has_internal_gap(self, mapper, var_g):
        """Return True if any I or D CIGAR segment falls strictly inside var_g's interval.

        This is distinct from _variant_spans_i_segment which checks endpoint CIGAR ops.
        Internal gaps occur when both endpoints are in normal (=) segments but a gap
        segment lies between them.
        """
        gaps = self._gap_segments_within_pos_g(mapper, var_g.posedit.pos)
        return bool(gaps["I"] or gaps["D"])

    def _get_altered_tx_sequence(self, strand, mapper, pos_g, var_g, tx_pos):
        """Compute transcript-level edit for variants at alignment gaps (uncertain path).

        Correctly handles I-segment positions (genomic bases with no transcript equivalent)
        by filtering them out when building the transcript reference and alt sequences.

        Returns (NARefAlt, adjusted_tx_pos).
        """
        # Fetch genomic reference for pos_g
        seq = list(self.hdp.get_seq(var_g.ac, pos_g.start.base - 1, pos_g.end.base))

        # Find gap segment positions (0-based offsets into seq)
        gaps = self._gap_segments_within_pos_g(mapper, pos_g)
        i_offsets = gaps["I"]
        # NOTE: gaps["D"] (transcript-only bases with no genomic counterpart) are detected but
        # not yet used for sequence reconstruction. Variants spanning D-segments are routed here
        # via _variant_has_internal_gap, and the deletion/delins mapping happens to be correct
        # because D-segment bases are not present in the fetched genomic seq. Full D-segment
        # splicing (inserting transcript-only bases into tx_ref_str) is a future enhancement.

        # Variant boundaries within seq (0-based)
        var_start = var_g.posedit.pos.start.base - pos_g.start.base
        var_end = var_g.posedit.pos.end.base - pos_g.start.base + 1
        edit = var_g.posedit.edit

        # Build transcript ref = all genomic bases excluding I-segment positions
        tx_ref_str = "".join(seq[j] for j in range(len(seq)) if j not in i_offsets)

        # Build transcript alt by applying variant to prefix/suffix (filtered of I-bases)
        if edit.type == "ins":
            # Insertion between two positions: prefix includes var_start base
            prefix_tx = [seq[j] for j in range(var_start + 1) if j not in i_offsets]
            suffix_tx = [seq[j] for j in range(var_start + 1, len(seq)) if j not in i_offsets]
            tx_alt_str = "".join(prefix_tx) + (edit.alt or "") + "".join(suffix_tx)
        elif edit.type == "sub":
            prefix_tx = [seq[j] for j in range(var_start) if j not in i_offsets]
            suffix_tx = [seq[j] for j in range(var_start + 1, len(seq)) if j not in i_offsets]
            tx_alt_str = "".join(prefix_tx) + (edit.alt or "") + "".join(suffix_tx)
        elif edit.type == "del":
            prefix_tx = [seq[j] for j in range(var_start) if j not in i_offsets]
            suffix_tx = [seq[j] for j in range(var_end, len(seq)) if j not in i_offsets]
            tx_alt_str = "".join(prefix_tx) + "".join(suffix_tx)
        elif edit.type in ("delins", "dup", "inv", "identity"):
            prefix_tx = [seq[j] for j in range(var_start) if j not in i_offsets]
            if edit.type == "delins":
                # For delins, include the adjacent I-seg base at exactly var_end in the suffix.
                # This base is in the genomic reference but absent from the transcript; including
                # it allows the prefix/suffix trimming step to cancel the double gap and produce
                # the correct minimal transcript edit (e.g. 6-base delins → single SNV).
                suffix_tx = [seq[j] for j in range(var_end, len(seq)) if j not in i_offsets or j == var_end]
                ins_seq = edit.alt or ""
            else:
                # For dup/inv/identity the duplicated/inverted/copied sequence must be
                # transcript-only bases — exclude any I-segment positions from the variant range.
                suffix_tx = [seq[j] for j in range(var_end, len(seq)) if j not in i_offsets]
                clean_variant_seq = "".join(seq[j] for j in range(var_start, var_end) if j not in i_offsets)
                if edit.type == "dup":
                    ins_seq = clean_variant_seq * 2
                elif edit.type == "inv":
                    ins_seq = reverse_complement(clean_variant_seq)
                else:  # identity
                    ins_seq = clean_variant_seq
            tx_alt_str = "".join(prefix_tx) + ins_seq + "".join(suffix_tx)
        else:
            msg = f"_get_altered_tx_sequence: unsupported edit type {edit.type!r}"
            raise HGVSUnsupportedOperationError(msg)

        # Apply strand: convert from genomic order to transcript order
        if strand == -1:
            tx_ref_str = reverse_complement(tx_ref_str)
            tx_alt_str = reverse_complement(tx_alt_str)

        # Trim common prefix and suffix to find the minimal edit
        n_prefix = 0
        while n_prefix < len(tx_ref_str) and n_prefix < len(tx_alt_str) and tx_ref_str[n_prefix] == tx_alt_str[n_prefix]:
            n_prefix += 1

        n_suffix = 0
        max_suffix = min(len(tx_ref_str) - n_prefix, len(tx_alt_str) - n_prefix)
        while n_suffix < max_suffix and tx_ref_str[-(n_suffix + 1)] == tx_alt_str[-(n_suffix + 1)]:
            n_suffix += 1

        ref_trimmed = tx_ref_str[n_prefix: len(tx_ref_str) - n_suffix if n_suffix else None]
        alt_trimmed = tx_alt_str[n_prefix: len(tx_alt_str) - n_suffix if n_suffix else None]

        adjusted_pos = copy.deepcopy(tx_pos)
        adjusted_pos.start.base += n_prefix
        adjusted_pos.end.base -= n_suffix

        return hgvs.edit.NARefAlt(ref=ref_trimmed, alt=alt_trimmed), adjusted_pos

    def _update_gene_symbol(self, var, symbol):
        if not symbol:
            symbol = self.hdp.get_tx_identity_info(var.ac).get("hgnc", None)
        var.gene = symbol
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
