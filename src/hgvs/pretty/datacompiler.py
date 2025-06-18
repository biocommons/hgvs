from typing import List, Tuple

from bioutils.normalize import normalize
from bioutils.sequences import aa1_to_aa3_lut

import hgvs
import hgvs.utils.altseq_to_hgvsp as altseq_to_hgvsp
import hgvs.utils.altseqbuilder as altseqbuilder
from hgvs.exceptions import HGVSInvalidIntervalError
from hgvs.pretty.models import (
    PositionDetail,
    PrettyConfig,
    ProteinData,
    VariantCoords,
    VariantData,
)
from hgvs.sequencevariant import SequenceVariant
from hgvs.utils.reftranscriptdata import RefTranscriptData


class DataCompiler:
    """
    DataCompiler is a class responsible for compiling and processing sequence variant data for pretty printing.
    Attributes:
        config (PrettyConfig): Configuration object containing various settings and data providers.
    Methods:
        __init__(config: PrettyConfig):
            Initializes the DataCompiler with the given configuration.
        get_shuffled_variant(var_g: SequenceVariant, direction: int) -> VariantCoords:
            Takes a sequence variant and returns VariantCoords that have been shuffled accordingly.
        get_position_and_state(sv: hgvs.sequencevariant.SequenceVariant) -> Tuple[int, int, str, str]:
            Get the start, end, ref, alt details of a sequence variant.
        _get_exon_nr(tx_exons, genomic_pos) -> Tuple[int, str]:
            Determines the exon number and feature type (exon or intron) for a given genomic position.
        _get_prot_alt(tx_ac: str, strand: int, reference_data, ref_base, c_interval) -> ProteinData:
            Retrieves the protein alteration data for a given transcript and coding interval.
        data(var_g: SequenceVariant, var_c_or_n: SequenceVariant = None, display_start: int = None, display_end: int = None) -> VariantData:
            Takes a sequence variant and provides all the data needed for pretty printing. This is the main method of the class.
        _backfill_gap_in_ref(var_c_or_n, tx_seq, tx_exons, mapper, reference_data, pdata, cig, prev_c_pos, prev_n_pos):
            Fills in gaps in the reference sequence for regions that have been deleted.
        _populate_with_n_c(var_c_or_n, tx_seq, tx_exons, mapper, reference_data, pdata, cig, n_interval, c_interval):
            Populates the PositionDetail object with nucleotide and coding information.
    """

    def __init__(self, config: PrettyConfig):
        self.config = config

    def get_shuffled_variant(
        self, var_g: SequenceVariant, direction: int
    ) -> VariantCoords:
        """Takes a sequence variant and returns VariantCoords that have been shuffled accordingly."""

        # get shuffled representation:
        if direction == 5:
            shuffle_direction = "LEFTSHUFFLE"
        elif direction == 3:
            shuffle_direction = "RIGHTSHUFFLE"
        else:
            shuffle_direction = "EXPAND"

        chrom_seq = self.config.hdp.get_seq(var_g.ac)

        start, end, ref, alt = self.get_position_and_state(var_g)

        if var_g.posedit.edit.type == "identity":
            return VariantCoords(start, end, ref, alt)

        shuffled_interval, shuffled_alleles = normalize(
            chrom_seq,
            interval=(start, end),
            alleles=(None, alt),
            mode=shuffle_direction,
        )

        return VariantCoords(
            shuffled_interval[0],
            shuffled_interval[1],
            shuffled_alleles[0],
            shuffled_alleles[1],
        )

    def _get_start_end(self, var):
        if isinstance(var.posedit.pos, hgvs.location.BaseOffsetInterval):
            s = var.posedit.pos.start
        elif isinstance(var.posedit.pos, hgvs.location.Interval):
            s = var.posedit.pos.start.start
            if not s.base:
                s = var.posedit.pos.start.end
        else:
            s = var.posedit.pos.start

        if isinstance(var.posedit.pos.end, hgvs.location.BaseOffsetPosition):
            e = var.posedit.pos.end
        elif isinstance(var.posedit.pos.end, hgvs.location.Interval):
            e = var.posedit.pos.end.end
            if not e.base:
                e = var.posedit.pos.end.start
        else:
            e = var.posedit.pos.end
        return s, e

    def get_position_and_state(
        self, sv: hgvs.sequencevariant.SequenceVariant
    ) -> Tuple[int, int, str, str]:
        """
        Get the details of a sequence variant.

        Args:
            sv (hgvs.sequencevariant.SequenceVariant): The sequence variant object.

        Returns:
            tuple: A tuple containing the start position, end position, and state of the variant.

        Raises:
            ValueError: If the HGVS variant type is unsupported.
        """

        if sv.posedit.edit.type == "ins":
            start = sv.posedit.pos.start.base
            end = sv.posedit.pos.start.base
            ref = sv.posedit.edit.ref
            alt = sv.posedit.edit.alt

        elif sv.posedit.edit.type in ("sub", "del", "delins", "identity"):
            if isinstance(sv.posedit.pos.start, hgvs.location.Interval):
                start = sv.posedit.pos.start.end.base - 1
                end = sv.posedit.pos.end.start.base
            else:
                start = sv.posedit.pos.start.base - 1
                end = sv.posedit.pos.end.base
            if sv.posedit.edit.type == "identity":
                alt = self.config.hdp.get_seq(sv.ac, start, end)
                ref = alt
            else:
                alt = sv.posedit.edit.alt or ""
                ref = sv.posedit.edit.ref

        elif sv.posedit.edit.type == "dup":
            start, end = self._get_start_end(sv)
            # start = sv.posedit.pos.start.base - 1
            # end = sv.posedit.pos.end.base
            ref = self.config.hdp.get_seq(sv.ac, start.base - 1, end.base)
            alt = ref + ref

        else:
            raise ValueError(f"HGVS variant type {sv.posedit.edit.type} is unsupported")

        return start, end, ref, alt

    def _get_exon_nr(self, tx_exons, genomic_pos) -> Tuple[int, str]:
        i = -1
        for ex in tx_exons:
            i += 1
            exon_nr = ex["ord"] + 1

            if ex["alt_start_i"] < genomic_pos and ex["alt_end_i"] >= genomic_pos:
                return (exon_nr, "exon")

            if i > 0:
                if ex["alt_strand"] > 0:
                    if (
                        tx_exons[i - 1]["alt_end_i"] < genomic_pos
                        and tx_exons[i]["alt_start_i"] >= genomic_pos
                    ):
                        return (exon_nr, "intron")
                else:
                    if (
                        tx_exons[i]["alt_start_i"] < genomic_pos
                        and tx_exons[i - 1]["alt_end_i"] >= genomic_pos
                    ):
                        return (i, "intron")

        return (-1, "no-overlap")

    def _get_prot_alt(
        self, tx_ac: str, strand: int, reference_data, ref_base, c_interval
    ) -> ProteinData:
        edit = hgvs.edit.NARefAlt(ref=ref_base, alt=ref_base)
        var_c = hgvs.sequencevariant.SequenceVariant(
            ac=tx_ac, type="c", posedit=hgvs.posedit.PosEdit(c_interval, edit)
        )
        builder = altseqbuilder.AltSeqBuilder(var_c, reference_data)
        all_alt_data = builder.build_altseq()

        alt_data = all_alt_data[0]

        if not alt_data.variant_start_aa:
            return None

        aa_index = alt_data.variant_start_aa - 1
        aa = reference_data.aa_sequence[aa_index]

        protaltseqbuilder = altseq_to_hgvsp.AltSeqToHgvsp(reference_data, alt_data)
        var_p = protaltseqbuilder.build_hgvsp()

        is_init_met = False
        if protaltseqbuilder._is_init_met:
            is_init_met = True

        # convert to three letter code
        tlc = aa1_to_aa3_lut[aa]
        is_stop_codon = False
        if tlc == "Ter":
            is_stop_codon = True

        c_pos = int(c_interval.start.base) - 1
        c3 = (c_pos) % 3
        aa_char = tlc[c3]
        if strand < 0 and not self.config.reverse_display:
            aa_char = tlc[2 - c3]

        return ProteinData(
            c_pos, aa, tlc, aa_char, var_p, is_init_met, is_stop_codon, aa_index
        )

    def data(
        self,
        var_g: SequenceVariant,
        var_c_or_n: SequenceVariant = None,
        display_start: int = None,
        display_end: int = None,
    ) -> VariantData:
        """
        Takes a sequence variant and provides all the data needed for pretty printing.
        Args:
            var_g (SequenceVariant): The genomic sequence variant.
            var_c_or_n (SequenceVariant, optional): The coding or non-coding sequence variant. Defaults to None.
            display_start (int, optional): The start position for display. Defaults to None.
            display_end (int, optional): The end position for display. Defaults to None.
        Returns:
            VariantData: An object containing all the necessary data for pretty printing the variant.
        """

        start, end, ref, alt = self.get_position_and_state(var_g)

        ls = self.get_shuffled_variant(var_g, 5)
        rs = self.get_shuffled_variant(var_g, 3)
        fs = self.get_shuffled_variant(var_g, 0)

        if ls.start < start:
            start = ls.start
        if rs.end > end:
            end = rs.end

        if not display_start or display_start > start:
            seq_start = start - self.config.padding_left
        else:
            seq_start = display_start

        if not display_end or display_end < end:
            seq_end = end + self.config.padding_right
        else:
            seq_end = display_end

        if var_c_or_n is not None:
            tx_ac = var_c_or_n.ac
        else:
            tx_ac = ""  # can't show transcript , since there is none.

        alt_ac = var_g.ac
        alt_aln_method = "splign"

        if tx_ac:
            tx_exons = self.config.hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
            tx_exons = sorted(tx_exons, key=lambda e: e["ord"])
        else:
            tx_exons = []

        chrom_seq = self.config.hdp.get_seq(var_g.ac)
        disp_seq = chrom_seq[seq_start:seq_end]

        if tx_ac:
            tx_seq = self.config.hdp.get_seq(tx_ac)

            mapper = self.config.assembly_mapper._fetch_AlignmentMapper(
                tx_ac=tx_ac, alt_ac=var_g.ac, alt_aln_method="splign"
            )

        else:
            tx_seq = ""
            mapper = None

        # we don't know the protein ac, get it looked up:
        pro_ac = None
        if var_c_or_n and var_c_or_n.type == "c":
            var_p = self.config.assembly_mapper.c_to_p(var_c_or_n)
            reference_data = RefTranscriptData(self.config.hdp, tx_ac, pro_ac)
        else:
            var_p = None
            reference_data = None

        position_details: List[PositionDetail] = []
        prev_mapped_pos = None
        prev_c_pos = -1
        prev_n_pos = -1

        for chromosome_pos in range(seq_start + 1, seq_end + 1):
            exon_nr, feat = self._get_exon_nr(tx_exons, chromosome_pos)

            pdata = PositionDetail(
                chromosome_pos=chromosome_pos, exon_nr=exon_nr, variant_feature=feat
            )
            position_details.append(pdata)

            pdata.ref = chrom_seq[chromosome_pos - 1]

            if not mapper:
                continue
            try:
                gr = chromosome_pos - mapper.gc_offset - 1
                pdata.alignment_pos = gr
                mapped_pos, mapped_pos_offset, cig = mapper.cigarmapper.map_ref_to_tgt(
                    gr, '"start"'
                )
            except HGVSInvalidIntervalError:
                continue

            pdata.mapped_pos = mapped_pos
            pdata.mapped_pos_offset = mapped_pos_offset
            pdata.cigar_ref = cig

            if prev_mapped_pos:
                while mapped_pos - prev_mapped_pos > 1:
                    prev_mapped_pos = prev_mapped_pos + 1

                    pdata.mapped_pos = prev_mapped_pos

                    # a region in ref that has been deleted. Fill in gaps.
                    self._backfill_gap_in_ref(
                        var_c_or_n,
                        tx_seq,
                        tx_exons,
                        mapper,
                        reference_data,
                        pdata,
                        cig,
                        prev_c_pos,
                        prev_n_pos,
                    )

                    exon_nr, feat = self._get_exon_nr(tx_exons, chromosome_pos)
                    # prev_mapped_pos += 1
                    pdata = PositionDetail(
                        chromosome_pos=chromosome_pos,
                        exon_nr=exon_nr,
                        variant_feature=feat,
                    )
                    position_details.append(pdata)

                    pdata.alignment_pos = gr
                    pdata.mapped_pos = prev_mapped_pos
                    pdata.mapped_pos_offset = mapped_pos_offset
                    pdata.cigar_ref = cig
                    pdata.ref = chrom_seq[chromosome_pos - 1]

                    if mapper.strand > 0:
                        prev_c_pos += 1
                        prev_n_pos += 1
                    else:
                        prev_c_pos -= 1
                        prev_n_pos -= 1

            prev_mapped_pos = mapped_pos

            g_interval = hgvs.location.Interval(
                start=hgvs.location.SimplePosition(chromosome_pos),
                end=hgvs.location.SimplePosition(chromosome_pos),
            )

            try:
                n_interval = mapper.g_to_n(g_interval)
                if var_c_or_n.type == "c":
                    c_interval = mapper.n_to_c(n_interval)
                else:
                    c_interval = None
            except hgvs.exceptions.HGVSInvalidIntervalError:
                # we are beyond the transcript space, can't set any of the other values.
                continue

            pdata.n_interval = n_interval
            if c_interval is not None:
                pdata.c_interval = c_interval
                c_pos = int(c_interval.start.base)
                prev_c_pos = c_pos
            else:
                prev_c_pos = -1

            n_pos = int(n_interval.start.base)
            prev_n_pos = n_pos

            self._populate_with_n_c(
                var_c_or_n,
                tx_seq,
                tx_exons,
                mapper,
                reference_data,
                pdata,
                cig,
                n_interval,
                c_interval,
            )

        # print(position_details[0].get_header()+"\t"+ProteinData.get_header())
        # for p in position_details:
        #     print(f"{p}\t{p.protein_data}\t")

        strand = mapper.strand if mapper else 1

        if strand < 0 and self.config.reverse_display:
            pd = reversed(position_details)
            position_details = list(pd)

        is_rna = var_c_or_n and var_c_or_n.ac.startswith(
            "NR_"
        )  # not sure how to check this for ENSTs

        vd = VariantData(
            seq_start,
            seq_end,
            ls,
            rs,
            fs,
            disp_seq,
            tx_seq,
            mapper,
            var_g,
            strand,
            var_c_or_n,
            var_p,
            position_details,
            is_rna=is_rna,
        )

        return vd

    def _backfill_gap_in_ref(
        self,
        var_c_or_n,
        tx_seq,
        tx_exons,
        mapper,
        reference_data,
        pdata,
        cig,
        prev_c_pos,
        prev_n_pos,
    ):
        pdata.chromosome_pos = None
        pdata.ref = None

        if mapper.strand > 0:
            pdata.c_pos = prev_c_pos + 1
            pdata.n_pos = prev_n_pos + 1
        else:
            pdata.c_pos = prev_c_pos - 1
            pdata.n_pos = prev_n_pos - 1
        pdata.c_offset = 0
        pdata.cigar_ref = "D"
        pdata.tx = tx_seq[pdata.n_pos]

        n_interval = hgvs.location.Interval(
            start=hgvs.location.BaseOffsetPosition(base=pdata.n_pos, offset=0),
            end=hgvs.location.BaseOffsetPosition(base=pdata.n_pos, offset=0),
        )
        c_interval = mapper.n_to_c(n_interval, strict_bounds=False)
        pdata.c_interval = c_interval

        self._populate_with_n_c(
            var_c_or_n,
            tx_seq,
            tx_exons,
            mapper,
            reference_data,
            pdata,
            cig,
            n_interval,
            c_interval,
        )

    def _populate_with_n_c(
        self,
        var_c_or_n,
        tx_seq,
        tx_exons,
        mapper,
        reference_data,
        pdata,
        cig,
        n_interval,
        c_interval,
    ):
        n_pos = int(n_interval.start.base)

        if c_interval:
            c_pos = int(c_interval.start.base)
            pdata.c_pos = c_pos
            pdata.c_offset = c_interval.start.offset
        else:
            c_pos = None
        # print(f"{n_interval} {c_pos} {c_interval.start.datum}")
        # set transcript_feature:
        tx_ac = var_c_or_n.ac

        pdata.n_pos = n_pos

        pdata.tx = tx_seq[pdata.n_pos - 1]

        coding = True
        if var_c_or_n.type == "n":  # rna coding can't be in protein space
            coding = False
        if cig == "N" or pdata.c_offset != 0:
            coding = False

        if cig == "I":
            coding = False

        ref_base = pdata.ref
        if coding:
            if not ref_base or pdata.tx:
                ref_base = pdata.tx
            prot_data = self._get_prot_alt(
                tx_ac, mapper.strand, reference_data, ref_base, c_interval
            )
            pdata.protein_data = prot_data
