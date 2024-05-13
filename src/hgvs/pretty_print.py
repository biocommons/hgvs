import math
from typing import Tuple

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from hgvs.enums import Datum
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig, VariantCoords, VariantData
from hgvs.sequencevariant import SequenceVariant
from hgvs.utils.reftranscriptdata import RefTranscriptData

TGREEN = "\033[32m"  # Green Text
TGREENBG = "\033[30;42m"
TRED = "\033[31m"  # Red Text
TREDBG = "\033[30;41m"
TBLUE = "\033[34m"  # Blue Text
TBLUEBG = "\033[30;44m"
TPURPLE = "\033[35m"  # Purple Text
TPURPLEBG = "\033[30;45m"
TYELLOW = "\033[33m"  # Yellow Text
TYELLOWBG = "\033[30;43m"

ENDC = "\033[m"  # reset to the defaults


class PrettyPrint:
    """A class that provides a pretty display of the genomic context of a variant."""

    def __init__(
        self,
        hdp: hgvs.dataproviders.interface.Interface,
        default_assembly: str = "GRCh37",
        padding_left: int = 20,
        padding_right: int = 20,
        useColor=False,
        showLegend=True,
        infer_hgvs_c=True,
    ):
        """
        :param hdp: HGVS Data Provider Interface-compliant instance
           (see :class:`hgvs.dataproviders.interface.Interface`)

        :param padding: spacing left and right of the variant for display purposes.
        """
        am37: AssemblyMapper = AssemblyMapper(hdp, assembly_name="GRCh37")
        am38: AssemblyMapper = AssemblyMapper(hdp, assembly_name="GRCh38")

        self.config = PrettyConfig(
            hdp,
            am37,
            am38,
            padding_left,
            padding_right,
            default_assembly,
            useColor,
            showLegend,
            infer_hgvs_c,
        )

    # def _display_position_info(self, data:VariantData)->str:

    #     seq_start = data.display_start
    #     seq_end = data.display_end

    #     line_start = f"{seq_start + 1:,} "
    #     line_end = f"{seq_end:,}"
    #     ls = len(line_start)
    #     le = len(line_end)

    #     total_chars = seq_end - seq_start + 1

    #     gap = total_chars - ls - le -1

    #     line = line_start
    #     p = 0
    #     while p < gap:
    #         line += ' '
    #         p+=1

    #     line += line_end
    #     return line

    def _display_position_info(self, data: VariantData) -> str:

        count = -1
        var_seq = ""
        for pdata in data.position_details:
            count += 1
            g = pdata.chromosome_pos

            total_chars = len(var_seq)

            if total_chars > count:
                continue

            if not g:
                var_seq += " "
                continue

            if g % 10 == 0:
                var_seq += f"{g:,} "

            else:
                var_seq += " "

        return var_seq.rstrip()

    def _display_ruler(self, data: VariantData) -> str:
        """Draw a position indicator that diplays | every 10 bases and a . every 5

        seq_start/end in interbase
        """

        ruler = ""

        for pd in data.position_details:
            p = pd.chromosome_pos
            if not p:
                ruler += "_"
                continue

            if p % 10 == 0:
                ruler += "|"
            elif p % 5 == 0:
                ruler += "."
            else:
                ruler += " "
        return ruler

    def _display_seq(self, data: VariantData) -> str:
        """colors the ref sequences with adenine (A, green), thymine (T, red), cytosine (C, yellow), and guanine (G, blue)"""

        var_seq = ""
        for p in data.position_details:
            c = p.ref

            if not c:
                var_seq += "."
                continue

            if self.config.useColor:
                if c == "A":
                    var_seq += TGREEN
                elif c == "T":
                    var_seq += TRED
                elif c == "C":
                    var_seq += TYELLOW
                elif c == "G":
                    var_seq += TBLUE

            var_seq += c

            if self.config.useColor:
                var_seq += ENDC

        return var_seq

    def _display_shuffled_range(
        self, var_g: SequenceVariant, data: VariantData, vc: VariantCoords
    ) -> str:

        seq_start = data.display_start
        seq_end = data.display_end

        split_char = "|"

        start = vc.start
        end = vc.end

        # map interbase coordinates to 1-based display coords:
        start = start + 1

        if var_g.posedit.edit.type == "sub":
            end = start
            if len(vc.alt) == 1:
                split_char = vc.alt

        l = end - start + 1
        if var_g.posedit.edit.type == "ins" and l == 0:
            start = start - 1
            end = end + 1
            split_char = "^"

        if var_g.posedit.edit.type == "del":
            split_char = ""
            if self.config.useColor:
                split_char = TRED
            split_char += "x"
            if self.config.useColor:
                split_char += ENDC

        elif var_g.posedit.edit.type == "identity":
            split_char = "="
        # print(l,start, end , vc)

        if start < seq_start:
            # raise ValueError(f"Can't create shuffled representation, since start {start} < seq_start {seq_start} ")
            return ""
        if end > seq_end:
            return ""
            # raise ValueError(f"Can't create shuffled representation, since end {end} > seq_end {seq_end} ")

        var_str = ""
        in_range = False
        for pdata in data.position_details:
            p = pdata.chromosome_pos

            if not p:
                if in_range:
                    var_str += "-"
                else:
                    var_str += " "
                continue

            if p == start:
                var_str += split_char
                in_range = True
            elif p == end:
                var_str += split_char
                in_range = False
            elif p > end and in_range:
                in_range = False
                var_str += " "
            elif in_range:
                var_str += "-"
            else:
                var_str += " "

        return var_str

    def _infer_hgvs_c(self, var_g: SequenceVariant) -> SequenceVariant:
        if self.config.default_assembly == "GRCh37":
            am = self.config.am37
        else:
            am = self.config.am38

        transcripts = am.relevant_transcripts(var_g)
        if transcripts:
            tx_ac = transcripts[0]
        else:
            return None

        var_c = am.g_to_c(var_g, tx_ac)
        return var_c

    def _display_tx_ref_disagree_str(self, data: VariantData) -> str:
        """show differences between tx and ref genome, if there are any"""

        if not data.var_c:
            return ""

        var_str = ""

        counter = -1
        for p in range(data.display_start + 1, data.display_end + 1):
            counter += 1
            pdata = data.position_details[counter]
            c_offset = pdata.c_offset
            if not pdata.mapped_pos:
                var_str += " "
                continue

            cig = pdata.cigar_ref
            if cig == "=":
                var_str += " "
            elif cig == "N" or c_offset != 0:
                var_str += " "
                continue
            else:
                # an alignment issue, show cigar string
                if self.config.useColor:
                    var_str += TRED + cig + ENDC
                else:
                    var_str += cig

        if var_str.isspace():
            return ""

        if self.config.showLegend:
            return f"tx ref dif: " + var_str

        return var_str

    def _display_transcript_str(self, data: VariantData) -> str:
        """If transcript info is available show the details of the tx for the region."""
        if not data.var_c:
            return ""

        mapper = data.alignmentmapper

        orientation = "->"
        if mapper.strand < 0:
            orientation = "<-"

        var_str = ""
        if self.config.showLegend:
            var_str = f"tx seq {orientation} : "

        first_c = None
        last_c = None

        tx_seq = data.tx_seq

        coding = False
        c_pos = None

        counter = -1
        for pdata in data.position_details:
            p = pdata.chromosome_pos
            counter += 1

            if not pdata.mapped_pos:
                var_str += " "
                continue

            n_pos = pdata.n_pos
            c_pos = pdata.c_pos
            cig = pdata.cigar_ref
            c_offset = pdata.c_offset

            if cig == "N" or c_offset != 0:
                coding = False
                var_str += " "
                continue

            # if we get here we are coding...
            last_c = f"c.{c_pos}"
            if not coding:

                if not first_c:
                    first_c = last_c

                coding = True

            if cig == "=":
                c3 = (c_pos - 1) % 3
                bg_col = math.ceil(c_pos / 3) % 2
                # print(p, c_pos, c3, bg_col )
                if n_pos >= 0:
                    base = pdata.tx
                    if self.config.useColor and c_pos > 0:
                        if bg_col:
                            var_str += TPURPLE + base + ENDC
                        else:
                            var_str += TYELLOW + base + ENDC
                    else:
                        if c_pos < 0:
                            var_str += base.lower()
                        else:
                            var_str += base
                    continue

            elif cig == "X" or cig == "D":
                # for mismatches and tx-insertions show sequence
                var_str += pdata.tx
                continue

            else:
                var_str += "-"

        return var_str

    def _display_protein_str(self, data: VariantData) -> str:
        if not data.var_c:
            return ""

        mapper = data.alignmentmapper

        var_str = ""
        for pdata in data.position_details:

            p = pdata.chromosome_pos

            if not pdata.mapped_pos:
                var_str += " "
                continue

            ref_base = pdata.ref

            c_interval = pdata.c_interval
            if not c_interval:
                var_str += " "
                continue

            cig = pdata.cigar_ref
            c_offset = pdata.c_offset
            if cig == "N" or c_offset != 0:
                var_str += " "
                continue

            if cig == "I":
                var_str += "-"
                continue

            if not ref_base or pdata.tx:
                ref_base = pdata.tx

            protein_data = pdata.protein_data

            if not protein_data:
                var_str += " "
                continue

            aa_char = protein_data.aa_char

            # color init met and stop codon if using color:
            if protein_data.is_init_met:
                if self.config.useColor:
                    aa_char = TGREEN + aa_char + ENDC
                var_str += aa_char
                continue
            if protein_data.is_stop_codon:
                if self.config.useColor:
                    aa_char = TRED + aa_char + ENDC
                var_str += aa_char
                continue

            if not protein_data.var_p.posedit:
                var_str += " "
                continue

            var_str += aa_char
            continue

        if self.config.showLegend:
            legend = "aa seq -> : "
            if mapper.strand < 0:
                legend = "aa seq <- : "
            return legend + var_str

        return var_str

    def _display_c_pos(self, data: VariantData) -> str:
        """show the position of the transcript seq"""
        var_str = ""

        count = -1
        for pdata in data.position_details:
            count += 1
            if not pdata.mapped_pos:
                var_str += " "
                continue

            c_pos = pdata.c_pos

            if c_pos is None:                
                var_str += " "
                continue

            if len(var_str) > count:
                continue

            if (c_pos + 1) % 10 == 0:
                # if pdata.c_interval.start.datum == Datum.CDS_END:
                #     var_str += "*"
                var_str += f"{pdata.c_interval} "
                continue

            elif c_pos == 0:
                var_str += f"{pdata.c_interval} "
                continue
            var_str += " "

        if self.config.showLegend:
            return "          : " + var_str.rstrip()
        return var_str.rstrip()

    def _display_c_pos_ruler(self, data: VariantData) -> str:
        """show the position of the transcript seq"""

        var_str = ""

        count = -1
        for pdata in data.position_details:
            count += 1
            if not pdata.mapped_pos:
                var_str += " "
                continue

            c_pos = pdata.c_pos
            
            if c_pos is None:
                var_str += " "
                continue


            if len(var_str) > count:
                continue

            if (c_pos + 1) % 10 == 0:
                var_str += "|"
                continue

            elif (c_pos + 1) % 5 == 0:
                var_str += "."
                continue
            elif c_pos == 0:
                var_str += "|"

            var_str += " "

        if self.config.showLegend:
            return "tx pos    : " + var_str
        return var_str

    def _map_to_chrom(self, sv: SequenceVariant) -> SequenceVariant:
        """maps a variant to chromosomal coords, if needed."""
        if self.config.default_assembly == "GRCh37":
            am = self.config.am37
        else:
            am = self.config.am38

        if sv.type == "c":
            return am.c_to_g(sv)
        elif sv.type == "n":
            return am.n_to_g(sv)
        elif sv.type == "t":
            return am.t_to_g(sv)

    def _colorize_hgvs(self, hgvs_str: str) -> str:
        if not self.config.useColor:
            return hgvs_str

        spl = hgvs_str.split(":")
        var_str = TPURPLE + spl[0] + ENDC
        var_str += ":"

        sec = spl[1].split(".")
        var_str += TYELLOW + sec[0] + ENDC
        var_str += "."
        var_str += sec[1]

        return var_str

    def display(
        self, sv: SequenceVariant, display_start: int = None, display_end: int = None
    ) -> str:
        """Takes a variant and prints the genomic context around it."""

        var_c = None
        if sv.type == "g":
            var_g = sv
        elif sv.type == "c":
            var_g = self._map_to_chrom(sv)
            var_c = sv
        elif sv.type in ["c", "n"]:
            # map back to genome
            var_g = self._map_to_chrom(sv)

        if not var_c:
            if self.config.infer_hgvs_c:
                var_c = self._infer_hgvs_c(var_g)

        dc = DataCompiler(config=self.config)

        data = dc.data(var_g, var_c, display_start, display_end)

        pos_info = self._display_position_info(data)

        ruler = self._display_ruler(data)

        seq = self._display_seq(data)

        left_shuffled_var = data.left_shuffled
        right_shuffled_var = data.right_shuffled
        fully_justified_var = data.fully_justified

        if self.config.showLegend:
            head = "hgvs      : "
            posh = "          : "
            rule = "chrom pos : "
            seqe = "seq    -> : "
            regi = "region    : "
            refa = "ref>alt   : "
        else:
            head = posh = rule = seqe = regi = refa = ""

        if self.config.useColor:
            var_g_print = self._colorize_hgvs(str(var_g))
        else:
            var_g_print = str(var_g)

        var_str = head + var_g_print + "\n"
        if data.var_c:
            if self.config.useColor:
                var_c_print = self._colorize_hgvs(str(var_c))
            else:
                var_c_print = str(var_c)
            var_str += head + var_c_print + "\n"

        var_str += posh + pos_info + "\n" + rule + ruler + "\n" + seqe + seq + "\n"

        tx_ref_disagree_str = self._display_tx_ref_disagree_str(data)
        if tx_ref_disagree_str:
            var_str += tx_ref_disagree_str + "\n"

        left_shuffled_str = self._display_shuffled_range(sv, data, left_shuffled_var)
        right_shuffled_str = self._display_shuffled_range(sv, data, right_shuffled_var)

        if left_shuffled_str != right_shuffled_str:

            fully_justified_str = self._display_shuffled_range(sv, data, fully_justified_var)
            var_str += regi + fully_justified_str + "\n"

        else:
            var_str += regi + left_shuffled_str + "\n"

        protein_str = self._display_protein_str(data)
        if protein_str:
            var_str += protein_str + "\n"

        transcript_str = self._display_transcript_str(data)
        if transcript_str:
            var_str += transcript_str + "\n"

        c_pos_count_str = self._display_c_pos(data)
        c_pos_str = self._display_c_pos_ruler(data)
        if c_pos_str:
            var_str += c_pos_str + "\n"
        if c_pos_count_str:
            var_str += c_pos_count_str + "\n"

        if left_shuffled_str != right_shuffled_str:
            # TODO: detect repeats?
            fully_justified_ref = fully_justified_var.ref
            fully_justified_alt = fully_justified_var.alt
            var_str += refa + f"{fully_justified_ref}>{fully_justified_alt}"

        return var_str
