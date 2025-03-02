import math
from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ProtSeqRenderer(BasicRenderer):
    """
    ProtSeqRenderer is a class that extends BasicRenderer to provide functionality for rendering protein sequences
    with optional color coding for specific amino acids such as the initiation methionine and stop codons.
    Methods
    -------
    legend():
        Returns a string representing the legend for the protein sequence, indicating the direction of the sequence.
    display(data: VariantData) -> str:
        Takes a VariantData object and returns a string representation of the protein sequence with optional color coding.
        The sequence is constructed based on the position details and various attributes of the VariantData object.
    """

    def legend(self):
        if self.orientation > 0:
            arrow = "->"
        elif self.orientation < 0 and self.config.reverse_display:
            arrow = "->"
        else:
            arrow = "<-"

        return f"aa seq {arrow} : "

    def display(self, data: VariantData) -> str:
        if not data.var_c_or_n:
            return ""

        from hgvs.pretty.console.constants import ENDC, COLOR_MAP

        var_str = ""
        for pdata in data.position_details:
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
                if self.config.use_color:
                    aa_char = COLOR_MAP["init_met"] + aa_char + ENDC
                var_str += aa_char
                continue
            if protein_data.is_stop_codon:
                if self.config.use_color:
                    aa_char = COLOR_MAP["stop_codon"] + aa_char + ENDC
                var_str += aa_char
                continue

            if not protein_data.var_p.posedit:
                var_str += " "
                continue

            if self.config.use_color:
                c_pos = pdata.c_pos
                if not c_pos:
                    c_pos = pdata.n_pos
                bg_col = math.ceil((c_pos) / 3) % 2
                if bg_col:
                    var_str += COLOR_MAP["codon1"] + aa_char + ENDC
                else:
                    var_str += COLOR_MAP["codon2"] + aa_char + ENDC
                continue

            var_str += aa_char
            continue

        return var_str
