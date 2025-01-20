from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ProtSeqRenderer(BasicRenderer):

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

        from hgvs.pretty.console.constants import ENDC, TGREEN, TRED

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

        return var_str
