from hgvs.pretty.models import VariantData
from hgvs.pretty.renderer.renderer import BasicRenderer


class ProtRulerRenderer(BasicRenderer):

    def legend(self) -> str:
        return "          : "

    def display(self, data: VariantData) -> str:
        """show the position of the protein seq"""
        var_str = ""

        count = -1
        prev_aa = -1
        for pdata in data.position_details:
            count += 1

            if not pdata.mapped_pos:
                var_str += " "
                prev_aa = -1
                continue

            prot_data = pdata.protein_data

            if not prot_data or prot_data.aa_pos < 0:
                var_str += " "
                prev_aa = -1
                continue

            if len(var_str) > count:
                continue

            aa_pos = prot_data.aa_pos

            if aa_pos == prev_aa:
                var_str += " "
                continue

            prev_aa = aa_pos

            if (aa_pos + 1) % 10 == 0:
                var_str += f"{aa_pos+1} "
                continue

            elif aa_pos == 0:
                var_str += f"{aa_pos+1} "
                continue

            var_str += " "

        return var_str
