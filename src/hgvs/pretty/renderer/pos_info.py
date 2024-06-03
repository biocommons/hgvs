from hgvs.pretty.models import VariantData
from hgvs.pretty.renderer.renderer import BasicRenderer


class ChrPositionInfo(BasicRenderer):

    def legend(self):
        return "          : "

    def display(self, data: VariantData) -> str:

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
