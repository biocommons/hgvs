from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ChrRuler(BasicRenderer):

    def legend(self):
        """returns the legend for this category of display"""
        return "chrom pos : "

    def display(self, data: VariantData) -> str:
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
