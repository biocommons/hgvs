from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ChrPositionInfo(BasicRenderer):
    """
    ChrPositionInfo is a class that extends BasicRenderer to provide functionality
    for rendering chromosome position every 10 bases .
    Methods
    -------
    legend() -> str
        Returns a string representing the legend for the chromosome position display.
    display(data: VariantData) -> str
        Takes a VariantData object and returns a formatted string representing
        the chromosome positions. Positions that are multiples of 10 are displayed
        with commas, while others are represented by spaces.
    """

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
