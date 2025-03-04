from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class TxMappingRenderer(BasicRenderer):
    """
    TxMappingRenderer is a class that extends BasicRenderer to provide
    functionality for rendering transcript position information by printing a | eveery 10 bases and a . every 5 bases.

    Methods
    -------
    legend():
        Returns a string representing the legend for the transcript position.
    display(data: VariantData) -> str:
        Generates a string representation of the transcript sequence positions
        based on the provided VariantData. The output string uses specific
        characters to denote different positions and offsets within the
        transcript sequence.
    """

    def legend(self):
        return "tx pos    : "

    def display(self, data: VariantData) -> str:
        """show the position of the transcript seq"""

        var_str = ""

        count = -1
        prev_c_pos = ""
        for pdata in data.position_details:
            count += 1
            if not pdata.mapped_pos:
                var_str += " "
                prev_c_pos = ""
                continue

            c_pos = pdata.c_pos
            if c_pos is None and pdata.n_pos:
                c_pos = pdata.n_pos

            if c_pos is None:
                var_str += " "
                prev_c_pos = c_pos
                continue

            if pdata.cigar_ref == "N":
                var_str += " "
                prev_c_pos = c_pos
                continue

            if len(var_str) > count:
                prev_c_pos = c_pos
                continue

            if (c_pos) % 10 == 0:
                var_str += "|"
                prev_c_pos = c_pos
                continue

            elif (c_pos) % 5 == 0:
                var_str += "."
                prev_c_pos = c_pos
                continue

            elif (
                prev_c_pos
                and prev_c_pos == c_pos
                and pdata.c_offset > 0
                and (pdata.c_offset % 10) == 0
            ):
                var_str += "^"
                prev_c_pos = c_pos
                continue
            elif (
                prev_c_pos
                and prev_c_pos == c_pos
                and pdata.c_offset > 0
                and (pdata.c_offset % 5) == 0
            ):
                var_str += "."
                prev_c_pos = c_pos
                continue

            elif c_pos == 1:
                var_str += "|"

            var_str += " "

        return var_str
