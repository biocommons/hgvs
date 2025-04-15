import math

from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class TxAligRenderer(BasicRenderer):
    """
    TxAligRenderer is a class that extends BasicRenderer to provide functionality for rendering
    the transcript sequence, optional with alternating colors per codon.
    It includes methods to generate a legend and display data with the transcript sequence.

    Methods
    -------
    legend() -> str:
        Generates a legend string indicating the orientation of the transcript sequence.
    display(data: VariantData) -> str:
        Generates a string representation of the variant data, showing details of the transcript
        for the specified region. It uses color coding if enabled in the configuration.
    """

    def legend(self) -> str:
        if self.orientation > 0:
            orientation = "->"
        elif self.orientation < 0 and self.config.reverse_display:
            orientation = "->"
        else:
            orientation = "<-"
        return f"tx seq {orientation} : "

    def display(self, data: VariantData) -> str:
        """If transcript info is available show the details of the tx for the region."""
        if not data.var_c_or_n:
            return ""

        from hgvs.pretty.console.constants import ENDC, COLOR_MAP

        var_str = ""

        first_c = None
        last_c = None

        coding = False
        c_pos = None

        counter = -1
        for pdata in data.position_details:
            counter += 1

            if pdata.mapped_pos is None:
                var_str += " "
                continue

            n_pos = pdata.n_pos
            c_pos = pdata.c_pos
            cig = pdata.cigar_ref
            c_offset = pdata.c_offset

            if cig == "N" or (c_offset is not None and c_offset != 0):
                coding = False
                var_str += " "
                continue

            if c_pos:
                # if we get here we are coding...
                last_c = f"c.{c_pos}"
            else:
                last_c = f"n.{n_pos}"
                c_pos = n_pos
            if not coding:
                if not first_c:
                    first_c = last_c

                coding = True

            if cig == "=":
                if c_pos:
                    bg_col = math.ceil((c_pos) / 3) % 2
                elif n_pos:
                    bg_col = math.ceil((n_pos) / 3) % 2
                else:
                    var_str += " "
                    continue

                # print(p, c_pos, c3, bg_col )
                if n_pos >= 0:
                    base = pdata.tx
                    if self.config.use_color and c_pos > 0:
                        if bg_col:
                            var_str += COLOR_MAP["codon1"] + base + ENDC
                        else:
                            var_str += COLOR_MAP["codon2"] + base + ENDC
                    else:
                        if c_pos is None or c_pos < 0:
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
