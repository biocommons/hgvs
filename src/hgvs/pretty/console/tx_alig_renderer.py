import math

from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class TxAligRenderer(BasicRenderer):

    def legend(self) -> str:
        orientation = "->"
        if self.orientation < 0:
            orientation = "<-"
        return f"tx seq {orientation} : "

    def display(self, data: VariantData) -> str:
        """If transcript info is available show the details of the tx for the region."""
        if not data.var_c_or_n:
            return ""

        from hgvs.pretty.console.constants import ENDC, TPURPLE, TYELLOW

        var_str = ""

        first_c = None
        last_c = None

        coding = False
        c_pos = None

        counter = -1
        for pdata in data.position_details:

            counter += 1

            if not pdata.mapped_pos:
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
                    bg_col = math.ceil((c_pos + 1) / 3) % 2
                elif n_pos:
                    bg_col = math.ceil((n_pos + 1) / 3) % 2
                else:
                    var_str += " "
                    continue

                # print(p, c_pos, c3, bg_col )
                if n_pos >= 0:
                    base = pdata.tx
                    if self.config.useColor and c_pos > 0:
                        if bg_col:
                            var_str += TPURPLE + base + ENDC
                        else:
                            var_str += TYELLOW + base + ENDC
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
