from hgvs.pretty.models import VariantData
from hgvs.pretty.renderer.renderer import BasicRenderer


class TxRefDisagreeRenderer(BasicRenderer):
    """Display tx-ref-disagree positions"""

    def legend(self) -> str:

        return f"tx ref dif: "

    def display(self, data: VariantData) -> str:
        """show differences between tx and ref genome, if there are any"""

        if not data.var_c_or_n:
            return ""

        from hgvs.pretty_print import ENDC, TRED

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

        return var_str
