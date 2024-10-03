from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class TxRefDisagreeRenderer(BasicRenderer):
    """Display tx-ref-disagree positions"""

    def legend(self) -> str:
        return "tx ref dif: "

    def display(self, data: VariantData) -> str:
        """show differences between tx and ref genome, if there are any"""

        if not data.var_c_or_n:
            return ""

        from hgvs.pretty.console.constants import ENDC, TRED

        var_str = ""

        for pdata in data.position_details:            
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
