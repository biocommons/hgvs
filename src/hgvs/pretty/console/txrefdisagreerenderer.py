from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class TxRefDisagreeRenderer(BasicRenderer):
    """
    TxRefDisagreeRenderer is a class that extends BasicRenderer to provide functionality for rendering
    differences between transcript (tx) and reference (ref) genome sequences aka "tx-ref-disagree".
    Methods:
        legend() -> str:
            Returns a string representing the legend for the renderer.
        display(data: VariantData) -> str:
            Displays differences between the transcript and reference genome sequences based on the provided
            VariantData. If there are no differences, an empty string is returned.
            Args:
                data (VariantData): The data containing variant information and position details.
            Returns:
                str: A string representation of the differences, with optional color coding if enabled in the
                configuration.
    """

    def legend(self) -> str:
        return "tx ref dif: "

    def display(self, data: VariantData) -> str:
        """show differences between tx and ref genome, if there are any"""

        if not data.var_c_or_n:
            return ""

        from hgvs.pretty.console.constants import ENDC, COLOR_MAP

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
                if self.config.use_color:
                    var_str += COLOR_MAP["tx_ref_disagree"] + cig + ENDC
                else:
                    var_str += cig

        if var_str.isspace():
            return ""

        return var_str
