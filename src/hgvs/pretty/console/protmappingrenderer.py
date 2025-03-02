from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ProtMappingRenderer(BasicRenderer):
    """
    A renderer class for displaying protein position by printing a | every 10 and a . every 5 positions.

    Methods
    -------
    legend() -> str
        Returns a string representing the legend for the protein mapping display.

    display(data: VariantData) -> str
        Generates a string representation of the protein mapping positions based on the provided VariantData.
    """

    def legend(self):
        return "aa pos    : "

    def display(self, data: VariantData) -> str:
        """show the position of the transcript seq"""

        var_str = ""

        count = -1
        for pdata in data.position_details:
            count += 1
            if not pdata.mapped_pos:
                var_str += " "
                continue

            prot_data = pdata.protein_data

            if not prot_data or prot_data.aa_pos < 0:
                var_str += " "
                continue

            if len(var_str) > count:
                continue

            aa_pos = prot_data.aa_pos
            if (aa_pos + 1) % 10 == 0:
                var_str += "|"
                continue

            elif (aa_pos + 1) % 5 == 0:
                var_str += "."
                continue

            var_str += " "

        return var_str
