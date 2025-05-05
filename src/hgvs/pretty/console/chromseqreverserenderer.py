from bioutils.sequences import reverse_complement

from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ChromReverseSeqRendered(BasicRenderer):
    """
    ChromReverseSeqRendered is a renderer class that extends BasicRenderer to provide
    functionality for rendering reverse complement sequences with optional color coding.
    Methods
    -------
    legend() -> str:
        Returns a string representing the legend for the sequence orientation.
    display(data: VariantData) -> str:
        Colors the reference sequences with adenine (A, green), thymine (T, red), cytosine (C, yellow),
        and guanine (G, blue) based on the provided VariantData. The sequences are displayed in reverse
        complement form.
    """

    def legend(self) -> str:
        if self.orientation < 0 and self.config.reverse_display:
            arrow = "->"
        else:
            arrow = "<-"

        return f"seq    {arrow} : "

    def display(self, data: VariantData) -> str:
        """colors the ref sequences with adenine (A, green), thymine (T, red), cytosine (C, yellow), and guanine (G, blue)"""
        from hgvs.pretty.console.constants import ENDC, COLOR_MAP

        var_seq = ""
        for p in data.position_details:
            c = reverse_complement(p.ref)

            if not c:
                var_seq += "."
                continue

            if self.config.use_color:
                var_seq += COLOR_MAP.get(c, COLOR_MAP["N"])

            var_seq += c

            if self.config.use_color:
                var_seq += ENDC

        return var_seq
