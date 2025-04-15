from hgvs.pretty.models import VariantCoords, VariantData
from hgvs.pretty.console.renderer import BasicRenderer
from hgvs.sequencevariant import SequenceVariant


class RegionImpacted(BasicRenderer):
    """
    A class used to represent the region impacted by a variant in a genomic sequence.

    Attributes
    ----------
    var_g : SequenceVariant
        The genomic variant to be displayed.
    vc : VariantCoords
        The coordinates of the variant.
    Methods
    -------
    legend() -> str
        Returns a string representing the legend for the variant display.
    display(data: VariantData) -> str
        Returns a string representing the shuffled variant display based on the provided data.
    """

    def __init__(
        self, config, orientation: int, var_g: SequenceVariant, vc: VariantCoords
    ) -> None:
        super().__init__(config, orientation)

        self.var_g = var_g
        self.vc = vc

    def legend(self) -> str:
        return "region    : "

    def display(self, data: VariantData) -> str:
        from hgvs.pretty.console.constants import ENDC, COLOR_MAP

        seq_start = data.display_start
        seq_end = data.display_end

        split_char = "|"

        start = self.vc.start
        end = self.vc.end

        # map interbase coordinates to 1-based display coords:
        start = start + 1

        if self.var_g.posedit.edit.type == "sub":
            end = start
            if len(self.vc.alt) == 1:
                split_char = self.vc.alt
        vlen = end - start + 1
        if self.var_g.posedit.edit.type == "ins" and vlen == 0:
            start = start - 1
            end = end + 1
            split_char = "^"

        if self.var_g.posedit.edit.type == "del":
            split_char = ""
            if self.config.use_color:
                split_char = COLOR_MAP["del"]
            split_char += "x"
            if self.config.use_color:
                split_char += ENDC

        elif self.var_g.posedit.edit.type == "identity":
            split_char = "="

        if start < seq_start:
            # raise ValueError(f"Can't create shuffled representation, since start {start} < seq_start {seq_start} ")
            return ""
        if end > seq_end:
            return ""
            # raise ValueError(f"Can't create shuffled representation, since end {end} > seq_end {seq_end} ")

        var_str = ""
        in_range = False

        reverse_display = False
        if self.orientation < 0 and self.config.reverse_display:
            reverse_display = True

        for pdata in data.position_details:
            p = pdata.chromosome_pos
            if not p:
                if in_range:
                    var_str += "-"
                else:
                    var_str += " "
                continue

            if not reverse_display and p == start:
                var_str += split_char
                in_range = True
            elif reverse_display and p == end:
                var_str += split_char
                in_range = True
            elif not reverse_display and p == end:
                var_str += split_char
                in_range = False
            elif reverse_display and p == start:
                var_str += split_char
                in_range = False
            elif not reverse_display and p > end and in_range:
                in_range = False
                var_str += " "
            elif reverse_display and p < start and in_range:
                in_range = False
                var_str += " "
            elif in_range:
                var_str += "-"
            else:
                var_str += " "

        return var_str
