from hgvs.pretty.models import VariantData
from hgvs.pretty.console.renderer import BasicRenderer


class ChromSeqRendered(BasicRenderer):

    def legend(self) -> str:
        return "seq    -> : "

    def display(self, data: VariantData) -> str:
        """colors the ref sequences with adenine (A, green), thymine (T, red), cytosine (C, yellow), and guanine (G, blue)"""
        from hgvs.pretty.console.constants import ENDC, TBLUE, TGREEN, TRED, TYELLOW

        var_seq = ""
        for p in data.position_details:
            c = p.ref

            if not c:
                var_seq += "."
                continue

            if self.config.useColor:
                if c == "A":
                    var_seq += TGREEN
                elif c == "T":
                    var_seq += TRED
                elif c == "C":
                    var_seq += TYELLOW
                elif c == "G":
                    var_seq += TBLUE

            var_seq += c

            if self.config.useColor:
                var_seq += ENDC

        return var_seq
