from hgvs.pretty.models import VariantData
from hgvs.pretty.renderer.renderer import BasicRenderer




class TxRulerRenderer(BasicRenderer):
        
    def legend(self)->str:
        return "          : "

    def display(self, data: VariantData) -> str:
        """show the position of the transcript seq"""
        var_str = ""

        rna_coding = data.var_c_or_n.ac.startswith("NR_")

        count = -1
        for pdata in data.position_details:
            count += 1
            if not pdata.mapped_pos:
                var_str += " "
                continue

            c_pos = pdata.c_pos
            interval = pdata.c_interval
            if rna_coding:
                c_pos = pdata.n_pos
                interval = pdata.n_interval

            if c_pos is None:                
                var_str += " "
                continue

            if len(var_str) > count:
                continue

            if (c_pos + 1) % 10 == 0:
                # if pdata.c_interval.start.datum == Datum.CDS_END:
                #     var_str += "*"
                var_str += f"{interval} "
                continue

            elif c_pos == 0:
                var_str += f"{interval} "
                continue
            var_str += " "

        return var_str.rstrip()