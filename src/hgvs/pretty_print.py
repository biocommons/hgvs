from typing import List

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig
from hgvs.pretty.renderer.chrom_seq_renderer import ChromSeqRendered
from hgvs.pretty.renderer.chrom_seq_reverse_renderer import ChromReverseSeqRendered
from hgvs.pretty.renderer.pos_info import ChrPositionInfo
from hgvs.pretty.renderer.prot_mapping_renderer import ProtMappingRenderer
from hgvs.pretty.renderer.prot_ruler_renderer import ProtRulerRenderer
from hgvs.pretty.renderer.prot_seq_renderer import ProtSeqRenderer
from hgvs.pretty.renderer.ruler import ChrRuler
from hgvs.pretty.renderer.shuffled_variant import ShuffledVariant
from hgvs.pretty.renderer.tx_alig_renderer import TxAligRenderer
from hgvs.pretty.renderer.tx_mapping_renderer import TxMappingRenderer
from hgvs.pretty.renderer.tx_pos import TxRulerRenderer
from hgvs.pretty.renderer.tx_ref_disagree_renderer import TxRefDisagreeRenderer
from hgvs.sequencevariant import SequenceVariant

TGREEN = "\033[32m"  # Green Text
TGREENBG = "\033[30;42m"
TRED = "\033[31m"  # Red Text
TREDBG = "\033[30;41m"
TBLUE = "\033[34m"  # Blue Text
TBLUEBG = "\033[30;44m"
TPURPLE = "\033[35m"  # Purple Text
TPURPLEBG = "\033[30;45m"
TYELLOW = "\033[33m"  # Yellow Text
TYELLOWBG = "\033[30;43m"

ENDC = "\033[m"  # reset to the defaults


class PrettyPrint:
    """A class that provides a pretty display of the genomic context of a variant."""

    def __init__(
        self,
        hdp: hgvs.dataproviders.interface.Interface,
        default_assembly: str = "GRCh37",
        padding_left: int = 20,
        padding_right: int = 20,
        useColor=False,
        showLegend=True,
        infer_hgvs_c=True,
        all=False,
        show_reverse_strand=False,
    ):
        """
        :param hdp: HGVS Data Provider Interface-compliant instance
           (see :class:`hgvs.dataproviders.interface.Interface`)

        :param padding: spacing left and right of the variant for display purposes.
        """
        am37: AssemblyMapper = AssemblyMapper(hdp, assembly_name="GRCh37")
        am38: AssemblyMapper = AssemblyMapper(hdp, assembly_name="GRCh38")

        self.config = PrettyConfig(
            hdp,
            am37,
            am38,
            padding_left,
            padding_right,
            default_assembly,
            useColor,
            showLegend,
            infer_hgvs_c,
            all,
            show_reverse_strand,
        )

    def _get_assembly_mapper(self) -> AssemblyMapper:
        if self.config.default_assembly == "GRCh37":
            am = self.config.am37
        else:
            am = self.config.am38

        return am

    def _get_all_transcripts(self, var_g) -> List[str]:
        am = self._get_assembly_mapper()

        transcripts = am.relevant_transcripts(var_g)

        return transcripts

    def _infer_hgvs_c(self, var_g: SequenceVariant, tx_ac: str = None) -> SequenceVariant:

        if not tx_ac:
            transcripts = self._get_all_transcripts(var_g)
            if transcripts:
                tx_ac = transcripts[0]
            else:
                return None

        am = self._get_assembly_mapper()

        if tx_ac.startswith("NR_"):
            var_n = am.g_to_n(var_g, tx_ac)
            return var_n
        var_c = am.g_to_c(var_g, tx_ac)
        return var_c

    def _map_to_chrom(self, sv: SequenceVariant) -> SequenceVariant:
        """maps a variant to chromosomal coords, if needed."""
        if self.config.default_assembly == "GRCh37":
            am = self.config.am37
        else:
            am = self.config.am38

        if sv.type == "c":
            return am.c_to_g(sv)
        elif sv.type == "n":
            return am.n_to_g(sv)
        elif sv.type == "t":
            return am.t_to_g(sv)

    def _colorize_hgvs(self, hgvs_str: str) -> str:
        if not self.config.useColor:
            return hgvs_str

        spl = hgvs_str.split(":")
        var_str = TPURPLE + spl[0] + ENDC
        var_str += ":"

        sec = spl[1].split(".")
        var_str += TYELLOW + sec[0] + ENDC
        var_str += "."
        var_str += sec[1]

        return var_str

    def display(
        self,
        sv: SequenceVariant,
        tx_ac: str = None,
        display_start: int = None,
        display_end: int = None,
    ) -> str:
        """Takes a variant and prints the genomic context around it."""

        var_c_or_n = None
        if sv.type == "g":
            var_g = sv
            if tx_ac is not None:
                var_c_or_n = self._infer_hgvs_c(var_g, tx_ac=tx_ac)
        elif sv.type == "c":
            var_g = self._map_to_chrom(sv)
            var_c_or_n = sv
        elif sv.type == "n":
            # map back to genome
            var_g = self._map_to_chrom(sv)
            var_c_or_n = sv

        data_compiler = DataCompiler(config=self.config)

        if self.config.all:
            # get all overlapping transcripts

            response = ""
            tx_acs = self._get_all_transcripts(var_g)
            print(f"displaying {len(tx_acs)} alternative transcripts")
            for tx_ac in tx_acs:
                var_c_or_n = self._infer_hgvs_c(var_g, tx_ac)
                response += self.create_repre(
                    var_g, var_c_or_n, display_start, display_end, data_compiler
                )
                response += "\n---\n"
            return response
        else:

            if not var_c_or_n:
                if self.config.infer_hgvs_c:
                    var_c_or_n = self._infer_hgvs_c(var_g)

            return self.create_repre(var_g, var_c_or_n, display_start, display_end, data_compiler)

    def create_repre(
        self,
        var_g: SequenceVariant,
        var_c_or_n: SequenceVariant,
        display_start: int,
        display_end: int,
        data_compiler: DataCompiler,
    ):
        data = data_compiler.data(var_g, var_c_or_n, display_start, display_end)

        left_shuffled_var = data.left_shuffled
        right_shuffled_var = data.right_shuffled
        fully_justified_var = data.fully_justified

        if self.config.showLegend:
            head = "hgvs      : "
            refa = "ref>alt   : "
        else:
            head = refa = ""

        if self.config.useColor:
            var_g_print = self._colorize_hgvs(str(var_g))
        else:
            var_g_print = str(var_g)

        var_str = head + var_g_print + "\n"
        if data.var_c_or_n:
            if self.config.useColor:
                var_c_print = self._colorize_hgvs(str(var_c_or_n))
            else:
                var_c_print = str(var_c_or_n)
            var_str += head + var_c_print + "\n"

        renderer_cls = [ChrPositionInfo, ChrRuler, ChromSeqRendered]

        if self.config.show_reverse_strand:
            renderer_cls.append(ChromReverseSeqRendered)

        renderer_cls.append(TxRefDisagreeRenderer)

        renderers = []
        for cls in renderer_cls:
            r = cls(self.config, data.strand)
            renderers.append(r)

        for renderer in renderers:
            d = ""
            if self.config.showLegend:
                d += renderer.legend()
            str_results = renderer.display(data)

            if str_results:
                var_str += d + str_results + "\n"

        left_shuffled_renderer = ShuffledVariant(self.config, data.strand, var_g, left_shuffled_var)
        left_shuffled_str = left_shuffled_renderer.display(data)

        right_shuffled_renderer = ShuffledVariant(
            self.config, data.strand, var_g, right_shuffled_var
        )
        right_shuffled_str = right_shuffled_renderer.display(data)
        if self.config.showLegend:
            shuffled_seq_header = left_shuffled_renderer.legend()
        else:
            shuffled_seq_header = ""

        if left_shuffled_str != right_shuffled_str:
            fully_justified_renderer = ShuffledVariant(
                self.config, data.strand, var_g, fully_justified_var
            )
            fully_justified_str = fully_justified_renderer.display(data)

            # var_str += shuffled_seq_header + left_shuffled_str + "\n"
            # var_str += shuffled_seq_header + right_shuffled_str + "\n"
            var_str += shuffled_seq_header + fully_justified_str + "\n"

        else:
            var_str += shuffled_seq_header + left_shuffled_str + "\n"

        renderers_cls = [TxAligRenderer, TxMappingRenderer, TxRulerRenderer, ProtSeqRenderer, ProtMappingRenderer, ProtRulerRenderer ]
        for cls in renderers_cls:
            renderer = cls(self.config, data.strand)

            d = ""
            if self.config.showLegend:
                d += renderer.legend()
            str_results = renderer.display(data)

            if str_results:
                var_str += d + str_results + "\n"

        if left_shuffled_str != right_shuffled_str:
            # TODO: detect repeats?
            fully_justified_ref = fully_justified_var.ref
            fully_justified_alt = fully_justified_var.alt
            var_str += refa + f"{fully_justified_ref}>{fully_justified_alt}"

        return var_str
