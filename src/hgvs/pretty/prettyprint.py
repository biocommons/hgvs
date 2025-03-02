from typing import List, Tuple

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from hgvs.pretty.console.renderer import colorize_hgvs
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig
from hgvs.pretty.console.chromseqrenderer import ChromSeqRendered
from hgvs.pretty.console.chromseqreverserenderer import ChromReverseSeqRendered
from hgvs.pretty.console.posinfo import ChrPositionInfo
from hgvs.pretty.console.protmappingrenderer import ProtMappingRenderer
from hgvs.pretty.console.protrulerrenderer import ProtRulerRenderer
from hgvs.pretty.console.protseqrenderer import ProtSeqRenderer
from hgvs.pretty.console.ruler import ChrRuler
from hgvs.pretty.console.shuffledvariant import RegionImpacted
from hgvs.pretty.console.txaligrenderer import TxAligRenderer
from hgvs.pretty.console.txmappingrenderer import TxMappingRenderer
from hgvs.pretty.console.txpos import TxRulerRenderer
from hgvs.pretty.console.txrefdisagreerenderer import TxRefDisagreeRenderer
from hgvs.repeats import RepeatAnalyser
from hgvs.sequencevariant import SequenceVariant


class PrettyPrint:
    """A class to handle the pretty printing of HGVS variants with various configurations.

    Attributes:
        config (PrettyConfig): Configuration object for pretty printing.
    Methods:
        __init__(hdp, assembly_mapper, padding_left=20, padding_right=20, use_color=False, show_legend=True, infer_hgvs_c=True, all=False, show_reverse_strand=False, alt_aln_method='splign', reverse_display=True):
            Initializes the PrettyPrint object with the given parameters.
        _get_all_transcripts(var_g) -> List[str]:
            Retrieves all relevant transcripts for a given genomic variant.
        _infer_hgvs_c(var_g: SequenceVariant, tx_ac: str = None) -> SequenceVariant:
            Infers the HGVS coding or non-coding variant from a genomic variant.
        _map_to_chrom(sv: SequenceVariant) -> SequenceVariant:
            Maps a variant to chromosomal coordinates if needed.
        get_hgvs_names(sv: SequenceVariant, tx_ac: str = None) -> Tuple[SequenceVariant, SequenceVariant]:
            Retrieves the HGVS names for a given sequence variant.
        display(sv: SequenceVariant, tx_ac: str = None, display_start: int = None, display_end: int = None) -> str:
            Takes a variant and prints the genomic context around it.
        create_repre(var_g: SequenceVariant, var_c_or_n: SequenceVariant, display_start: int, display_end: int, data_compiler: DataCompiler) -> str:
            Creates a representation of the variant with the given parameters.
    """

    def __init__(
        self,
        hdp: hgvs.dataproviders.interface.Interface,
        assembly_mapper: AssemblyMapper,
        padding_left: int = 20,
        padding_right: int = 20,
        use_color=False,
        show_legend=True,
        infer_hgvs_c=True,
        all=False,
        show_reverse_strand=False,
        alt_aln_method="splign",
        reverse_display=True,
    ):
        """
        :param hdp: HGVS Data Provider Interface-compliant instance
           (see :class:`hgvs.dataproviders.interface.Interface`)

        :param padding: spacing left and right of the variant for display purposes.
        """

        self.config = PrettyConfig(
            hdp=hdp,
            assembly_mapper=assembly_mapper,
            padding_left=padding_left,
            padding_right=padding_right,
            use_color=use_color,
            show_legend=show_legend,
            infer_hgvs_c=infer_hgvs_c,
            all=all,
            show_reverse_strand=show_reverse_strand,
            alt_aln_method=alt_aln_method,
            reverse_display=reverse_display,
        )

    def _get_all_transcripts(self, var_g) -> List[str]:
        transcripts = self.config.assembly_mapper.relevant_transcripts(var_g)

        return transcripts

    def _infer_hgvs_c(
        self, var_g: SequenceVariant, tx_ac: str = None
    ) -> SequenceVariant:
        if not tx_ac:
            transcripts = self._get_all_transcripts(var_g)
            if transcripts:
                tx_ac = transcripts[0]
            else:
                return None

        am = self.config.assembly_mapper

        if tx_ac.startswith("NR_"):
            var_n = am.g_to_n(var_g, tx_ac)
            return var_n
        var_c = am.g_to_c(var_g, tx_ac)
        return var_c

    def _map_to_chrom(self, sv: SequenceVariant) -> SequenceVariant:
        """maps a variant to chromosomal coords, if needed."""
        am = self.config.assembly_mapper

        if sv.type == "c":
            return am.c_to_g(sv)
        elif sv.type == "n":
            return am.n_to_g(sv)
        elif sv.type == "t":
            return am.t_to_g(sv)
        elif sv.type == "r":
            return am.r

    def get_hgvs_names(
        self, sv: SequenceVariant, tx_ac: str = None
    ) -> Tuple[SequenceVariant, SequenceVariant]:
        var_c_or_n = None
        if sv.type == "g":
            var_g = sv
            if tx_ac is not None:
                var_c_or_n = self._infer_hgvs_c(var_g, tx_ac=tx_ac)
        elif sv.type == "c":
            var_g = self._map_to_chrom(sv)
            var_c_or_n = sv
        elif sv.type == "r":
            var_g = self._map_to_chrom(sv)
            var_c_or_n = sv
        elif sv.type == "n":
            # map back to genome
            var_g = self._map_to_chrom(sv)
            var_c_or_n = sv

        return var_g, var_c_or_n

    def display(
        self,
        sv: SequenceVariant,
        tx_ac: str = None,
        display_start: int = None,
        display_end: int = None,
    ) -> str:
        """Takes a variant and prints the genomic context around it in a string representation."""

        var_g, var_c_or_n = self.get_hgvs_names(sv, tx_ac)

        self.data_compiler = DataCompiler(config=self.config)

        if self.config.all:
            # get all overlapping transcripts

            response = ""
            tx_acs = self._get_all_transcripts(var_g)
            for tx_ac in tx_acs:
                var_c_or_n = self._infer_hgvs_c(var_g, tx_ac)
                response += self.create_repre(
                    var_g, var_c_or_n, display_start, display_end, self.data_compiler
                )
                response += "\n---\n"
            return response
        else:
            if not var_c_or_n:
                if self.config.infer_hgvs_c:
                    var_c_or_n = self._infer_hgvs_c(var_g)

            return self.create_repre(
                var_g, var_c_or_n, display_start, display_end, self.data_compiler
            )

    def create_repre(
        self,
        var_g: SequenceVariant,
        var_c_or_n: SequenceVariant,
        display_start: int,
        display_end: int,
        data_compiler: DataCompiler,
    ):
        """
        Creates a string representation of sequence variants with optional colorization and legends.

        The string representation is built by calling various renderers to display the genomic, coding, and protein
        sequences, as well as the alignment of the sequences. The display is formatted with optional colorization.

        Configuration is possible through the PrettyConfig object, which has been passed into the constructor.

        Args:
            var_g (SequenceVariant): Genomic sequence variant.
            var_c_or_n (SequenceVariant): Coding or non-coding sequence variant.
            display_start (int): Start position for display.
            display_end (int): End position for display.
            data_compiler (DataCompiler): Compiler for sequence variant data.
        Returns:
            str: Formatted string representation of the sequence variants.
        """

        data = data_compiler.data(var_g, var_c_or_n, display_start, display_end)

        left_shuffled_var = data.left_shuffled
        right_shuffled_var = data.right_shuffled
        fully_justified_var = data.fully_justified

        if self.config.show_legend:
            head = "hgvs_g    : "
            head_c = "hgvs_c    : "
            head_n = "hgvs_n    : "
            head_p = "hgvs_p    : "
            refa = "ref>alt   : "
        else:
            head = head_c = head_p = refa = ""

        if self.config.use_color:
            var_g_print = colorize_hgvs(str(var_g))
        else:
            var_g_print = str(var_g)

        var_str = head + var_g_print + "\n"
        if data.var_c_or_n:
            if self.config.use_color:
                var_c_print = colorize_hgvs(str(var_c_or_n))
            else:
                var_c_print = str(var_c_or_n)
            if data.var_c_or_n.type == "c":
                var_str += head_c
            elif data.var_c_or_n.type == "n":
                var_str += head_n
            var_str += var_c_print + "\n"

        if data.var_p:
            if self.config.use_color:
                var_p_print = colorize_hgvs(str(data.var_p))
            else:
                var_p_print = str(data.var_p)
            var_str += head_p + var_p_print + "\n"

        renderer_cls = [ChrPositionInfo, ChrRuler, ChromSeqRendered]

        # if we show reverse strand transcripts in forward facing orientation, always
        # show both forward and reverse strand sequences.
        if (
            self.config.show_reverse_strand
            or self.config.reverse_display
            and data.strand < 0
        ):
            renderer_cls.append(ChromReverseSeqRendered)

        renderer_cls.append(TxRefDisagreeRenderer)

        renderers = []
        for cls in renderer_cls:
            r = cls(self.config, data.strand)
            renderers.append(r)

        for renderer in renderers:
            d = ""
            if self.config.show_legend:
                d += renderer.legend()
            str_results = renderer.display(data)

            if str_results:
                var_str += d + str_results + "\n"

        left_shuffled_renderer = RegionImpacted(
            self.config, data.strand, var_g, left_shuffled_var
        )
        left_shuffled_str = left_shuffled_renderer.display(data)

        right_shuffled_renderer = RegionImpacted(
            self.config, data.strand, var_g, right_shuffled_var
        )
        right_shuffled_str = right_shuffled_renderer.display(data)
        if self.config.show_legend:
            shuffled_seq_header = left_shuffled_renderer.legend()
        else:
            shuffled_seq_header = ""

        if left_shuffled_str != right_shuffled_str:
            fully_justified_renderer = RegionImpacted(
                self.config, data.strand, var_g, fully_justified_var
            )
            fully_justified_str = fully_justified_renderer.display(data)
            if self.config.show_all_shuffleable_regions:
                var_str += shuffled_seq_header + left_shuffled_str + "\n"
                var_str += shuffled_seq_header + right_shuffled_str + "\n"
            var_str += shuffled_seq_header + fully_justified_str + "\n"

        else:
            var_str += shuffled_seq_header + left_shuffled_str + "\n"

        renderers_cls = [
            TxAligRenderer,
            TxMappingRenderer,
            TxRulerRenderer,
            ProtSeqRenderer,
            ProtMappingRenderer,
            ProtRulerRenderer,
        ]
        for cls in renderers_cls:
            renderer = cls(self.config, data.strand)

            d = ""
            if self.config.show_legend:
                d += renderer.legend()
            str_results = renderer.display(data)

            if str_results:
                var_str += d + str_results + "\n"

        if left_shuffled_str != right_shuffled_str:
            # TODO: detect repeats?
            fully_justified_ref = fully_justified_var.ref
            fully_justified_alt = fully_justified_var.alt

            ra = RepeatAnalyser(fully_justified_var)
            if ra.is_repeat:
                var_str += refa + str(ra)
            else:
                var_str += refa + f"{fully_justified_ref}>{fully_justified_alt}"

        return var_str
