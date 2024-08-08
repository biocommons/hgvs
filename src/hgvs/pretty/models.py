from dataclasses import dataclass
from typing import List

import hgvs
from hgvs.alignmentmapper import AlignmentMapper
from hgvs.assemblymapper import AssemblyMapper
from hgvs.location import Interval
from hgvs.sequencevariant import SequenceVariant


@dataclass(eq=True, repr=True, frozen=True, order=True)
class VariantCoords:
    """Representation of a variant in one of the shuffled representation"""

    start: int
    end: int
    ref: str
    alt: str


@dataclass(eq=True, frozen=True, order=True)
class ProteinData:
    c_pos: int = None
    aa: str = None  # one letter AA code
    tlc: str = None  # three letter AA code
    aa_char: str = None  # the character from the tlc that maps to the c_pos
    var_p: SequenceVariant = None
    is_init_met: bool = False
    is_stop_codon: bool = False
    aa_pos: int = -1

    @classmethod
    def get_header(cls) -> str:
        return f"c_pos\ttlc\taa_char\tc_pos\taa_pos"

    def __repr__(self) -> str:
        return f"{self.c_pos}\t{self.tlc}\t{self.aa_char}\t{self.c_pos %3}\t{self.aa_pos}"


@dataclass(eq=True, frozen=False, order=True)
class PositionDetail:

    chromosome_pos: int = None
    alignment_pos: int = None
    mapped_pos: int = None
    mapped_pos_offset: int = None
    cigar_ref: str = None
    n_pos: int = None
    c_pos: int = None
    c_offset: int = None
    ref: str = None
    tx: str = None
    n_interval: Interval = None
    c_interval: Interval = None
    protein_data: ProteinData = None
    exon_nr: int = None
    variant_feature: str = None

    @classmethod
    def get_header(cls) -> str:
        return (
            "alignment_pos\tc_pos\tc_offset\tc_interval\tcigar_ref\tchromosome_pos\tref\tmapped_pos\t"
            "mapped_pos_offset\tn_pos\ttx\tvariant_feature\texon_nr"
        )

    def __repr__(self) -> str:
        return (
            f"{self.alignment_pos}\t{self.c_pos}\t{self.c_offset}\t{self.c_interval}\t"
            + f"{self.cigar_ref}\t"
            + f"{self.chromosome_pos}\t{self.ref}\t{self.mapped_pos}\t{self.mapped_pos_offset}\t"
            + f"\t{self.n_pos}\t{self.tx}"
            + f"\t{self.variant_feature} {self.exon_nr}"
        )


@dataclass(eq=True, repr=True, frozen=True, order=True)
class VariantData:
    """A data container for knowledge we have about a variant."""

    display_start: int
    display_end: int
    left_shuffled: VariantCoords
    right_shuffled: VariantCoords
    fully_justified: VariantCoords
    display_seq: str
    tx_seq: str
    alignmentmapper: AlignmentMapper
    var_g: SequenceVariant
    strand: int
    var_c_or_n: SequenceVariant = None
    var_p: SequenceVariant = None
    position_details: List[PositionDetail] = None
    all: bool = False


@dataclass
class PrettyConfig:
    """acontainer for various configurations."""

    hdp: hgvs.dataproviders.interface.Interface
    am37: AssemblyMapper
    am38: AssemblyMapper
    padding_left: int = 20
    padding_right: int = 20
    default_assembly: str = "GRCh37"
    useColor: bool = False
    showLegend: bool = True
    showAllShuffleableRegions = False
    infer_hgvs_c: bool = True
    all: bool = False  # print all possible hgvs_c (for all UTA transcripts)
    show_reverse_strand: bool = False  # show the reverse strand sequence for the chromosome
