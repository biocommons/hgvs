from dataclasses import dataclass

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

    def __repr__(self) -> str:
        return f"{self.c_pos}\t{self.tlc}\t{self.aa_char}\t{self.c_pos %3}"


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

    def __repr__(self) -> str:
        return (
            f"{self.mapped_pos}\t{self.c_pos}\t{self.c_offset}\t{self.c_interval}\t"
            + f"{self.cigar_ref}\t"
            + f"{self.chromosome_pos}\t{self.ref}\t{self.alignment_pos}\t{self.mapped_pos_offset}\t"
            + f"\t{self.n_pos}\t{self.tx}"
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
    var_c: SequenceVariant = None
    position_details: list[PositionDetail] = None


@dataclass
class PrettyConfig:
    """acontainer for various configurations."""

    hdp: hgvs.dataproviders.interface.Interface
    am37: AssemblyMapper
    am38: AssemblyMapper
    padding_left: int
    padding_right: int
    default_assembly: str = "GRCh37"
    useColor: bool = False
    showLegend: bool = True
    infer_hgvs_c: bool = True
