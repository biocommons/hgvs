r"""translate between HGVS and other formats

Writing VCF rows, with large del/dup/inv as symbolic alleles::

    babelfish = Babelfish(hdp, assembly_name="GRCh38", symbolic_alt_min_length=50)
    for hgvs_string in ["NC_000017.11:g.43045712_43045713del", "NC_000017.11:g.43044295_43125483dup"]:
        vc = babelfish.hgvs_to_vcf_coordinate(parser.parse(hgvs_string))
        chrom = babelfish.vcf_chrom(vc)
        info = ";".join(f"{k}={v}" for k, v in babelfish.vcf_info(vc).items()) or "."
        print("\t".join([chrom, str(vc.position), hgvs_string, vc.ref, vc.alt, ".", ".", info]))

    17  43045711  NC_000017.11:g.43045712_43045713del  GTA  G      .  .  .
    17  43044294  NC_000017.11:g.43044295_43125483dup  G    <DUP>  .  .  SVTYPE=DUP;SVLEN=81189;END=43125483

and reading them back::

    # NC_000017.11:g.43044295_43125483dup
    babelfish.vcf_to_g_hgvs("17", 43044294, "G", "<DUP>", svlen=81189)

A VCFCoordinate holds the sequence accession, as the rest of hgvs does. Babelfish converts it to
and from the VCF CHROM name: the chromosome name ("17") by default, or the accession itself with
chrom_style="accession". The VCF header needs ##ALT lines for DEL/DUP/INV and ##INFO lines for
SVTYPE, SVLEN and END. SVLEN follows vcf_version (default 4.3, where a deletion's SVLEN is negative).
"""

from dataclasses import dataclass, replace

from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement

import hgvs
import hgvs.normalizer
from hgvs.edit import Dup, Inv, NARefAlt
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
from hgvs.utils.position import get_start_end_interbase

# Positions are 1-based, as in VCF and HGVS. An _i suffix is 0-based interbase, as hdp.get_seq takes

# HGVS edit type <-> VCF symbolic ALT allele
SYMBOLIC_DEL = "<DEL>"
SYMBOLIC_DUP = "<DUP>"
SYMBOLIC_INV = "<INV>"
SYMBOLIC_ALTS = {"del": SYMBOLIC_DEL, "dup": SYMBOLIC_DUP, "inv": SYMBOLIC_INV}
SYMBOLIC_ALT_EDITS = {SYMBOLIC_DEL: NARefAlt, SYMBOLIC_DUP: Dup, SYMBOLIC_INV: Inv}
# Normalizing fetches the whole sequence and holds several copies of it (~8 bytes/base: 1.8GB
# for a 223Mb dup), so the largest symbolic alleles are written as called instead
DEFAULT_NORMALIZE_MAX_LENGTH = 1_000_000
# SVLEN of a deletion is negative up to VCF 4.3, and positive from 4.4
DEFAULT_VCF_VERSION = "4.3"
# VCF CHROM as the chromosome name ("17") or the sequence accession ("NC_000017.11")
CHROM_STYLES = ("name", "accession")


def _check_vcf_version(vcf_version):
    """(major, minor) of a VCF version string like "4.3" """
    try:
        major, minor = (int(v) for v in vcf_version.split("."))
    except (AttributeError, ValueError):
        msg = f"vcf_version must be like '4.3', not {vcf_version!r}"
        raise ValueError(msg) from None
    return major, minor


@dataclass(frozen=True)
class VCFCoordinate:
    """A VCF record's POS, REF, ALT and END (VCF INFO END - the last base affected) on a
    sequence accession (ac) rather than a CHROM name, which Babelfish.vcf_chrom provides.
    END is derived for explicit alleles, and required for a symbolic allele (eg <DEL>), whose
    ref is the preceding (padding) base."""

    ac: str
    position: int
    ref: str
    alt: str
    end: int | None = None

    def __post_init__(self):
        if self.is_symbolic:
            if self.end is None:
                msg = f"Symbolic allele {self.alt} requires end or svlen"
                raise ValueError(msg)
            if self.end <= self.position:
                msg = f"Symbolic allele {self.alt} end ({self.end}) must be after position ({self.position})"
                raise ValueError(msg)
        else:
            end = self.position + len(self.ref) - 1
            if self.end is None:
                object.__setattr__(self, "end", end)  # frozen
            elif self.end != end:
                msg = f"{self.ac}:{self.position} {self.ref} ends at {end}, not {self.end}"
                raise ValueError(msg)

    @classmethod
    def from_svlen(cls, ac, position, ref, alt, svlen):
        """Either sign of SVLEN, as a record doesn't say which VCF version it follows"""
        return cls(ac, position, ref, alt, end=position + abs(svlen))

    @property
    def is_symbolic(self):
        return self.alt.startswith("<")

    @property
    def svlen(self):
        """Length of a symbolic allele, otherwise None. Always positive: Babelfish.vcf_info
        gives SVLEN with the sign of its vcf_version"""
        if not self.is_symbolic:
            return None
        return self.end - self.position

    def as_tuple(self):
        return self.ac, self.position, self.ref, self.alt, self.end

    def format(self):
        if self.is_symbolic:
            return f"{self.ac}:{self.position}-{self.end} {self.alt}"
        return f"{self.ac}:{self.position} {self.ref}>{self.alt}"

    def __str__(self):
        return self.format()


def _common_prefix_len(a: str, b: str) -> int:
    """Return the length of the common character prefix of two strings.

    Replaces the deprecated os.path.commonprefix for sequence strings.
    """
    return next(
        (i for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y), min(len(a), len(b))
    )


class Babelfish:
    def __init__(
        self,
        hdp,
        assembly_name,
        symbolic_alt_min_length=None,
        normalize_max_length=None,
        *,
        vcf_version=DEFAULT_VCF_VERSION,
        chrom_style="name",
    ):
        """symbolic_alt_min_length: a del, dup or inv spanning at least this many bases is
        written to and read from VCF as a symbolic allele (<DEL>, <DUP>, <INV>), and smaller
        symbolic alleles are read as explicit sequence. Structural variants are conventionally
        >= 50bp. None (the default) never writes symbolic alleles, and reads them as given.

        normalize_max_length: a symbolic allele spanning at least this many bases is not
        normalized: shuffling the full sequence of a large structural variant is
        impractically slow, and structural variant callers don't left-align symbolic alleles
        either. Must be at least symbolic_alt_min_length. Defaults to the larger of
        symbolic_alt_min_length and DEFAULT_NORMALIZE_MAX_LENGTH.

        vcf_version: of VCF written, which sets the sign of a deletion's SVLEN from vcf_info
        (negative up to VCF 4.3, positive from 4.4). SVLEN of either sign is read.

        chrom_style: VCF CHROM written by hgvs_to_vcf and vcf_chrom - "name" (the default) for
        the chromosome name without a "chr" prefix (eg "17"), or "accession" for the sequence
        accession (eg "NC_000017.11"). Either is read.
        """
        _check_vcf_version(vcf_version)
        if chrom_style not in CHROM_STYLES:
            msg = f"chrom_style must be one of {CHROM_STYLES}, not {chrom_style!r}"
            raise ValueError(msg)
        if symbolic_alt_min_length is not None:
            if normalize_max_length is None:
                normalize_max_length = max(symbolic_alt_min_length, DEFAULT_NORMALIZE_MAX_LENGTH)
            elif normalize_max_length < symbolic_alt_min_length:
                msg = "normalize_max_length must be at least symbolic_alt_min_length"
                raise ValueError(msg)
        elif normalize_max_length is not None:
            msg = "normalize_max_length requires symbolic_alt_min_length"
            raise ValueError(msg)
        self.assembly_name = assembly_name
        self.symbolic_alt_min_length = symbolic_alt_min_length
        self.normalize_max_length = normalize_max_length
        self.vcf_version = vcf_version
        self.chrom_style = chrom_style
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(
            hdp, cross_boundaries=False, shuffle_direction=5, validate=False
        )
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)
        # We need to accept accessions as chromosome names, so add them pointing at themselves
        self.name_to_ac_map.update({ac: ac for ac in self.name_to_ac_map.values()})

    def hgvs_to_vcf(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, pos, ref, alt, type) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix), unless chrom_style="accession".

        Symbolic alleles have no END here - use hgvs_to_vcf_coordinate
        """
        vc, typ = self._hgvs_to_vcf(var_g)
        return self.vcf_chrom(vc), vc.position, vc.ref, vc.alt, typ

    def hgvs_to_vcf_coordinate(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to a VCFCoordinate. A del, dup or inv of at least
        symbolic_alt_min_length bases (see __init__) is returned as a symbolic allele.
        """
        vc, _ = self._hgvs_to_vcf(var_g)
        return vc

    def vcf_chrom(self, vc):
        """VCF CHROM for a VCFCoordinate in this Babelfish's chrom_style (see __init__)"""
        if self.chrom_style == "accession":
            return vc.ac
        return self.ac_to_name_map[vc.ac]

    def vcf_info(self, vc):
        """VCF INFO fields a symbolic allele needs (SVTYPE, SVLEN, END), empty otherwise.
        A deletion's SVLEN is negative before vcf_version 4.4"""
        if not vc.is_symbolic:
            return {}
        svlen = vc.svlen
        if vc.alt.upper() == SYMBOLIC_DEL and _check_vcf_version(self.vcf_version) < (4, 4):
            svlen = -svlen
        return {"SVTYPE": vc.alt.strip("<>"), "SVLEN": svlen, "END": vc.end}

    def vcf_coordinate(self, chrom, position, ref, alt, *, end=None, svlen=None):
        """VCFCoordinate from a VCF record. chrom is a chromosome name or accession. Symbolic
        alleles <DEL>, <DUP> and <INV> require end (VCF INFO END - the last base affected) or
        svlen (VCF INFO SVLEN, of either sign)"""
        ac = self.name_to_ac_map[chrom]
        if svlen is not None and end is None:
            return VCFCoordinate.from_svlen(ac, position, ref, alt, svlen)
        return VCFCoordinate(ac, position, ref, alt, end=end)

    def _hgvs_to_vcf(self, var_g):
        """(VCFCoordinate, HGVS edit type)"""
        if var_g.type != "g":
            msg = f"Expected g. variant, got {var_g}"
            raise RuntimeError(msg)

        if self._too_big_to_normalize(self.symbolic_alt_length(var_g)):
            return self._hgvs_to_vcf_symbolic(var_g)

        vleft = self.hn.normalize(var_g)
        if self._is_symbolic_length(self.symbolic_alt_length(vleft)):
            return self._hgvs_to_vcf_symbolic(vleft)
        return self._hgvs_to_vcf_explicit(vleft)

    def _hgvs_to_vcf_explicit(self, vleft):
        """(VCFCoordinate, HGVS edit type) with explicit ref/alt, for a normalized variant, or
        one too big to normalize at the start of the sequence"""
        # We are taking the inner interval, but plan on implementing INFO fields in issue #788
        start_i, end_i = get_start_end_interbase(vleft.posedit.pos, outer_confidence=False)
        typ = vleft.posedit.edit.type

        if typ == "dup":
            if start_i == 0:
                # The copy is inserted before the original
                copy = self.hdp.get_seq(vleft.ac, 0, end_i)
                return self._vcf_padded_after(vleft.ac, 0, 0, copy), typ
            start_i -= 1
            alt = self.hdp.get_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
        elif typ == "inv":
            # Unnormalized (see _hgvs_to_vcf_symbolic) an inv may have no ref
            ref = vleft.posedit.edit.ref or self.hdp.get_seq(vleft.ac, start_i, end_i)
            alt = reverse_complement(ref)
        else:
            alt = vleft.posedit.edit.alt or ""

            if typ in ("del", "ins"):
                if typ == "ins":
                    # ins coordinates (only) exclude left position
                    start_i += 1
                    end_i -= 1
                if start_i == 0:
                    return self._vcf_padded_after(vleft.ac, start_i, end_i, alt), typ
                # Left anchored
                start_i -= 1
                ref = self.hdp.get_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = "."
        return VCFCoordinate(vleft.ac, start_i + 1, ref, alt), typ

    def _vcf_padded_after(self, ac, start_i, end_i, alt):
        """VCFCoordinate replacing start_i:end_i with alt at the start of the sequence, which
        VCF pads with the base after rather than the base before"""
        ref = self.hdp.get_seq(ac, start_i, end_i + 1)
        return VCFCoordinate(ac, start_i + 1, ref, alt + ref[-1])

    @staticmethod
    def symbolic_alt_length(var_g):
        """Length of a del, dup or inv (ie that can be represented as a symbolic allele),
        otherwise None. Also None when an explicit sequence doesn't fit the interval, which
        is left for the normalizer to reject."""
        edit = var_g.posedit.edit
        if edit.type not in SYMBOLIC_ALTS:
            return None
        start_i, end_i = get_start_end_interbase(var_g.posedit.pos, outer_confidence=False)
        length = end_i - start_i
        if edit.ref_s is not None and len(edit.ref_s) != length:
            return None
        return length

    def _is_symbolic_length(self, length):
        min_length = self.symbolic_alt_min_length
        return length is not None and min_length is not None and length >= min_length

    def _too_big_to_normalize(self, length):
        max_length = self.normalize_max_length
        return length is not None and max_length is not None and length >= max_length

    def _hgvs_to_vcf_symbolic(self, var_g):
        start_i, end_i = get_start_end_interbase(var_g.posedit.pos, outer_confidence=False)
        typ = var_g.posedit.edit.type
        padding_i = start_i - 1  # The base before the variant
        if padding_i < 0:
            # No padding base before a variant at the start of the sequence
            return self._hgvs_to_vcf_explicit(var_g)
        ref = self.hdp.get_seq(var_g.ac, padding_i, start_i)
        return VCFCoordinate(var_g.ac, padding_i + 1, ref, SYMBOLIC_ALTS[typ], end=end_i), typ

    def vcf_to_g_hgvs(self, chrom, position, ref, alt, *, end=None, svlen=None):
        """chrom is a chromosome name or accession. Symbolic alleles <DEL>, <DUP> and <INV>
        require end (VCF INFO END - the last base affected) or svlen (VCF INFO SVLEN).
        See vcf_coordinate_to_g_hgvs"""
        vc = self.vcf_coordinate(chrom, position, ref, alt, end=end, svlen=svlen)
        return self.vcf_coordinate_to_g_hgvs(vc)

    def vcf_coordinate_to_g_hgvs(self, vc):
        """With symbolic_alt_min_length set (see __init__), symbolic alleles below
        normalize_max_length are expanded to explicit sequence and normalized, and explicit
        del, dup or inv at or above it are read as symbolic without normalizing, so each
        variant has one HGVS form. A del, dup or inv of at least symbolic_alt_min_length is
        returned without its sequence.
        """
        ac = vc.ac
        # VCF spec https://samtools.github.io/hts-specs/VCFv4.1.pdf
        # says for REF/ALT "Each base must be one of A,C,G,T,N (case insensitive)"
        vc = replace(vc, ref=vc.ref.upper(), alt=vc.alt.upper())

        if vc.is_symbolic:
            self._check_symbolic_alt(vc.alt)
            if self.symbolic_alt_min_length is None or self._too_big_to_normalize(vc.svlen):
                return self._symbolic_g_hgvs(ac, vc.alt, vc.position + 1, vc.end)
            vc = self.vcf_coordinate_as_explicit(vc)
        elif self.normalize_max_length is not None:
            symbolic = self.vcf_coordinate_as_symbolic(vc, min_length=self.normalize_max_length)
            if symbolic.is_symbolic:
                return self._symbolic_g_hgvs(ac, symbolic.alt, symbolic.position + 1, symbolic.end)

        n = Normalizer(self.hdp)
        var_g = n.normalize(self._explicit_vcf_to_g_hgvs(ac, vc.position, vc.ref, vc.alt))
        if self._is_symbolic_length(self.symbolic_alt_length(var_g)):
            # Drop the sequence, as if read from a symbolic allele
            start_i, end_i = get_start_end_interbase(var_g.posedit.pos, outer_confidence=False)
            symbolic_alt = SYMBOLIC_ALTS[var_g.posedit.edit.type]
            var_g = self._symbolic_g_hgvs(ac, symbolic_alt, start_i + 1, end_i)
        return var_g

    @staticmethod
    def _explicit_vcf_to_g_hgvs(ac, position, ref, alt):
        """Unnormalized g. variant for an explicit (non-symbolic) uppercase ref/alt"""
        if ref != alt:
            # Strip common prefix
            if len(alt) > 1 and len(ref) > 1:
                # Character-wise common prefix (see _common_prefix_len).
                lp = _common_prefix_len(ref, alt)
                if lp > 0:
                    ref = ref[lp:]
                    alt = alt[lp:]
                    position += lp
            elif alt == ".":
                alt = ref

        if ref == "":  # Insert
            # Insert uses coordinates around the insert point.
            start = position - 1
            end = position
        else:
            start = position
            end = position + len(ref) - 1

        return SequenceVariant(
            ac=ac,
            type="g",
            posedit=PosEdit(
                Interval(
                    start=SimplePosition(start),
                    end=SimplePosition(end),
                    uncertain=False,
                ),
                NARefAlt(ref=ref or None, alt=alt or None, uncertain=False),
            ),
        )

    def vcf_coordinate_as_explicit(self, vc):
        """Explicit ref/alt for a symbolic <DEL>, <DUP> or <INV>, otherwise vc. Not normalized"""
        if not vc.is_symbolic:
            return vc
        alt = vc.alt.upper()
        self._check_symbolic_alt(alt)
        padding_i = vc.position - 1
        seq = self.hdp.get_seq(vc.ac, padding_i, vc.end)  # padding base onwards
        if alt == SYMBOLIC_DEL:
            return VCFCoordinate(vc.ac, vc.position, seq, seq[0])
        if alt == SYMBOLIC_DUP:
            return VCFCoordinate(vc.ac, vc.position, seq[0], seq)
        # An explicit inversion has no padding base
        inverted = seq[1:]
        return VCFCoordinate(vc.ac, vc.position + 1, inverted, reverse_complement(inverted))

    def vcf_coordinate_as_symbolic(self, vc, min_length=None):
        """<DEL>, <DUP> or <INV> for an explicit del, dup or inv of at least min_length
        (default symbolic_alt_min_length) bases, otherwise vc. Not normalized"""
        if min_length is None:
            min_length = self.symbolic_alt_min_length
        if min_length is None:
            msg = "vcf_coordinate_as_symbolic requires min_length or symbolic_alt_min_length"
            raise ValueError(msg)
        if vc.is_symbolic:
            return vc

        ac, position = vc.ac, vc.position
        ref, alt = vc.ref.upper(), vc.alt.upper()
        if len(alt) == 1 and len(ref) > min_length and ref[0] == alt:
            return VCFCoordinate(ac, position, alt, SYMBOLIC_DEL, end=position + len(ref) - 1)

        if len(ref) == 1 and len(alt) > min_length and alt[0] == ref:
            return self._insertion_as_symbolic_dup(vc) or vc

        if len(ref) == len(alt) >= min_length and ref == reverse_complement(alt):
            # Symbolic alts are padded with the base before
            padding_position = position - 1
            padding_i = padding_position - 1
            if padding_i < 0:
                # No padding base before an inversion at the start of the sequence
                return vc
            padding = self.hdp.get_seq(ac, padding_i, padding_position)
            return VCFCoordinate(
                ac, padding_position, padding, SYMBOLIC_INV, end=position + len(ref) - 1
            )
        return vc

    def _insertion_as_symbolic_dup(self, vc):
        """<DUP> if the inserted sequence duplicates the bases either side, otherwise None"""
        ac, position = vc.ac, vc.position
        inserted = vc.alt[1:].upper()
        length = len(inserted)
        insertion_i = position  # Between the ref base at position and the next
        # Left aligned, the inserted copy precedes the bases it duplicates
        if self.hdp.get_seq(ac, insertion_i, insertion_i + length) == inserted:
            return VCFCoordinate(ac, position, vc.ref.upper(), SYMBOLIC_DUP, end=position + length)
        # Otherwise a preceding copy ends at position, padded by the base at padding_position
        padding_position = position - length
        padding_i = padding_position - 1
        if padding_i < 0:
            # No room for a preceding copy this close to the start of the sequence
            return None
        padded = self.hdp.get_seq(ac, padding_i, insertion_i)
        if padded[1:] == inserted:
            return VCFCoordinate(ac, padding_position, padded[0], SYMBOLIC_DUP, end=position)
        return None

    @staticmethod
    def _check_symbolic_alt(alt):
        if alt not in SYMBOLIC_ALT_EDITS:
            msg = f"Unsupported symbolic allele {alt}"
            raise ValueError(msg)

    @staticmethod
    def _symbolic_g_hgvs(ac, alt, start, end):
        """g. del, dup or inv of start_end (1-based) without its sequence, from a symbolic alt"""
        edit = SYMBOLIC_ALT_EDITS[alt](ref="", uncertain=False)
        return SequenceVariant(
            ac=ac,
            type="g",
            posedit=PosEdit(
                Interval(
                    start=SimplePosition(start),
                    end=SimplePosition(end),
                    uncertain=False,
                ),
                edit,
            ),
        )
