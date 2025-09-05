"""translate between HGVS and other formats"""

import os

from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement

import hgvs
import hgvs.normalizer
from hgvs.edit import NARefAlt
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
from hgvs.utils.position import get_start_end_interbase


class Babelfish:
    def __init__(self, hdp, assembly_name):
        self.assembly_name = assembly_name
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

        converts a single hgvs allele to (chr, pos, ref, alt) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).
        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))

        vleft = self.hn.normalize(var_g)

        # We are taking the inner interval, but plan on implementing INFO fields in issue #788
        start_i, end_i = get_start_end_interbase(
            vleft.posedit.pos, outer_confidence=False
        )

        chrom = self.ac_to_name_map[vleft.ac]

        typ = vleft.posedit.edit.type

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
        elif typ == "inv":
            ref = vleft.posedit.edit.ref
            alt = reverse_complement(ref)
        else:
            alt = vleft.posedit.edit.alt or ""

            if typ in ("del", "ins"):
                if typ == "ins":
                    # ins coordinates (only) exclude left position
                    start_i += 1
                    end_i -= 1
                # Left anchored
                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = "."
        return chrom, start_i + 1, ref, alt, typ

    def vcf_to_g_hgvs(self, chrom, position, ref, alt):
        # VCF spec https://samtools.github.io/hts-specs/VCFv4.1.pdf
        # says for REF/ALT "Each base must be one of A,C,G,T,N (case insensitive)"
        ref = ref.upper()
        alt = alt.upper()

        ac = self.name_to_ac_map[chrom]

        if ref != alt:
            # Strip common prefix
            if len(alt) > 1 and len(ref) > 1:
                pfx = os.path.commonprefix([ref, alt])
                lp = len(pfx)
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

        var_g = SequenceVariant(
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
        n = Normalizer(self.hdp)
        return n.normalize(var_g)
