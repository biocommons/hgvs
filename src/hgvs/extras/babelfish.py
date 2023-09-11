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


def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        # ins coordinates (only) exclude left position
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)


class Babelfish:
    def __init__(self, hdp, assembly_name):
        self.assembly_name = assembly_name
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(hdp, cross_boundaries=False, shuffle_direction=5, validate=False)
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)

    def hgvs_to_vcf(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, pos, ref, alt) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).
        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))

        vleft = self.hn.normalize(var_g)

        (start_i, end_i) = _as_interbase(vleft.posedit)

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

            if typ in ("del", "ins"):  # Left anchored
                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = "."
        return chrom, start_i + 1, ref, alt, typ

    def vcf_to_g_hgvs(self, chrom, position, ref, alt):
        ac = self.name_to_ac_map[chrom]

        # Strip common prefix
        if len(alt) > 1 and len(ref) > 1:
            pfx = os.path.commonprefix([ref, alt])
            lp = len(pfx)
            if lp > 0:
                ref = ref[lp:]
                alt = alt[lp:]
                position += lp

        if ref == "":  # Insert
            # Insert uses coordinates around the insert point.
            start = position - 1
            end = position
        else:
            start = position
            end = position + len(ref) - 1

        if alt == ".":
            alt = ref

        var_g = SequenceVariant(
            ac=ac,
            type="g",
            posedit=PosEdit(
                Interval(start=SimplePosition(start), end=SimplePosition(end), uncertain=False),
                NARefAlt(ref=ref or None, alt=alt or None, uncertain=False),
            ),
        )
        n = Normalizer(self.hdp)
        return n.normalize(var_g)


if __name__ == "__main__":
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """

    import hgvs.easy
    import hgvs.normalizer

    babelfish38 = Babelfish(hgvs.easy.hdp, assembly_name="GRCh38")
    hnl = hgvs.normalizer.Normalizer(hgvs.easy.hdp, cross_boundaries=False, shuffle_direction=5, validate=False)

    def _h2v(h):
        return babelfish38.hgvs_to_vcf(hgvs.easy.parser.parse(h))

    def _vp(h):
        v = hgvs.easy.parser.parse(h)
        vl = hnl.normalize(v)
        return (v, vl)

    for h in (
        # Non-variation
        "NC_000006.12:g.49949407=",
        # SNV
        "NC_000006.12:g.49949407A>T",
        # delins
        "NC_000006.12:g.49949413_49949414delinsCC",
        # del
        "NC_000006.12:g.49949415del",
        "NC_000006.12:g.49949413del",
        "NC_000006.12:g.49949414del",
        "NC_000006.12:g.49949413_49949414del",
        # ins
        "NC_000006.12:g.49949413_49949414insC",
        "NC_000006.12:g.49949414_49949415insC",
        "NC_000006.12:g.49949414_49949415insCC",
        # ins (dup)
        "NC_000006.12:g.49949413_49949414insA",
        "NC_000006.12:g.49949414_49949415insA",
        "NC_000006.12:g.49949414_49949415insAA",
    ):
        print('assert _h2v("{h}") == {res}'.format(res=str(_h2v(h)), h=h))
