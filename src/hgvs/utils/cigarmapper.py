"""provides CIGAR-based mapping for hgvs

Heads up: UTA CIGAR strings are relative to the transcript; that is,
the transcript is the reference and the genome is the target sequence.
hgvs, and this module specifically, swaps this nomenclature so that
the genome is ref and the transcript is the target.

UTA CIGAR strings are relative to the transcript.  Specifically:

  * D is a deletion in the genome relative to the transcript
  * I is an insertion in the genome relative to the transcript
  * N is an intronic insertion relative to the transcript

The conceptual swap occurs in the _map function.

"""


import re

from hgvs.exceptions import HGVSInvalidIntervalError


class CIGARMapper:
    """provides coordinate mapping between two sequences whose alignment
    is given by a CIGAR string

    CIGAR is about alignments between positions in two sequences.  It
    is base-centric.

    Unfortunately, base-centric coordinate systems require additional
    complexity to refer to zero-width positions.

    This code uses interbase intervals.  Interbase positions are
    zero-width boundaries between bases.  They often look similar to
    zero-based, right open coordinates. (But don't call them that.  It
    upsets me deeply.)  The most important difference is that zero
    width intervals neatly represent insertions between bases (or
    before or after the sequence).

    """

    def __init__(self, cigar):
        self.cigar = cigar
        self.ref_pos, self.tgt_pos, self.cigar_op = _parse_cigar(self.cigar)

    def __repr__(self):
        return f"CIGARMapper(cigar='{self.cigar}')"

    @property
    def ref_len(self):
        return self.ref_pos[-1]

    @property
    def tgt_len(self):
        return self.tgt_pos[-1]

    def map_ref_to_tgt(self, pos, end, strict_bounds=True):
        return self._map(from_pos=self.ref_pos, to_pos=self.tgt_pos, pos=pos, end=end, strict_bounds=strict_bounds)

    def map_tgt_to_ref(self, pos, end, strict_bounds=True):
        return self._map(from_pos=self.tgt_pos, to_pos=self.ref_pos, pos=pos, end=end, strict_bounds=strict_bounds)

    def _map(self, from_pos, to_pos, pos, end, strict_bounds):
        """Map position between aligned segments

        Positions in this function are 0-based, base-counting.
        """

        if strict_bounds and (pos < 0 or pos > from_pos[-1]):
            raise HGVSInvalidIntervalError("Position is beyond the bounds of transcript record")

        # find aligned segment to use as basis for mapping
        # okay for pos to be before first element or after last
        for pos_i in range(len(self.cigar_op)):
            if pos < from_pos[pos_i + 1]:
                break

        if self.cigar_op[pos_i] in "=MX":
            mapped_pos = to_pos[pos_i] + (pos - from_pos[pos_i])
            mapped_pos_offset = 0

        elif self.cigar_op[pos_i] in "DI":
            mapped_pos = to_pos[pos_i]
            if end == "start":
                mapped_pos -= 1
            mapped_pos_offset = 0

        elif self.cigar_op[pos_i] == "N":
            if pos - from_pos[pos_i] + 1 <= from_pos[pos_i + 1] - pos:
                mapped_pos = to_pos[pos_i] - 1
                mapped_pos_offset = pos - from_pos[pos_i] + 1
            else:
                mapped_pos = to_pos[pos_i]
                mapped_pos_offset = -(from_pos[pos_i + 1] - pos)

        return mapped_pos, mapped_pos_offset, self.cigar_op[pos_i]


def _parse_cigar(cigar):
    """For a given CIGAR string, return the start positions of each
    aligned segment in ref and tgt, and a list of CIGAR operators.

    """

    cigar_re = re.compile(r"(?P<len>\d+)?(?P<op>[=DIMNX])")
    advance_ref_ops = "=MXIN"
    advance_tgt_ops = "=MXD"

    ces = [m.groupdict() for m in cigar_re.finditer(cigar)]
    cigar_len = len(ces)

    ref_pos = [None] * cigar_len
    tgt_pos = [None] * cigar_len
    cigar_op = [None] * cigar_len
    ref_cur = tgt_cur = 0
    for i, ce in enumerate(ces):
        ref_pos[i] = ref_cur
        tgt_pos[i] = tgt_cur
        cigar_op[i] = ce["op"]
        step = int(ce["len"] or 1)
        # Note: these cases are
        if ce["op"] in advance_ref_ops:
            ref_cur += step
        if ce["op"] in advance_tgt_ops:
            tgt_cur += step
    ref_pos.append(ref_cur)
    tgt_pos.append(tgt_cur)
    return ref_pos, tgt_pos, cigar_op


if __name__ == "__main__":
    # cigar = "2=2N=X=2N=I=2N=D="
    cigar = "3=2N=X=3N=I=D="
    cm = CIGARMapper(cigar)
