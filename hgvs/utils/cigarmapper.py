from hgvs.exceptions import HGVSInvalidIntervalError

from hgvs.utils import parse_cigar


class CIGARMapper:
    """provides coordinate mapping between two sequences whose alignment
    is given by a CIGAR string
    
    CIGAR is about alignments between bases.  It is base-centric.
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
        self.ref_pos, self.tgt_pos, self.cigar_op = parse_cigar(self.cigar)


    def __repr__(self):
        return f"CIGARMapper(cigar='{self.cigar}')"


    @property
    def ref_len(self):
        return self.ref_pos[-1]


    @property
    def tgt_len(self):
        return self.tgt_pos[-1]


    def map_ref_to_tgt(self, pos, end):
        return self._map(from_pos=self.ref_pos, to_pos=self.tgt_pos, pos=pos, end=end)


    def map_tgt_to_ref(self, pos, end):
        return self._map(from_pos=self.tgt_pos, to_pos=self.ref_pos, pos=pos, end=end)


    def _map(self, from_pos, to_pos, pos, end):
        """Map position between aligned segments

        Positions in this function are 0-based, base-counting. 
        """
        pos_i = -1
        while pos_i < len(self.cigar_op) and pos >= from_pos[pos_i + 1]:
            pos_i += 1

        if pos_i == -1 or pos_i == len(self.cigar_op):
            raise HGVSInvalidIntervalError("Position is beyond the bounds of transcript record")
        if self.cigar_op[pos_i] in "=MX":
            mapped_pos = to_pos[pos_i] + (pos - from_pos[pos_i])
            mapped_pos_offset = 0
        elif self.cigar_op[pos_i] in "DI":
            mapped_pos = to_pos[pos_i]
            if base == "start":
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
