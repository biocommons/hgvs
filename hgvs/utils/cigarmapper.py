import itertools
import re



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



    # i'base   0 1 2 3 4    5    6 7 8 9 0 1 2 3 4  
    # base      0 1 2 3 4       5 6 7 8 9 0 1 2 3    
    # ref       A B C D E       I J K L M N O P Q    len=14
    #          |= = = = =|I I I|= = = =|D D|= = =|   
    #          |    5=   | 3I  |   4=  | 2D| 3=  |   5= 3I 4= 2D 3=
    # tgt       A B C D E F G H I J K L     O P Q    len=15
    # base      0 1 2 3 4 5 6 7 8 9 0 1     2 3 4  
    # i'base   0 1 2 3 4 5 6 7 8 9 0 1   2   3 4 5 

    """
    
    def __init__(self, cigar):
        self.cigar = cigar
        self._process_cigar()


    def _process_cigar(self):
        """build data structures used for mapping
        """

        cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")
        advance_ref = "=MXD"
        advance_tgt = "=MXIN"

        cigar_elems = list(m.groupdict() for m in cigar_re.finditer(self.cigar)) 
        self.ref_lengths = [int(ce["len"]) if ce["op"] in advance_ref else 0 for ce in cigar_elems] 
        self.tgt_lengths = [int(ce["len"]) if ce["op"] in advance_tgt else 0 for ce in cigar_elems] 
        self.ref_pos = [0] + list(itertools.accumulate(self.ref_lengths))
        self.tgt_pos = [0] + list(itertools.accumulate(self.tgt_lengths))
        self.ref_length = self.ref_pos[-1]
        self.tgt_length = self.tgt_pos[-1]
        self.cigar_ops = [ce["op"] for ce in cigar_elems]



    def _map(self, from_pos, to_pos, pos, base):
        """Map position between aligned sequences

        Positions in this function are 0-based.
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
            if base == "start":
                mapped_pos = to_pos[pos_i] - 1
            elif base == "end":
                mapped_pos = to_pos[pos_i]
            mapped_pos_offset = 0
        elif self.cigar_op[pos_i] == "N":
            if pos - from_pos[pos_i] + 1 <= from_pos[pos_i + 1] - pos:
                mapped_pos = to_pos[pos_i] - 1
                mapped_pos_offset = pos - from_pos[pos_i] + 1
            else:
                mapped_pos = to_pos[pos_i]
                mapped_pos_offset = -(from_pos[pos_i + 1] - pos)

        return mapped_pos, mapped_pos_offset, self.cigar_op[pos_i]

    





cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")

class AlignmentMapper:
    def __init__(self, cigar):
            self.cigar = cigar
            self.ref_pos, self.tgt_pos, self.cigar_op = self._parse_cigar(self.cigar)
            self.tgt_len = self.tgt_pos[-1]

    def _parse_cigar(self, cigar):
        """For a given CIGAR string, return the start positions of
        each aligned segment in ref and tgt, and a list of CIGAR operators.
        """
        ces = [m.groupdict() for m in cigar_re.finditer(cigar)]
        ref_pos = [None] * len(ces)
        tgt_pos = [None] * len(ces)
        cigar_op = [None] * len(ces)
        ref_cur = tgt_cur = 0
        for i, ce in enumerate(ces):
            ref_pos[i] = ref_cur
            tgt_pos[i] = tgt_cur
            cigar_op[i] = ce["op"]
            step = int(ce["len"])
            if ce["op"] in "=MINX":
                ref_cur += step
            if ce["op"] in "=MDX":
                tgt_cur += step
        ref_pos.append(ref_cur)
        tgt_pos.append(tgt_cur)
        return ref_pos, tgt_pos, cigar_op

    def _map(self, from_pos, to_pos, pos, base):
        """Map position between aligned sequences

        Positions in this function are 0-based.
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
            if base == "start":
                mapped_pos = to_pos[pos_i] - 1
            elif base == "end":
                mapped_pos = to_pos[pos_i]
            mapped_pos_offset = 0
        elif self.cigar_op[pos_i] == "N":
            if pos - from_pos[pos_i] + 1 <= from_pos[pos_i + 1] - pos:
                mapped_pos = to_pos[pos_i] - 1
                mapped_pos_offset = pos - from_pos[pos_i] + 1
            else:
                mapped_pos = to_pos[pos_i]
                mapped_pos_offset = -(from_pos[pos_i + 1] - pos)

        return mapped_pos, mapped_pos_offset, self.cigar_op[pos_i]




    

if __name__ == "__main__":

    am = AlignmentMapper("5= 3I 4= 2D 3=")
    cm = CIGARMapper("5= 3I 4= 2D 3=")
    
