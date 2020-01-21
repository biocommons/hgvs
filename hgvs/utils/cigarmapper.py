import itertools
import re

from hgvs.exceptions import HGVSInvalidIntervalError


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
        self._process_cigar()


    def _process_cigar(self):
        """build data structures used for mapping
        """

        cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")
        advance_ref = "=MXD"
        advance_tgt = "=MXIN"

        cigar_elems = list(m.groupdict() for m in cigar_re.finditer(self.cigar)) 
        self.cigar_ops = [ce["op"] for ce in cigar_elems]

        ref_lens = [int(ce["len"]) if ce["op"] in advance_ref else 0 for ce in cigar_elems] 
        tgt_lens = [int(ce["len"]) if ce["op"] in advance_tgt else 0 for ce in cigar_elems] 
        self.ref_pos = [0] + list(itertools.accumulate(ref_lens))
        self.tgt_pos = [0] + list(itertools.accumulate(tgt_lens))

        self.ref_len = self.ref_pos[-1] + 1
        self.tgt_len = self.tgt_pos[-1] + 1


    def _map(self, from_pos, to_pos, pos, base):
        """Map position between aligned segments

        Positions in this function are 0-based, base-counting. 
        """
        pos_i = -1
        while pos_i < len(self.cigar_ops) and pos >= to_pos[pos_i + 1]:
            pos_i += 1

        if pos_i == -1 or pos_i == len(self.cigar_ops):
            raise HGVSInvalidIntervalError("Position is beyond the bounds of transcript record")

        mapped_pos_offset = 0
        if self.cigar_ops[pos_i] in "=MX":
            mapped_pos = from_pos[pos_i] + (pos - to_pos[pos_i])
        elif self.cigar_ops[pos_i] in "DI":
            if base == "start":
                mapped_pos = from_pos[pos_i] - 1
            elif base == "end":
                mapped_pos = from_pos[pos_i]
        elif self.cigar_ops[pos_i] == "N":
            if pos - to_pos[pos_i] + 1 <= to_pos[pos_i + 1] - pos:
                mapped_pos = from_pos[pos_i] - 1
                mapped_pos_offset = pos - to_pos[pos_i] + 1
            else:
                mapped_pos = from_pos[pos_i]
                mapped_pos_offset = -(to_pos[pos_i + 1] - pos)

        return mapped_pos, mapped_pos_offset, self.cigar_ops[pos_i]

    def _map_dir(self, dir, pos, base):
        if dir == "rt":
            return self._map(self.ref_pos, self.tgt_pos, pos, base)
        if dir == "tr":
            return self._map(self.tgt_pos, self.ref_pos, pos, base)
        assert False, "Invalid `dir` argument"

    def map_ref_to_tgt(self, pos, base):
        return self._map_dir("rt", pos, base)

    def map_tgt_to_ref(self, pos, base):
        return self._map_dir("tr", pos, base)



class _AlignmentMapper:
    """this was ripped out of the AlignmentMapper for comparison purposes"""
    def __init__(self, cigar):
        self.cigar = cigar
        self.tgt_pos, self.ref_pos, self.cigar_op = self._parse_cigar(self.cigar)
        self.ref_len = self.ref_pos[-1]
        self.tgt_len = self.tgt_pos[-1]

    def _parse_cigar(self, cigar):
        """For a given CIGAR string, return the start positions of
        each aligned segment in ref and tgt, and a list of CIGAR operators.
        """
        cigar_re = re.compile(r"(?P<len>\d+)(?P<op>[=DIMNX])")
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
            # These two cases are swapped!
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


    def ref_to_target_interval(self, start, end):
        """projects a (start,end) interbase interval on the reference sequence
        onto the target sequence, with offset a start_offset 
        """
        pass
    




if __name__ == "__main__":

    # ref (tx)     A B C          F G H   J K L    M N    O P Q R S T U V 
    # tgt          A B C   D E    F G H I J K L    M N    O P Q     T U V 


    # i'base      0 1 2   3   4 5  6  7 8     9     0 1 2 3 4 5
    # base         0 1 2     3 4 5   6 7 8         9 0 1 2 3 4
    # ref (tx)     A B C     F G H   J K L         Q R S T U V 
    #              | | | - - | | | - | | |         | - - | | |
    #              = = = N N = = = I = = = N N N N = D D = X = 
    #                3=   2N    3=1I3=       4N      1=2D1X1=
    # tgt          A B C D E F G H I J K L M N O P Q - - T Z V 
    # base         0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6     7 8 9
    # i'base      0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6   7   8 9 0
    # ref len is 15 bases, target is 20
    # 

    cigar =  "3= 2N 3=1I3= 4N 1=2D1X1="

    am = _AlignmentMapper(cigar)
    print(str(am.ref_len) + " // " + str(am.ref_pos))
    print(str(am.tgt_len) + " // " + str(am.tgt_pos))

    cm = CIGARMapper(cigar)
    print(str(cm.ref_len) + " // " + str(cm.ref_pos))
    print(str(cm.tgt_len) + " // " + str(cm.tgt_pos))


    for refpos in range(cm.ref_len + 2):
        ams = am._map(am.tgt_pos, am.ref_pos, refpos, "start")
        ame = am._map(am.tgt_pos, am.ref_pos, refpos, "end")
        #cms = cm._map(cm.tgt_pos, cm.ref_pos, refpos, "start")
        #cme = cm._map(cm.tgt_pos, cm.ref_pos, refpos, "end")
        cms = cm.map_ref_to_tgt(refpos, "start")
        cme = cm.map_ref_to_tgt(refpos, "end")
        seck = "✔" if ams == ame else " "
        sck = "✔" if ams == cms else " "
        eck = "✔" if ame == cme else " "
        print(f"{refpos}:  ({cms},{cme})  ({ams},{ame})  ({seck},{sck},{eck})") 
        #print(f"{refpos}:  ({cms[0]} {cms[1]}  ({cme[0]} {cme[1]})   ({sck},{eck})") 
