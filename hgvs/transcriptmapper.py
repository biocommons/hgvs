"""
Maps transcription coordinates between g <-> r, c, p
"""
import hgvs.location
import hgvs.variant
import hgvs.posedit
import math
import re
from hgvs.intervalmapper import IntervalMapper
from hgvs.exceptions import *


class TranscriptMapper(object):
    __doc__ = """All coordinates are in HGVS format (human 1-based)"""

    def __init__(self, db, ac, ref='GRCH37.p10'):
        self.db = db
        self.ref = ref
        self.ac = ac
        self.tx_info = db.get_tx_info(self.ac)
        self.tx_exons = db.get_tx_exons(self.ac, ref)
        if self.tx_info is None or len(self.tx_exons) == 0:
            raise HGVSError("Couldn't build TranscriptMapper(ref={self.ref}, ac={self.ac})".format(self=self))
        self.strand = self.tx_info['strand']
        self.cds_start_i = self.tx_info['cds_start_i']
        self.cds_end_i = self.tx_info['cds_end_i']
        self.gc_offset = self.tx_exons[0]['g_start_i']
        self.cigar = build_tx_cigar(self.tx_exons, self.strand)
        self.im = IntervalMapper.from_cigar(self.cigar)

    def __str__(self):
        return '{self.__class__.__name__}: {self.ac} ~ {self.ref}; {self.strand_pm} strand; {n_exons} exons; ' \
               'offset={self.gc_offset}'.format(self=self, n_exons=len(self.tx_exons))

    @property
    def strand_pm(self):
        return (None if self.strand is None
                else '+' if self.strand == 1
                else '-' if self.strand == -1
                else '?')

    def hgvsg_to_hgvsr(self, g_interval):
        """ Convert hgvsg interval into an hgvsr interval
        Input: g_interval
        Output: r_interval
        """
        # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
        frs, fre = self.im.map_ref_to_tgt(g_interval.start - self.gc_offset, g_interval.end - self.gc_offset, max_extent=False)
        if self.strand == 1:
            pass
        elif self.strand == -1:
            frs, fre = self.im.tgt_len - fre, self.im.tgt_len - frs
        else:
            raise HGVSError('Code fell through strand check (g_to_r); should never get here.')
        # frs & fre will be equal except for cases with an insert in the transcript relative to the genome
        if frs == fre:
            # get BaseOffsetPosition information for r.
            if self.strand == -1:
                grs, gre = self.im.map_tgt_to_ref((self.im.tgt_len - fre), (self.im.tgt_len - frs), max_extent=False)
            elif self.strand == 1:
                grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
            grs, gre = grs + self.gc_offset, gre + self.gc_offset
            start_offset = hgvs_offset(g_interval.start, grs, gre)
            end_offset = hgvs_offset(g_interval.end, grs, gre)
        else:
            start_offset = 0
            end_offset = 0
        if self.strand == -1:
            start_offset = self.strand * start_offset
            end_offset = self.strand * end_offset
        return hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=ci_to_hgvs_coord(frs), offset=start_offset),
                    end=hgvs.location.BaseOffsetPosition(base=ci_to_hgvs_coord(fre), offset=end_offset))

    def hgvsr_to_hgvsg(self, r_interval):
        if self.strand == 1:
            frs = hgvs_coord_to_ci(r_interval.start.base)
            fre = hgvs_coord_to_ci(r_interval.end.base)
        elif self.strand == -1:
            frs = self.im.tgt_len - hgvs_coord_to_ci(r_interval.end.base)
            fre = self.im.tgt_len - hgvs_coord_to_ci(r_interval.start.base)
        else:
           raise HGVSError('Code fell through strand check (r_to_g); should never get here.')
        # returns the genomic range start (grs) and end (gre)
        start_offset, end_offset = self.strand * r_interval.start.offset, self.strand * r_interval.end.offset
        grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
        grs, gre = grs + self.gc_offset, gre + self.gc_offset
        gs = grs + start_offset if start_offset >= 0 else gre + start_offset
        ge = grs + end_offset if end_offset > 0 else gre + end_offset   # note > 0 accounts for 0-based offsets
        return hgvs.location.Interval(start=gs, end=ge)

    def hgvsr_to_hgvsc(self, r_interval):
        c_interval = hgvs.location.Interval(
                        start=hgvs.location.BaseOffsetPosition(base=r_interval.start.base - self.cds_start_i,
                                                               offset=r_interval.start.offset,
                                                               datum=hgvs.location.CDS_START),
                        end=hgvs.location.BaseOffsetPosition(base=r_interval.end.base - self.cds_start_i,
                                                             offset=r_interval.end.offset,
                                                             datum=hgvs.location.CDS_START))
        return c_interval

    def hgvsc_to_hgvsr(self, c_interval):
        r_interval = hgvs.location.Interval(
                        start=hgvs.location.BaseOffsetPosition(base=c_interval.start.base + self.cds_start_i,
                                                               offset=c_interval.start.offset,
                                                               datum=hgvs.location.SEQ_START),
                        end=hgvs.location.BaseOffsetPosition(base=c_interval.end.base + self.cds_start_i,
                                                             offset=c_interval.end.offset,
                                                             datum=hgvs.location.SEQ_START))
        return r_interval

    def hgvsg_to_hgvsc(self, g_interval):
        return self.hgvsr_to_hgvsc(self.hgvsg_to_hgvsr(g_interval))

    def hgvsc_to_hgvsg(self, c_interval):
        return self.hgvsr_to_hgvsg(self.hgvsc_to_hgvsr(c_interval))


def hgvs_offset(g_position, grs, gre):
    """ Calculates the hgvs coordinate offset from a given genomic position
    """
    mid = math.ceil((grs + gre) / 2)
    if g_position < mid:
        offset = g_position - grs
    else:
        offset = g_position - gre
    return offset


def ci_to_hgvs_coord(pos):
    """ Convert continuous interbase (right-open) coordinates (..,-2,-1,0,1,..) to
    discontinuous HGVS coordinates (..,-2,-1,1,2,..)
    """
    return pos + 1 if pos >= 0 else pos


def hgvs_coord_to_ci(pos):
    """convert start,end interval in inclusive, discontinuous HGVS coordinates
    (..,-2,-1,1,2,..) to continuous interbase (right-open) coordinates
    (..,-2,-1,0,1,..)"""
    def _cds_to_ci(c):
        assert c != 0, 'received CDS coordinate 0; expected ..,-2,-1,1,1,...'
        return c-1 if c>0 else c
    return _cds_to_ci(pos)


def build_tx_cigar(exons, strand):
    if len(exons) == 0:
        return None

    cigarelem_re = re.compile('\d+[DIMNX]')

    def _reverse_cigar(c):
        return ''.join(reversed(cigarelem_re.findall(c)))

    if strand == -1:
        for i in range(len(exons)):
            exons[i]['g_cigar'] = _reverse_cigar(exons[i]['g_cigar'])

    tx_cigar = [exons[0]['g_cigar']]  # exon 1
    for i in range(1, len(exons)):     # and intron + exon pairs thereafter
        tx_cigar += [str(exons[i]['g_start_i'] - exons[i - 1]['g_end_i']) + 'N',
                     exons[i]['g_cigar']]
    return ''.join(tx_cigar)


if __name__ == '__main__':
    from uta.db.transcriptdb import TranscriptDB

    ref = 'GRCh37.p10'
    ac = 'NM_145249.2'
    db = TranscriptDB()
    tm = TranscriptMapper(db,ac,ref)