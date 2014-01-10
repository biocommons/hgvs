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

    def __init__(self, bdi, ac, ref='GRCH37.p10'):
        self.bdi = bdi
        self.ref = ref
        self.ac = ac
        self.tx_info = bdi.get_tx_info(self.ac)
        self.tx_exons = bdi.get_tx_exons(self.ac, ref)
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
        assert self.strand in [1,-1], 'strand = '+str(self.strand)+'; must be 1 or -1'

        g_ci = hgvs_coord_to_ci(g_interval.start.base, g_interval.end.base)
        # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
        frs, fre = self.im.map_ref_to_tgt(g_ci[0] - self.gc_offset, g_ci[1] - self.gc_offset, max_extent=False)
        if self.strand == -1:
            frs, fre = self.im.tgt_len - fre, self.im.tgt_len - frs
        # get BaseOffsetPosition information for r.
        if self.strand == 1:
            grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
        elif self.strand == -1:
            grs, gre = self.im.map_tgt_to_ref((self.im.tgt_len - fre), (self.im.tgt_len - frs), max_extent=False)
        grs, gre = grs + self.gc_offset, gre + self.gc_offset
        # 0-width interval indicates an intron.  Need to calculate offsets but we're are in ci coordinates
        # requires adding 1 strategically to get the HGVS position (shift coordinates to 3' end of the ref nucleotide)
        if frs == fre:
            start_offset = hgvs_offset(g_ci[0] + 1, grs, gre + 1)
            end_offset = hgvs_offset(g_ci[1], grs, gre + 1)
        else:
            start_offset = hgvs_offset(g_ci[0], grs, gre)
            end_offset = hgvs_offset(g_ci[1], grs, gre)
        if self.strand == -1:
            start_offset, end_offset = self.strand * end_offset, self.strand * start_offset
        if start_offset > 0:
            frs -= 1
        if end_offset < 0:
            fre += 1
        return hgvs.location.Interval(
                    start=hgvs.location.BaseOffsetPosition(base=ci_to_hgvs_coord(frs, fre)[0], offset=start_offset),
                    end  =hgvs.location.BaseOffsetPosition(base=ci_to_hgvs_coord(frs, fre)[1], offset=end_offset),
                    uncertain=g_interval.uncertain)

    def hgvsr_to_hgvsg(self, r_interval):
        assert self.strand in [1,-1], 'strand = '+str(self.strand)+'; must be 1 or -1'

        if self.strand == 1:
            frs, fre = hgvs_coord_to_ci(r_interval.start.base, r_interval.end.base)
            start_offset, end_offset = r_interval.start.offset, r_interval.end.offset
        elif self.strand == -1:
            frs, fre = hgvs_coord_to_ci(r_interval.start.base, r_interval.end.base)
            fre, frs = self.im.tgt_len - frs, self.im.tgt_len - fre
            start_offset, end_offset = self.strand * r_interval.end.offset, self.strand * r_interval.start.offset

        # returns the genomic range start (grs) and end (gre)
        grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
        grs, gre = grs + self.gc_offset, gre + self.gc_offset
        gs, ge = grs + start_offset, gre + end_offset
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(ci_to_hgvs_coord(gs, ge)[0], uncertain=r_interval.start.uncertain),
            end  =hgvs.location.SimplePosition(ci_to_hgvs_coord(gs, ge)[1], uncertain=r_interval.end.uncertain),
            uncertain=r_interval.uncertain
            )

    def hgvsr_to_hgvsc(self, r_interval):
        # start
        if r_interval.start.base <= self.cds_start_i:
            cs = r_interval.start.base - (self.cds_start_i + 1)
            cs_datum = hgvs.location.CDS_START
        elif r_interval.start.base > self.cds_start_i and r_interval.start.base <= self.cds_end_i:
            cs = r_interval.start.base - self.cds_start_i
            cs_datum = hgvs.location.CDS_START
        else:
            cs = r_interval.start.base - self.cds_end_i
            cs_datum = hgvs.location.CDS_END
        # end
        if r_interval.end.base <= self.cds_start_i:
            ce = r_interval.end.base - (self.cds_start_i + 1)
            ce_datum = hgvs.location.CDS_START
        elif r_interval.end.base > self.cds_start_i and r_interval.end.base <= self.cds_end_i:
            ce = r_interval.end.base - self.cds_start_i
            ce_datum = hgvs.location.CDS_START
        else:
            ce = r_interval.end.base - self.cds_end_i
            ce_datum = hgvs.location.CDS_END

        c_interval = hgvs.location.Interval(
            start=hgvs.location.BaseOffsetPosition(base=cs, offset=r_interval.start.offset, datum=cs_datum),
            end  =hgvs.location.BaseOffsetPosition(base=ce, offset=r_interval.end.offset,   datum=ce_datum),
            uncertain=r_interval.uncertain
            )
        return c_interval

    def hgvsc_to_hgvsr(self, c_interval):
        # start
        if c_interval.start.datum == hgvs.location.CDS_START and c_interval.start.base < 0:
            rs = c_interval.start.base + self.cds_start_i + 1
        elif c_interval.start.datum == hgvs.location.CDS_START and c_interval.start.base > 0:
            rs = c_interval.start.base + self.cds_start_i
        elif c_interval.start.datum == hgvs.location.CDS_END:
            rs = c_interval.start.base + self.cds_end_i
        # end
        if c_interval.end.datum == hgvs.location.CDS_START and c_interval.end.base < 0:
            re = c_interval.end.base + self.cds_start_i + 1
        elif c_interval.end.datum == hgvs.location.CDS_START and c_interval.end.base > 0:
            re = c_interval.end.base + self.cds_start_i
        elif c_interval.end.datum == hgvs.location.CDS_END:
            re = c_interval.end.base + self.cds_end_i

        r_interval = hgvs.location.Interval(
            start=hgvs.location.BaseOffsetPosition(base=rs, offset=c_interval.start.offset, datum=hgvs.location.SEQ_START),
            end  =hgvs.location.BaseOffsetPosition(base=re, offset=c_interval.end.offset, datum=hgvs.location.SEQ_START),
            uncertain=c_interval.uncertain
            )
        return r_interval

    def hgvsg_to_hgvsc(self, g_interval):
        return self.hgvsr_to_hgvsc(self.hgvsg_to_hgvsr(g_interval))

    def hgvsc_to_hgvsg(self, c_interval):
        return self.hgvsr_to_hgvsg(self.hgvsc_to_hgvsr(c_interval))


def hgvs_offset(g_position, grs, gre):
    """ Calculates the HGVS coordinate offset from a given genomic position"""
    if g_position == grs or g_position == gre:
        return 0
    mid = math.ceil((grs + gre) / 2)
    if g_position < mid:
        offset = g_position - grs
    else:
        offset = g_position - gre
    return offset


def ci_to_hgvs_coord(s, e):
    """ Convert continuous interbase (right-open) coordinates (..,-2,-1,0,1,..) to
    discontinuous HGVS coordinates (..,-2,-1,1,2,..)
    """
    def _ci_to_hgvs(c):
        return c + 1 if c >= 0 else c
    return (None if s is None else _ci_to_hgvs(s),
            None if e is None else _ci_to_hgvs(e) - 1)

def hgvs_coord_to_ci(s, e):
    """convert start,end interval in inclusive, discontinuous HGVS coordinates
    (..,-2,-1,1,2,..) to continuous interbase (right-open) coordinates
    (..,-2,-1,0,1,..)"""
    def _hgvs_to_ci(c):
        assert c != 0, 'received CDS coordinate 0; expected ..,-2,-1,1,1,...'
        return c-1 if c>0 else c
    return (None if s is None else _hgvs_to_ci(s),
            None if e is None else _hgvs_to_ci(e) + 1)


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
    import bdi.sources.uta0
    ref = 'GRCh37.p10'
    ac = 'NM_145249.2'
    bdi = bdi.sources.uta0.UTA0('/tmp/uta-0.0.4.db')
    tm = TranscriptMapper(bdi,ac,ref)

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
