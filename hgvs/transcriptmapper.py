# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from bioutils.coordinates import strand_int_to_pm

import hgvs.intervalmapper
import hgvs.location
import hgvs.posedit
import hgvs.variant

from hgvs.exceptions import HGVSError
from hgvs.decorators.deprecated import deprecated
from hgvs.utils import build_tx_cigar


class TranscriptMapper(object):
    """Provides coordinate (not variant) mapping operations between
    genomic (g), rna (r), cds (c), and protein (p) coordinates.  All
    coordinates are 1-based inclusive, per the HGVS recommendations.
    All methods take :class:`hgvs.location.Interval` objects.

    :param hdp: HGVS Data Provider Interface-compliant instance (see :class:`hgvs.dataproviders.interface.Interface`)
    :param str tx_ac: string representing transcript accession (e.g., NM_000551.2)
    :param str alt_ac: string representing the reference sequence accession (e.g., NM_000551.3)
    :param str alt_aln_method: string representing the alignment method; valid values depend on data source

    """

    def __init__(self, hdp, tx_ac, alt_ac, alt_aln_method):
        self.hdp = hdp
        self.tx_ac = tx_ac
        self.alt_ac = alt_ac
        self.alt_aln_method = alt_aln_method
        if self.alt_aln_method != 'transcript':
            self.tx_info = hdp.get_tx_info(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if self.tx_info is None:
                raise HGVSError("TranscriptMapper(tx_ac={self.tx_ac}, "
                                "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                                "No transcript info".format(self=self))

            self.tx_exons = hdp.get_tx_exons(self.tx_ac, self.alt_ac, self.alt_aln_method)
            if self.tx_exons is None:
                raise HGVSError("TranscriptMapper(tx_ac={self.tx_ac}, "
                                "alt_ac={self.alt_ac}, alt_aln_method={self.alt_aln_method}): "
                                "No transcript exons".format(self=self))

            self.strand = self.tx_exons[0]['alt_strand']
            self.cds_start_i = self.tx_info['cds_start_i']
            self.cds_end_i = self.tx_info['cds_end_i']
            self.gc_offset = self.tx_exons[0]['alt_start_i']
            self.cigar = build_tx_cigar(self.tx_exons, self.strand)
            self.im = hgvs.intervalmapper.IntervalMapper.from_cigar(self.cigar)
            self.tgt_len = self.im.tgt_len
        else:
            # this covers the identity cases r <-> c
            self.tx_identity_info = hdp.get_tx_identity_info(self.tx_ac)
            self.cds_start_i = self.tx_identity_info['cds_start_i']
            self.cds_end_i = self.tx_identity_info['cds_end_i']
            self.tgt_len = sum(self.tx_identity_info['lengths'])


    def __str__(self):
        return '{self.__class__.__name__}: {self.tx_ac} ~ {self.alt_ac} ~ {self.alt_aln_method); ' \
               '{strand_pm} strand; {n_exons} exons; offset={self.gc_offset}'.format(
                   self=self, n_exons=len(self.tx_exons), strand_pm=strand_int_to_pm(self.strand))


    def g_to_r(self, g_interval):
        """convert a genomic (g.) interval to a transcript cDNA (r.) interval"""

        # This code is extremely convoluted. To begin with, it
        # confuses interbase intervals for a *single* hgvs position
        # with intervals of *two* hgvs positions.  This is the origin
        # of bug #205 (now closed). Although the bug is itself fixed,
        # the mappers really need an overhaul. See #208.

        def _hgvs_offset(g_position, grs, gre, strand):
            """ Calculates the HGVS coordinate offset from a given genomic position"""
            if g_position == grs or g_position == gre:
                return 0
            mid = (grs + gre) / 2
            if g_position < mid or (g_position == mid and strand == 1):
                offset = g_position - grs
            else:
                offset = g_position - gre
            return offset

        def map_g_to_r_pos(pos):
            g_ci = _hgvs_coord_to_ci(pos,pos)
            # frs, fre = (f)orward (r)na (s)tart & (e)nd; forward w.r.t. genome
            frs, fre = self.im.map_ref_to_tgt(g_ci[0] - self.gc_offset, g_ci[1] - self.gc_offset, max_extent=False)
            if self.strand == -1:
                frs, fre = self.tgt_len - fre, self.tgt_len - frs
            # get BaseOffsetPosition information for r.
            if self.strand == 1:
                grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
            elif self.strand == -1:
                grs, gre = self.im.map_tgt_to_ref((self.tgt_len - fre), (self.tgt_len - frs), max_extent=False)
            grs, gre = grs + self.gc_offset, gre + self.gc_offset
            # 0-width interval indicates an intron.  Need to calculate offsets but we're are in ci coordinates
            # requires adding 1 strategically to get the HGVS position (shift coordinates to 3' end of the ref nucleotide)
            if frs == fre:
                start_offset = _hgvs_offset(g_ci[0] + 1, grs, gre + 1, self.strand)
                end_offset = _hgvs_offset(g_ci[1], grs, gre + 1, self.strand)
            else:
                start_offset = _hgvs_offset(g_ci[0], grs, gre, self.strand)
                end_offset = _hgvs_offset(g_ci[1], grs, gre, self.strand)
            if self.strand == -1:
                start_offset, end_offset = self.strand * end_offset, self.strand * start_offset
            if start_offset > 0:
                frs -= 1
            if end_offset < 0:
                fre += 1
            start_base = _ci_to_hgvs_coord(frs, fre)[0]
            end_base = _ci_to_hgvs_coord(frs, fre)[1]
            return start_base, start_offset, end_base, end_offset

        start_bo = map_g_to_r_pos(g_interval.start.base)[0:2]
        end_bo = start_bo if g_interval.start.base == g_interval.end.base else map_g_to_r_pos(g_interval.end.base)[2:4]
        if self.strand == -1:
            start_bo, end_bo = end_bo, start_bo

        return hgvs.location.Interval(
            start=hgvs.location.BaseOffsetPosition(base=start_bo[0], offset=start_bo[1]),
            end=hgvs.location.BaseOffsetPosition(base=end_bo[0], offset=end_bo[1]),
            uncertain=g_interval.uncertain)


    def r_to_g(self, r_interval):
        """convert a transcript cDNA (r.) interval to a genomic (g.) interval"""

        assert self.strand in [1,-1], 'strand = '+str(self.strand)+'; must be 1 or -1'

        if self.strand == 1:
            frs, fre = _hgvs_coord_to_ci(r_interval.start.base, r_interval.end.base)
            start_offset, end_offset = r_interval.start.offset, r_interval.end.offset
        elif self.strand == -1:
            frs, fre = _hgvs_coord_to_ci(r_interval.start.base, r_interval.end.base)
            fre, frs = self.tgt_len - frs, self.tgt_len - fre
            start_offset, end_offset = self.strand * r_interval.end.offset, self.strand * r_interval.start.offset

        # returns the genomic range start (grs) and end (gre)
        grs, gre = self.im.map_tgt_to_ref(frs, fre, max_extent=False)
        grs, gre = grs + self.gc_offset, gre + self.gc_offset
        gs, ge = grs + start_offset, gre + end_offset
        return hgvs.location.Interval(
            start=hgvs.location.SimplePosition(_ci_to_hgvs_coord(gs, ge)[0], uncertain=r_interval.start.uncertain),
            end  =hgvs.location.SimplePosition(_ci_to_hgvs_coord(gs, ge)[1], uncertain=r_interval.end.uncertain),
            uncertain=r_interval.uncertain)


    def r_to_c(self, r_interval):
        """convert a transcript cDNA (r.) interval to a transcript CDS (c.) interval"""

        if r_interval.start.base <= 0:
            raise HGVSError("Coordinate out of bounds. Start position ({rs}) is <= 0.".
                            format(rs=r_interval.start.base))
        if r_interval.end.base > self.tgt_len:
            raise HGVSError("Coordinate out of bounds. End position ({re}) is > than transcript length ({len}).".
                            format(re=r_interval.end.base, len=self.tgt_len))
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
            uncertain=r_interval.uncertain)
        return c_interval


    def c_to_r(self, c_interval):
        """convert a transcript CDS (c.) interval to a transcript cDNA (r.) interval"""

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

        if rs <= 0:
            raise HGVSError("Coordinate out of bounds. Start position ({rs}) is <= 0.".format(rs=rs))
        if re > self.tgt_len:
            raise HGVSError("Coordinate out of bounds. End position ({re}) is > than transcript length ({len}).".
                            format(re=re, len=self.tgt_len))

        r_interval = hgvs.location.Interval(
            start=hgvs.location.BaseOffsetPosition(base=rs, offset=c_interval.start.offset, datum=hgvs.location.SEQ_START),
            end  =hgvs.location.BaseOffsetPosition(base=re, offset=c_interval.end.offset, datum=hgvs.location.SEQ_START),
            uncertain=c_interval.uncertain)
        return r_interval

    def g_to_c(self, g_interval):
        """convert a genomic (g.) interval to a transcript CDS (c.) interval"""
        return self.r_to_c(self.g_to_r(g_interval))

    def c_to_g(self, c_interval):
        """convert a transcript CDS (c.) interval to a genomic (g.) interval"""
        return self.r_to_g(self.c_to_r(c_interval))



    ############################################################################
    ## DEPRECATED METHODS
    @deprecated(use_instead='g_to_c(...)')
    def hgvsg_to_hgvsc(self,*args,**kwargs):
        return self.g_to_c(*args,**kwargs)
    @deprecated(use_instead='g_to_r(...)')
    def hgvsg_to_hgvsr(self,*args,**kwargs):
        return self.g_to_r(*args,**kwargs)
    @deprecated(use_instead='r_to_g(...)')
    def hgvsr_to_hgvsg(self,*args,**kwargs):
        return self.r_to_g(*args,**kwargs)
    @deprecated(use_instead='c_to_g(...)')
    def hgvsc_to_hgvsg(self,*args,**kwargs):
        return self.c_to_g(*args,**kwargs)
    @deprecated(use_instead='c_to_r(...)')
    def hgvsc_to_hgvsr(self,*args,**kwargs):
        return self.c_to_r(*args,**kwargs)
    @deprecated(use_instead='r_to_c(...)')
    def hgvsr_to_hgvsc(self,*args,**kwargs):
        return self.r_to_c(*args,**kwargs)
    @deprecated(use_instead='c_to_p(...)')
    def hgvsc_to_hgvsp(self,*args,**kwargs):
        return self.c_to_p(*args,**kwargs)



def _ci_to_hgvs_coord(s, e):
    """ Convert continuous interbase (right-open) coordinates (..,-2,-1,0,1,..) to
    discontinuous HGVS coordinates (..,-2,-1,1,2,..)
    """
    def _ci_to_hgvs(c):
        return c + 1 if c >= 0 else c
    return (None if s is None else _ci_to_hgvs(s),
            None if e is None else _ci_to_hgvs(e) - 1)

def _hgvs_coord_to_ci(s, e):
    """convert start,end interval in inclusive, discontinuous HGVS coordinates
    (..,-2,-1,1,2,..) to continuous interbase (right-open) coordinates
    (..,-2,-1,0,1,..)"""
    def _hgvs_to_ci(c):
        assert c != 0, 'received CDS coordinate 0; expected ..,-2,-1,1,1,...'
        return c-1 if c>0 else c
    return (None if s is None else _hgvs_to_ci(s),
            None if e is None else _hgvs_to_ci(e) + 1)


## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
