# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

"""
hgvs.normalizer
"""


import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.validator
import hgvs.posedit
import hgvs.location

from hgvs.exceptions import HGVSDataNotAvailableError, HGVSValidationError, HGVSNormalizationError


from collections import namedtuple
from itertools   import chain



class Normalizer(object):
    """Perform variant normalization
    """
    
    def __init__(self, hdp, direction=3, fill=True):
        assert direction==3 or direction==5, "The shuffling direction should be 3 (3' most) or 5 (5' most)."
        self.hdp = hdp
        self.direction = direction      #shuffling direction
        self.fill = fill                #fill in nucleotides or strip nucleotides for delins and dup
    

    
    #Code of the normalization utilities were imported from vgraph
    #https://github.com/bioinformed/vgraph
    def _trim_common_suffixes(self, strs, min_len=0):
        '''trim common suffixes'''
        
        if len(strs) < 2:
            return 0, strs
        
        rev_strs = [ s[::-1] for s in strs ]
        
        trimmed, rev_strs = self._trim_common_prefixes(rev_strs, min_len)
        
        if trimmed:
            strs = [ s[::-1] for s in rev_strs ]
        
        return trimmed, strs
    
    
    def _trim_common_prefixes(self, strs, min_len=0):
        '''trim common prefixes'''
        
        trimmed = 0
        
        if len(strs) > 1:
            s1 = min(strs)
            s2 = max(strs)
            
            for i in range(len(s1) - min_len):
                if s1[i] != s2[i]:
                    break
                trimmed = i + 1
        
        if trimmed > 0:
            strs = [ s[trimmed:] for s in strs ]
        
        return trimmed, strs
    
    
    
    def _normalize_alleles_left(self, ref, start, stop, alleles, bound, ref_step, shuffle=True):
        '''Normalize loci by removing extraneous reference padding'''
        
        normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')
        
        if len(alleles) < 2 or start <= 0 or stop <= 0:
            return normalized_alleles(start, stop, alleles)
        
        # STEP 1: Trim common suffix
        trimmed, alleles = self._trim_common_suffixes(alleles)
        stop -= trimmed
        
        # STEP 2: Trim common prefix
        trimmed, alleles = self._trim_common_prefixes(alleles)
        start += trimmed
        
        #assert bound <= start,'start={:d}, left bound={:d}'.format(start, bound)
        
        # STEP 3: While a null allele exists, left shuffle by prepending alleles
        #         with reference and trimming common suffixes
        while shuffle and '' in alleles and start > bound:
            step = min(ref_step, start - bound)
            
            r = ref[start-step : start].upper()
            new_alleles = [ r+a for a in alleles ]
            
            trimmed, new_alleles = self._trim_common_suffixes(new_alleles)
            
            if not trimmed:
                break
            
            start -= trimmed
            stop  -= trimmed
            
            if trimmed == step:
                alleles = new_alleles
            else:
                left    = step - trimmed
                alleles = [ a[left:] for a in new_alleles ]
                break
        
        return normalized_alleles(start, stop, tuple(alleles))
    
    
    def _normalize_alleles_right(self, ref, start, stop, alleles, bound, ref_step, shuffle=True):
        '''Normalize loci by removing extraneous reference padding'''
        
        normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')
        
        chrom_stop = len(ref)
        
        if len(alleles) < 2 or stop >= chrom_stop:
            return normalized_alleles(start, stop, alleles)
        
        # STEP 1: Trim common prefix
        trimmed, alleles = self._trim_common_prefixes(alleles)
        start += trimmed
        
        # STEP 2: Trim common suffix
        trimmed, alleles = self._trim_common_suffixes(alleles)
        stop -= trimmed
        
        #assert bound >= stop,'stop={:d}, right bound={:d}'.format(stop, bound)
        
        # STEP 3: While a null allele exists, right shuffle by appending alleles
        #         with reference and trimming common prefixes
        while shuffle and '' in alleles and stop < bound:
            step = min(ref_step, bound - stop)
            
            r = ref[stop:stop+step].upper()
            new_alleles = [ a+r for a in alleles ]
            
            trimmed, new_alleles = self._trim_common_prefixes(new_alleles)
            
            if not trimmed:
                break
            
            start += trimmed
            stop  += trimmed
            
            if trimmed == step:
                alleles = new_alleles
            else:
                left    = step - trimmed
                alleles = [ a[:-left] for a in new_alleles ]
                break
        
        return normalized_alleles(start, stop, tuple(alleles))
    #End of code from vgraph
    
    
    
    
    
    def _fetch_seq(self, type, ac, start, end):
        """Fetch reference sequence from hgvs data provider.
           
           The start posotion is 0 and the interval is half open
        """
        try:
            if type == 'r' or type == 'c' or type == 'n':
                seq = self.hdp.get_tx_seq(ac)
                
            elif type == 'g':
                raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))
                
            elif type == 'p':
                raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))
                
            elif type == 'm':
                raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))
                
        except TypeError:
            raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))
        
        start = start if start >= 0 else 0
        end = end if end <= len(seq) else len(seq)
        return seq[start:end]
    
    
    
    
    def _get_ref_alt(self, var):
        """Get reference allele and alternative allele of the variant
        """
        
        #Get reference allele
        if var.posedit.edit.type == 'ins':
            ref = ''
        elif var.posedit.edit.type == 'dup':
            if var.posedit.edit.seq:
                ref = self._fetch_seq(var.type, var.ac, var.posedit.pos.start.base-1, var.posedit.pos.end.base)
                #validate whether the ref of the var is the same as the reference sequence
                if var.posedit.edit.seq != ref:
                    raise HGVSValidationError(str(var) + ': ' + hgvs.validator.SEQ_ERROR_MSG)
            ref = ''
        else:
            ref = self._fetch_seq(var.type, var.ac, var.posedit.pos.start.base-1, var.posedit.pos.end.base)
            #validate whether the ref of the var is the same as the reference sequence
            if var.posedit.edit.ref_s is not None and var.posedit.edit.ref != '' and var.posedit.edit.ref != ref:
                raise HGVSValidationError(str(var) + ': ' + hgvs.validator.SEQ_ERROR_MSG)
        ref = ref.encode('ascii')
        
        #Get alternative allele
        if var.posedit.edit.type == 'sub' or var.posedit.edit.type == 'delins' or var.posedit.edit.type == 'ins':
            alt = var.posedit.edit.alt
        elif var.posedit.edit.type == 'del':
            alt = ''
        elif var.posedit.edit.type == 'dup':
            alt = var.posedit.edit.seq or self._fetch_seq(var.type, var.ac, var.posedit.pos.start.base-1, var.posedit.pos.end.base)
        elif var.posedit.edit.type == 'inv':
            alt = ref[::-1]
        alt = alt.encode('ascii')
        
        return ref, alt
    
    
    
    
    
    def _normalize_alleles(self, var):
        """Normalize the variant until it could not be shuffled
        """
        
        ref, alt = self._get_ref_alt(var)
        
        win_size = max(len(ref), len(alt)) * 2
        
        if self.direction == 3:
            if var.posedit.edit.type == 'ins':
                base  = var.posedit.pos.start.base
                start = 1
                stop  = 1
            elif var.posedit.edit.type == 'dup':
                base  = var.posedit.pos.end.base
                start = 1
                stop  = 1
            else:
                base  = var.posedit.pos.start.base
                start = 0
                stop  = var.posedit.pos.end.base - base + 1
            
            while True:
                ref_seq = self._fetch_seq(var.type, var.ac, base - 1, base -1 + win_size)
                if ref_seq == '':
                    break
                orig_start, orig_stop = start, stop
                start, stop, (ref, alt) = self._normalize_alleles_right(ref_seq, start, stop, (ref, alt), len(ref_seq), 1)
                if stop < len(ref_seq):
                    break
                #if stop at the end of the window, try to extend the shuffling to the right
                base += start - orig_start
                stop -= start - orig_start
                start = orig_start
            
        elif self.direction == 5:
            if var.posedit.edit.type == 'ins':
                base  = var.posedit.pos.end.base - win_size + 1
                start = win_size - 1
                stop  = win_size - 1
            elif var.posedit.edit.type == 'dup':
                base  = var.posedit.pos.end.base - win_size + 1
                start = win_size
                stop  = win_size
            else:
                base  = var.posedit.pos.end.base - win_size + 1
                start = var.posedit.pos.start.base - base
                stop  = win_size
            
            while True:
                ref_seq = self._fetch_seq(var.type, var.ac, base - 1, base -1 + win_size)
                if ref_seq == '':
                    break
                orig_start, orig_stop = start, stop
                start, stop, (ref, alt) = self._normalize_alleles_left(ref_seq, start, stop, (ref, alt), 0, 1)
                if start > 0:
                    break
                #if stop at the end of the window, try to extend the shuffling to the left
                base -= orig_stop - stop
                start += orig_stop - stop
                stop = orig_stop
        
        return base + start, base + stop, (ref, alt)
    
    
    
    
    
    def normalize(self, var):
        """Perform variants canonicalization
        """
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        
        if var.posedit.uncertain:
            return var
        
        start, end, (ref, alt) = self._normalize_alleles(var)
        
        ref_len = len(ref)
        alt_len = len(alt)
        
        
        #Generate normalized variant        
        if alt_len == ref_len:
            ref_start = start
            ref_end   = end - 1
            #substitution
            if start == end - 1:
                edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
            #delins
            else:
                if self.fill:
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
                else:
                    edit = hgvs.edit.NARefAlt(ref='', alt=alt)
        elif alt_len < ref_len:
            ref_start = start
            ref_end   = end - 1
            #del
            if alt_len == 0:
                if self.fill:
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=None)
                else:
                    edit = hgvs.edit.NARefAlt(ref='', alt=None)
            #delins
            else:
                if self.fill:
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
                else:
                    edit = hgvs.edit.NARefAlt(ref='', alt=alt)
        elif alt_len > ref_len:
            #ins or dup
            if ref_len == 0:
                left_seq  = self._fetch_seq(var.type, var.ac, start-alt_len-1, end-1) if self.direction == 3 else ''
                right_seq = self._fetch_seq(var.type, var.ac, start-1, start+alt_len-1) if self.direction == 5 else ''
                #dup
                if alt == left_seq:
                    ref_start = start - alt_len
                    ref_end   = end - 1
                    if self.fill:
                        edit = hgvs.edit.Dup(seq=alt)
                    else:
                        edit = hgvs.edit.Dup(seq='')
                elif alt == right_seq:
                    ref_start = start
                    ref_end   = start + alt_len - 1
                    if self.fill:
                        edit = hgvs.edit.Dup(seq=alt)
                    else:
                        edit = hgvs.edit.Dup(seq='')
                #ins
                else:
                    ref_start = start - 1
                    ref_end   = end
                    edit = hgvs.edit.NARefAlt(ref=None, alt=alt)
            #delins
            else:
                ref_start = start
                ref_end   = end -1
                if self.fill:
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
                else:
                    edit = hgvs.edit.NARefAlt(ref='', alt=alt)
        
        
        var_ac    = var.ac
        var_type  = var.type
        var_start = hgvs.location.SimplePosition(base=ref_start)
        var_end   = hgvs.location.SimplePosition(base=ref_end)
        var_pos   = hgvs.location.Interval(start=var_start, end=var_end)
        var_posedit = hgvs.posedit.PosEdit(pos=var_pos, edit=edit)
        return hgvs.variant.SequenceVariant(ac=var_ac, type=var_type, posedit=var_posedit)
    






if __name__ == '__main__':
    hgvsparser = hgvs.parser.Parser()
    var = hgvsparser.parse_hgvs_variant('NM_001166478.1:c.30_31insT')
    hdp = hgvs.dataproviders.uta.connect()
    norm = Normalizer(hdp, direction=5)
    res  = norm.normalize(var)
    print(str(var) + '    =>    ' + str(res))





## <LICENSE>
## Copyright 2015 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
