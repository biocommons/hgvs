# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

"""
hgvs.normalizer
"""


import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.validator
import hgvs.posedit
import hgvs.location

from hgvs.exceptions import HGVSDataNotAvailableError, HGVSValidationError, HGVSNormalizationError




class Normalizer(object):
    """Perform variant normalization
    """
    
    def __init__(self, hdp, direction=3, fill=True):
        assert direction==3 or direction==5, "The shuffling direction should be 3 (3' most) or 5 (5' most)."
        self.hdp = hdp
        self.hm  = hgvs.variantmapper.VariantMapper(self.hdp)
        self.direction = direction      #shuffling direction
        self.fill = fill                #fill in nucleotides or strip nucleotides for delins and dup
    
    
    
    
    def _get_ref_context(self, var):
        """Get reference sequence around the variant using dynamic window size
        """
        var_x = self.hm.c_to_r(var) if var.type == 'c' else var
        #var start
        start = var_x.posedit.pos.start.base
        #var length
        length = var_x.posedit.pos.end.base - var_x.posedit.pos.start.base + 1
        
        #Get reference sequence. The default window size is three times of the
        #var length, which is left and right extend by the var length.
        left_pos = start - length
        right_pos = start + length + length - 1
        reference = self._fetch_seq(var_x.type, var_x.ac, left_pos, right_pos)
        base = left_pos if left_pos > 1 else 1
        
        #Extend the reference window if the edge of the window locates in homopolymers or repeats
        if self.direction == 3:
            #For 3' end
            while True:
                left_pos = right_pos + 1
                right_pos = left_pos + length - 1
                tmp_seq = self._fetch_seq(var_x.type, var_x.ac, left_pos, right_pos)
                if tmp_seq == '': break
                tmp_seq_len = len(tmp_seq)
                ref_len = len(reference)
                i=0
                while i < tmp_seq_len:
                    if reference[ref_len-tmp_seq_len+i : ref_len] == tmp_seq[0 : tmp_seq_len-i]:
                        break
                    i += 1
                if i == tmp_seq_len: break
                reference = reference + tmp_seq[0 : tmp_seq_len-i]
                right_pos = left_pos + tmp_seq_len - i - 1
        if self.direction == 5:
            #For 5' end
            while True:
                right_pos = left_pos - 1
                left_pos = right_pos - length + 1
                tmp_seq = self._fetch_seq(var_x.type, var_x.ac, left_pos, right_pos)
                if tmp_seq == '': break
                tmp_seq_len = len(tmp_seq)
                i=0
                while i < tmp_seq_len:
                    if reference[0 : tmp_seq_len-i] == tmp_seq[i : tmp_seq_len]:
                        break
                    i += 1
                if i == tmp_seq_len: break
                reference = tmp_seq[i : tmp_seq_len] + reference
                left_pos = right_pos - tmp_seq_len + i + 1
            base = right_pos + 1 if right_pos + 1 > 1 else 1
        
        return list(reference), base
    
    
    
    
    
    
    
    def _build_alignment(self, var, reference, base):
        """Reconstruct the variant sequence against the reference sequence
        """
        input = reference[::]
        if var.posedit.edit.type != 'dup':
            pos = var.posedit.pos.start.base - base
        else:
            pos = var.posedit.pos.end.base - base
        
        if var.posedit.edit.type == 'sub' or var.posedit.edit.type == 'delins':
            ref_len = var.posedit.pos.end.base - var.posedit.pos.start.base + 1
            alt_len = len(var.posedit.edit.alt)
            #validate whether the ref of the var is the same as the reference sequence
            if var.posedit.edit.ref_s is not None and var.posedit.edit.ref != ''.join(reference[pos : pos+ref_len]):
                raise HGVSValidationError(str(var) + ': ' + hgvs.validator.SEQ_ERROR_MSG)
            #reconstruct the variant alignment
            if alt_len == ref_len:
                input[pos : pos+ref_len] = list(var.posedit.edit.alt)
            elif alt_len < ref_len:
                input[pos : pos+ref_len] = ['-'] * (ref_len - alt_len) + list(var.posedit.edit.alt)
            elif alt_len > ref_len:
                input[pos : pos+ref_len] = list(var.posedit.edit.alt)
                reference[pos : pos+1] = ['-'] * (alt_len - ref_len) + reference[pos : pos+1]
            
        elif var.posedit.edit.type == 'del':
            ref_len = var.posedit.edit.ref_n
            if ref_len is None:
                ref_len = var.posedit.pos.end.base - var.posedit.pos.start.base + 1
            #validate whether the ref of the var is the same as the reference sequence
            if var.posedit.edit.ref_s is not None and var.posedit.edit.ref_s != ''.join(reference[pos : pos+ref_len]):
                raise HGVSValidationError(str(var) + ': ' + hgvs.validator.SEQ_ERROR_MSG)
            #reconstruct the variant alignment
            input[pos : pos+ref_len] = ['-'] * ref_len
            
        elif var.posedit.edit.type == 'ins':
            alt_len = len(var.posedit.edit.alt)
            input[pos : pos+1] = input[pos : pos+1] + list(var.posedit.edit.alt)
            reference[pos : pos+1] = reference[pos :pos+1] + ['-'] * alt_len
        
        elif var.posedit.edit.type == 'dup':
            dup_len = var.posedit.pos.end.base - var.posedit.pos.start.base + 1
            if var.posedit.edit.seq is not None and var.posedit.edit.seq != '' and var.posedit.edit.seq != ''.join(reference[pos-dup_len+1 : pos+1]):
                raise HGVSValidationError(str(var) + ': ' + hgvs.validator.SEQ_ERROR_MSG)
            input[pos : pos+1] = input[pos : pos+1] + reference[pos-dup_len+1 : pos+1]
            reference[pos : pos+1] = reference[pos :pos+1] + ['-'] * dup_len
        
        return input
    
    
    
    
    
    
    def _fetch_seq(self, type, ac, start, end):
        """Fetch reference sequence from hgvs data provider
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
        
        start = start - 1 if start >= 1 else 0
        end = end if end <= len(seq) else len(seq)
        return seq[start:end]
    
    
    
    
    def _right_shuffle(self, input, reference):
        """Shuffling variants to 3' most position
        """
        i = 0
        while i < len(input)-1:
            if input[i] == '-':
                end = i + 1
                while end < len(input) and input[end] == input[i]:
                    end += 1
                if end == len(input):
                    break
                j = 0
                while end + j < len(input) and input[end + j] == reference[i + j]:
                    j += 1
                if j == 0:
                    i = end
                    continue
                tmp = input[end : end+j]
                input[i+j : end+j] = input[i : end]
                input[i : i+j] = tmp
                i = end + j
            else:
                i += 1
        i = 0
        while i < len(reference)-1:
            if reference[i] == '-':
                end = i + 1
                while end < len(reference) and reference[end] == reference[i]:
                    end += 1
                if end == len(reference):
                    break
                j = 0
                while end + j < len(reference) and reference[end + j] == input[i + j]:
                    j += 1
                if j == 0:
                    i = end
                    continue
                tmp = reference[end : end+j]
                reference[i+j : end+j] = reference[i : end]
                reference[i : i+j] = tmp
                i = end + j
            else:
                i += 1
        
        return input, reference
    
    
    
    
    
    def _left_shuffle(self, input, reference):
        """Shuffling variants to 5' most position
        """
        input.reverse()
        reference.reverse()
        input, reference = self._right_shuffle(input, reference)
        reference.reverse()
        input.reverse()
        return input, reference
    
    
    
    
    
    def normalize(self, var):
        """Perform variants canonicalization
        """
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        
        if var.posedit.uncertain:
            return var
        
        reference, base = self._get_ref_context(var)
        input = self._build_alignment(var, reference, base)
        if self.direction == 3:
            input, reference = self._right_shuffle(input, reference)
        elif self.direction == 5:
            input, reference = self._left_shuffle(input, reference)
        
        #Remove positions that both input and reference are null
        i = 0
        while i < len(reference):
            if input[i] == '-' and reference[i] == '-':
                del input[i]
                del reference[i]
                if i == 0:
                    base += 1
            else:
                i += 1
        
        
        #Left trim of common base pairs
        start = 0
        while start < len(reference) and input[start] == reference[start]:
            start += 1
        
        #Right trim of common base pairs
        end = len(reference) - 1
        while end >= start and input[end] == reference[end]:
            end -= 1
        
        #ref and alt are identity
        if start > end:
            return var
        
        ref = reference[start : end + 1]
        alt = input[start : end + 1]
        while len(ref) != 0 and ref[0] == '-':
            del ref[0]
        while len(ref) != 0 and ref[len(ref) - 1] == '-':
            del ref[len(ref) - 1]
        while len(alt) != 0 and alt[0] == '-':
            del alt[0]
        while len(alt) != 0 and alt[len(alt) - 1] == '-':
            del alt[len(alt) - 1]
        
        ref = ''.join(ref)
        alt = ''.join(alt)
        ref_len = len(ref)
        alt_len = len(alt)
        
        
        #Generate normalized variant        
        if alt_len == ref_len:
            ref_start = base + start
            ref_end   = base + end
            #substitution
            if start == end:
                edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
            #delins
            else:
                if self.fill:
                    edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
                else:
                    edit = hgvs.edit.NARefAlt(ref='', alt=alt)
        elif alt_len < ref_len:
            ref_start = base + start
            ref_end   = base + end
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
                #dup
                if alt == ''.join(reference[start-alt_len : start]):
                    ref_start = base + start - alt_len
                    ref_end   = base + start - 1
                    if self.fill:
                        edit = hgvs.edit.Dup(seq=alt)
                    else:
                        edit = hgvs.edit.Dup(seq='')
                #ins
                else:
                    ref_start = base + start -1
                    ref_end   = ref_start + 1
                    edit = hgvs.edit.NARefAlt(ref=None, alt=alt)
            #delins
            else:
                ref_start = base + start
                ref_end = ref_start + ref_len -1
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
    var = hgvsparser.parse_hgvs_variant('NM_001166478.1:c.31del')
    hdp = hgvs.dataproviders.uta.connect()
    norm = Normalizer(hdp, direction=3)
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
