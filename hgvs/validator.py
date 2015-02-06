# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

"""
hgvs.hgvsvalidator
"""

from hgvs.exceptions import HGVSValidationError, HGVSDataNotAvailableError
import hgvs.parser
import hgvs.variantmapper


BASE_RANGE_ERROR_MSG = 'base start position must be <= end position'
OFFSET_RANGE_ERROR_MSG = 'offset start must be <= end position'
INS_ERROR_MSG = 'insertion length must be 1'
DEL_ERROR_MSG = 'start and end position range must equal sequence deletion length'
SUB_ERROR_MSG = 'substitution cannot have the same ref and alt base'
AC_ERROR_MSG = 'Accession is not present in BDI database'
SEQ_ERROR_MSG = 'Variant reference does not agree with reference sequence'



class Validator(object):
    """invoke intrinsic and extrinsic validation"""
    def __init__(self, hdp):
        self._ivr = IntrinsicValidator()
        self._evr = ExtrinsicValidator(hdp)

    def validate(self, var):
        return self._ivr.validate(var) and self._evr.validate(var)


class IntrinsicValidator(object):
    """
    Attempts to determine if the HGVS name is internally consistent
    """

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self._start_lte_end(var)
        self._ins_length_is_one(var)
        self._del_length(var)
        self._sub_alt_is_not_ref(var)
        return True

    def _start_lte_end(self, var):
        if var.type == 'g':
            if var.posedit.pos.start.base > var.posedit.pos.end.base:
                raise HGVSValidationError(BASE_RANGE_ERROR_MSG)
        if var.type in ['c', 'r', 'p']:
            if var.posedit.pos.start.base > var.posedit.pos.end.base:
                raise HGVSValidationError(BASE_RANGE_ERROR_MSG)
            elif var.posedit.pos.start.base == var.posedit.pos.end.base:
                if var.posedit.pos.start.offset > var.posedit.pos.end.offset:
                    raise HGVSValidationError(OFFSET_RANGE_ERROR_MSG)
            if var.posedit.pos.start.datum > var.posedit.pos.end.datum:
                raise HGVSValidationError(BASE_RANGE_ERROR_MSG)
        return True

    def _ins_length_is_one(self, var):
        if var.posedit.edit.type == 'ins':
            if var.type == 'g':
                if (var.posedit.pos.end.base - var.posedit.pos.start.base) != 1:
                    raise HGVSValidationError(INS_ERROR_MSG)
            if var.type in ['c', 'r', 'p']:
                if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                        (var.posedit.pos.start.base + var.posedit.pos.start.offset)) != 1:
                    raise HGVSValidationError(INS_ERROR_MSG)
            return True

    def _del_length(self, var):
        if var.posedit.edit.type == 'del':
            del_len = len(var.posedit.edit.ref)
            if var.type == 'g':
                if (var.posedit.pos.end.base - var.posedit.pos.start.base + 1) != del_len:
                    raise HGVSValidationError(DEL_ERROR_MSG)
            if var.type in ['c', 'r', 'p']:
                if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                        (var.posedit.pos.start.base + var.posedit.pos.start.offset) + 1) != del_len:
                    raise HGVSValidationError(DEL_ERROR_MSG)
            return True

    def _sub_alt_is_not_ref(self, var):
        if var.posedit.edit.type == 'identity':
            if var.posedit.edit.ref == var.posedit.edit.alt:
                raise HGVSValidationError(SUB_ERROR_MSG)


class ExtrinsicValidator():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """

    def __init__(self, hdp):
        self.hdp = hdp
        self.hm = hgvs.variantmapper.VariantMapper(self.hdp)

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self._ref_is_valid(var)
        return True

    def _ref_is_valid(self, var):
        var_ref_seq = None
        if var.posedit.edit.type == 'dup':
            # Handle Dup and NADupN objects.
            var_ref_seq = getattr(var.posedit.edit, 'seq', None)
        else:
            # use reference sequence of original variant, even if later converted (eg c_to_r)
            var_ref_seq = getattr(var.posedit.edit, 'ref', None)

        if var_ref_seq is not None:
            var_x = self.hm.c_to_r(var) if var.type == 'c' else var
            ref_seq = self._fetch_seq(var_x.ac, var_x.posedit.pos.start.base - 1, var_x.posedit.pos.end.base)
            if ref_seq != var_ref_seq:
                raise HGVSValidationError(str(var) + ': ' + SEQ_ERROR_MSG)
        return True

    def _fetch_seq(self, ac, start_i, end_i):
        """fetch 0-based, right-open subsequence [start_i,end_i) of specified accession
        first using the hdpi instance.

        :raises: hgvs.exceptions.HGVSDataNotAvailableError if couldn't get sequence
        """
        try:
            return self.hdp.get_tx_seq(ac)[start_i:end_i]
        except TypeError:
            raise HGVSDataNotAvailableError("No sequence available for {ac}".format(ac=ac))



if __name__ == '__main__':
    hgvsparser = hgvs.parser.Parser()
    var1 = hgvsparser.parse_hgvs_variant('NM_001005405.2:r.2T>A')
    validate_ext = ExtrinsicValidator()
    validate_ext.validate(var1)


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
