# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
"""
hgvs.hgvsvalidator
"""

from .exceptions import HGVSValidationError, HGVSUnsupportedOperationError

import hgvs.parser
import hgvs.variantmapper

BASE_RANGE_ERROR_MSG = 'base start position must be <= end position'
OFFSET_RANGE_ERROR_MSG = 'offset start must be <= end position'
INS_ERROR_MSG = 'insertion length must be 1'
DEL_ERROR_MSG = 'Length implied by coordinates ({span_len})  must equal sequence deletion length ({del_len})'
AC_ERROR_MSG = 'Accession is not present in BDI database'
SEQ_ERROR_MSG = 'Variant reference ({var_ref_seq}) does not agree with reference sequence ({ref_seq})'


# TODO: #249: redisign validation interface for greater flexibility

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
        return True

    def _start_lte_end(self, var):
        if var.type == 'g':
            if var.posedit.pos.start.base > var.posedit.pos.end.base:
                raise HGVSValidationError(BASE_RANGE_ERROR_MSG)
        if var.type in 'cmnp' and var.posedit.pos:
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
            if var.type in 'cmnp':
                if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                    (var.posedit.pos.start.base + var.posedit.pos.start.offset)) != 1:
                    raise HGVSValidationError(INS_ERROR_MSG)
            return True

    def _del_length(self, var):
        if var.posedit.edit.type in ['del', 'delins']:
            ref_len = var.posedit.edit.ref_n
            if ref_len is None:
                return True

            if var.type in 'cnr':
                assert ((var.posedit.pos.start.offset == var.posedit.pos.end.offset == 0) or
                        (var.posedit.pos.start.base == var.posedit.pos.end.base))
                span_len = ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                            (var.posedit.pos.start.base + var.posedit.pos.start.offset) + 1)
            else:
                span_len = var.posedit.pos.end.base - var.posedit.pos.start.base + 1

            if span_len != ref_len:
                raise HGVSValidationError(DEL_ERROR_MSG.format(span_len=span_len, del_len=ref_len))
        return True


class ExtrinsicValidator():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """

    def __init__(self, hdp):
        self.hdp = hdp
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp)

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self._ref_is_valid(var)
        return True

    def _ref_is_valid(self, var):
        # use reference sequence of original variant, even if later converted (eg c_to_n)
        if (var.type in 'cnr'
            and var.posedit.pos is not None
            and (var.posedit.pos.start.offset != 0 or var.posedit.pos.end.offset != 0)):
            raise HGVSUnsupportedOperationError(
                "Cannot validate sequence of an intronic variant ({})".format(str(var)))

        var_ref_seq = getattr(var.posedit.edit, 'ref', None)

        if var_ref_seq == '':
            var_ref_seq = None

        if var_ref_seq:
            # ref_seq is digit, as in 'del6'
            try:
                int(var_ref_seq)
                var_ref_seq = None
            except ValueError:
                pass

        if var_ref_seq is not None:
            var_x = self.vm.c_to_n(var) if var.type == 'c' else var
            ref_seq = self.hdp.fetch_seq(var_x.ac, var_x.posedit.pos.start.base - 1, var_x.posedit.pos.end.base)
            if ref_seq != var_ref_seq:
                raise HGVSValidationError(str(var) + ': ' + SEQ_ERROR_MSG.format(ref_seq=ref_seq,
                                                                                 var_ref_seq=var_ref_seq))

        return True


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
