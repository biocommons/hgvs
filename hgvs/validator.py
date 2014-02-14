"""
hgvs.hgvsvalidator
"""
from hgvs.exceptions import HGVSValidationError

import hgvs.hgvsmapper
import hgvs.parser


BASE_RANGE_ERROR_MSG = 'base start position must be <= end position'
OFFSET_RANGE_ERROR_MSG = 'offset start must be <= end position'
INS_ERROR_MSG = 'insertion length must be 1'
DEL_ERROR_MSG = 'start and end position range must equal sequence deletion length'
SUB_ERROR_MSG = 'substitution cannot have the same ref and alt base'
AC_ERROR_MSG = 'Accession is not present in BDI database'
SEQ_ERROR_MSG = 'Ref variant does not agree with reference sequence/transcript'


def validate(var):
    """Validates a parsed HGVS variant using intrinsic and extrinsic methods"""
    return Validator().validate(var)


class Validator(object):
    """invoke intrinsic and extrinsic validation"""
    def __init__(self, bdi, mfdb):
        self.ivr = IntrinsicValidator()
        self.evr = ExtrinsicValidator(bdi,mfdb)

    def validate(self,var):
        return self.ivr.validate(var) and self.evr.validate(var)


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

    def __init__(self, bdi, mfdb):
        self.bdi = bdi
        self.mfdb = mfdb
        self.hm = hgvs.hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self._ac_is_valid(var)
        self._ref_is_valid(var)
        return True

    def _ac_is_valid(self, var):
        ac = self.bdi.get_tx_info(var.ac)
        if ac is None:
            raise HGVSValidationError(AC_ERROR_MSG)
        elif ac['ac'] != var.ac:
            raise HGVSValidationError(AC_ERROR_MSG)
        else:
            return True

    def _ref_is_valid(self, var):
        if var.posedit.edit.ref is not None:
            if var.type == 'c':
                var = self.hm.hgvsc_to_hgvsr(var)
            # fetch uses interbase coordinates
            seq = self.mfdb.fetch(var.ac, var.posedit.pos.start.base - 1, var.posedit.pos.end.base)
            if seq != var.posedit.edit.ref:
                raise HGVSValidationError(SEQ_ERROR_MSG)
        return True


if __name__ == '__main__':
    hgvsparser = hgvs.parser.Parser()
    var1 = hgvsparser.parse_hgvs_variant('NM_001005405.2:r.2T>A')
    validate_ext = ExtrinsicValidator()
    validate_ext.validate(var1)

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
