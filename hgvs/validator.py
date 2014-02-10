"""
hgvs.hgvsvalidator
"""
from hgvs.exceptions import HGVSValidationError
from bdi.sources import uta0
from bdi.multifastadb import MultiFastaDB

import hgvsmapper
import hgvs.parser


def validate(var):
    """Validates a parsed HGVS variant using intrinsic and extrinsic methods"""
    assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
    iv = IntrinsicValidator()
    ev = ExtrinsicValidator()
    if iv.validate(var) and ev.validate(var):
        return True
    else:
        return False


class IntrinsicValidator(object):
    """
    Attempts to determine if the HGVS name is internally consistent
    """
    BASE_RANGE_ERROR_MSG = 'base start position must be <= end position'
    OFFSET_RANGE_ERROR_MSG = 'offset start must be <= end position'
    INS_ERROR_MSG = 'insertion length must be 1'
    DEL_ERROR_MSG = 'start and end position range must equal sequence deletion length'
    SUB_ERROR_MSG = 'substitution cannot have the same ref and alt base'

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
                raise HGVSValidationError(self.BASE_RANGE_ERROR_MSG)
        if var.type in ['c', 'r', 'p']:
            if var.posedit.pos.start.base > var.posedit.pos.end.base:
                raise HGVSValidationError(self.BASE_RANGE_ERROR_MSG)
            elif var.posedit.pos.start.base == var.posedit.pos.end.base:
                if var.posedit.pos.start.offset > var.posedit.pos.end.offset:
                    raise HGVSValidationError(self.OFFSET_RANGE_ERROR_MSG)
            if var.posedit.pos.start.datum > var.posedit.pos.end.datum:
                raise HGVSValidationError(self.BASE_RANGE_ERROR_MSG)
        return True
    
    def _ins_length_is_one(self, var):
        if var.posedit.edit.type == 'ins':
            if var.type == 'g':
                if (var.posedit.pos.end.base - var.posedit.pos.start.base) != 1:
                    raise HGVSValidationError(self.INS_ERROR_MSG)
            if var.type in ['c', 'r', 'p']:
                if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                        (var.posedit.pos.start.base + var.posedit.pos.start.offset)) != 1:
                    raise HGVSValidationError(self.INS_ERROR_MSG)
            return True
        
    def _del_length(self, var):
        if var.posedit.edit.type == 'del':
            del_len = len(var.posedit.edit.ref)
            if var.type == 'g':
                if (var.posedit.pos.end.base - var.posedit.pos.start.base + 1) != del_len:
                    raise HGVSValidationError(self.DEL_ERROR_MSG)
            if var.type in ['c', 'r', 'p']:
                if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                        (var.posedit.pos.start.base + var.posedit.pos.start.offset) + 1) != del_len:
                    raise HGVSValidationError(self.DEL_ERROR_MSG)
            return True

    def _sub_alt_is_not_ref(self, var):
        if var.posedit.edit.type == 'identity':
            if var.posedit.edit.ref == var.posedit.edit.alt:
                raise HGVSValidationError(self.SUB_ERROR_MSG)


class ExtrinsicValidator():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """
    AC_ERROR_MSG = 'Accession is not present in BDI database'
    SEQ_ERROR_MSG = 'Ref variant does not agree with reference sequence/transcript'

    def __init__(self, bdi=None, mfdb=None):
        # optional args for bdi and mfdb (users may already have them)
        if bdi is None:
            self.bdi = uta0.connect()
        else:
            self.bdi = bdi
        if mfdb is None:
            # specify path to data files
            db_dir = ['tests/data/sample_data']  # local sample data - replace if necessary with path to real data
            self.mfdb = MultiFastaDB(db_dir, use_meta_index=True)
        else:
            self.mfdb = mfdb
        self.hm = hgvsmapper.HGVSMapper(self.bdi, cache_transcripts=True)

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self._ac_is_valid(var)
        self._ref_is_valid(var)
        return True

    def _ac_is_valid(self, var):
        ac = self.bdi.get_tx_info(var.ac)
        if ac is None:
            raise HGVSValidationError(self.AC_ERROR_MSG)
        elif ac['ac'] != var.ac:
            raise HGVSValidationError(self.AC_ERROR_MSG)
        else:
            return True

    def _ref_is_valid(self, var):
        if var.posedit.edit.ref is not None:
            if var.type == 'c':
                var = self.hm.hgvsc_to_hgvsr(var)
            # fetch uses interbase coordinates
            seq = self.mfdb.fetch(var.ac, var.posedit.pos.start.base - 1, var.posedit.pos.end.base)
            if seq != var.posedit.edit.ref:
                raise HGVSValidationError(self.SEQ_ERROR_MSG)
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