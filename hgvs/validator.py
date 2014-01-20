"""
hgvs.hgvsvalidator
"""

import hgvs.parser


class Validate(object):
    """Validates a parsed HGVS variant using intrinsic and extrinsic methods"""
    pass


class IntrinsicValidation(object):
    """
    Attempts to determine if the HGVS name is internally consistent
    """
    BASE_RANGE_ERROR_MSG = 'base start position must be <= end position'
    OFFSET_RANGE_ERROR_MSG = 'offset start must be <= end position'
    INS_ERROR_MSG = 'insertion length must be 1'
    DEL_ERROR_MSG = 'start and end position range must equal sequence deletion length'

    def __init__(self):
        self.var = None

    def validate(self, var):
        assert isinstance(var, hgvs.variant.SequenceVariant), 'variant must be a parsed HGVS sequence variant object'
        self.var = var
        self._start_lte_end()
        self._ins_length_is_one()
        self._del_length()
        return True

    def _start_lte_end(self):
        if self.var.type == 'g':
            if self.var.posedit.pos.start.base > self.var.posedit.pos.end.base:
                raise Exception(self.BASE_RANGE_ERROR_MSG)
        if self.var.type in ['c', 'r', 'p']:
            if self.var.posedit.pos.start.base > self.var.posedit.pos.end.base:
                raise Exception(self.BASE_RANGE_ERROR_MSG)
            elif self.var.posedit.pos.start.base == self.var.posedit.pos.end.base:
                if self.var.posedit.pos.start.offset > self.var.posedit.pos.end.offset:
                    raise Exception(self.OFFSET_RANGE_ERROR_MSG)
            if self.var.posedit.pos.start.datum > self.var.posedit.pos.end.datum:
                raise Exception(self.BASE_RANGE_ERROR_MSG)
        return True
    
    def _ins_length_is_one(self):
        if self.var.posedit.edit.type == 'ins':
            if self.var.type == 'g':
                if (self.var.posedit.pos.end.base - self.var.posedit.pos.start.base) != 1:
                    raise Exception(self.INS_ERROR_MSG)
            if self.var.type in ['c', 'r', 'p']:
                if ((self.var.posedit.pos.end.base + self.var.posedit.pos.end.offset) -
                        (self.var.posedit.pos.start.base + self.var.posedit.pos.start.offset)) != 1:
                    raise Exception(self.INS_ERROR_MSG)
            return True
        
    def _del_length(self):
        if self.var.posedit.edit.type == 'del':
            del_len = len(self.var.posedit.edit.ref)
            if self.var.type == 'g':
                if (self.var.posedit.pos.end.base - self.var.posedit.pos.start.base + 1) != del_len:
                    raise Exception(self.DEL_ERROR_MSG)
            if self.var.type in ['c', 'r', 'p']:
                if ((self.var.posedit.pos.end.base + self.var.posedit.pos.end.offset) -
                        (self.var.posedit.pos.start.base + self.var.posedit.pos.start.offset) + 1) != del_len:
                    raise Exception(self.DEL_ERROR_MSG)
            return True



class ExtrinsicValidation():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """
    pass
    #def __init__(self):
    #    self.hp = hgvs.parser.Parser()
    #    self.ivalid = IntrinsicValidation()
    #
    #def valid_ac(self, hgvs):
    #    if self.ivalid.valid_parse(hgvs):
    #        var = self.hp.parse_hgvs_variant(hgvs)
    #        r = requests.get('http://uta.locusdev.net/api/v0/transcripts/fetch_seq?ac={ac}&start=0&end=10'
    #                        .format(ac=var.ac))
    #        if r.status_code == 200:
    #            return True
    #        else:
    #            return False
    #    else:
    #        return False
    #
    #def valid_ref(self, hgvs):
    #    if self.ivalid.valid_parse(hgvs):
    #        var = self.hp.parse_hgvs_variant(hgvs)
    #        if len(var.posedit.edit.ref) >= 1:
    #            r = requests.get('http://uta.locusdev.net/api/v0/transcripts/fetch_seq?ac={ac}&start={start}&end={end}'
    #                            .format(ac=var.ac, start=var.posedit.pos.start.base, end=var.posedit.pos.end.base))
    #            if r.status_code == 200:
    #                return True
    #            else:
    #                return False



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