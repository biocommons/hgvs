"""
hgvs.hgvsvalidator
"""

import requests
import hgvs.parser
import re


class IntrinsicValidation():
    """
    Attempts to determine if the HGVS name is internally consistent
    """
    def __init__(self):
        self.hp = hgvs.parser.Parser()

    def valid_parse(self, hgvs):
        try:
            self.hp.parse_hgvs_variant(hgvs)
        except Exception, e:
            #return 'Failed to parse {hgvs}: {error}'.format(hgvs=hgvs, error=e)
            return False
        return True

    def end_gt_start(self, hgvs):
        if self.valid_parse(hgvs):
            var = self.hp.parse_hgvs_variant(hgvs)
            if var.type == 'g':
                if var.posedit.pos.end.base >= var.posedit.pos.start.base:
                    return True
                else:
                    return False
            if var.type in ['c', 'r', 'p']:
                if (var.posedit.pos.end.base + var.posedit.pos.end.offset) >= \
                        (var.posedit.pos.start.base + var.posedit.pos.start.offset):
                    return True
                else:
                    return False

    def valid_ins(self, hgvs):
        if self.valid_parse(hgvs):
            var = self.hp.parse_hgvs_variant(hgvs)
            #? can I get the edit type somewhere?
            ins_re = re.compile('^ins')
            if len(ins_re.findall(str(var.posedit.edit))) == 1:
                # ins length should be 1
                if var.type == 'g':
                    if (var.posedit.pos.end.base - var.posedit.pos.start.base) == 1:
                        return True
                    else:
                        return False
                if var.type in ['c', 'r', 'p']:
                    if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                            (var.posedit.pos.start.base + var.posedit.pos.start.offset)) == 1:
                        return True
                    else:
                        return False
            else:
                return False

    def valid_del(self, hgvs):
        if self.valid_parse(hgvs):
            var = self.hp.parse_hgvs_variant(hgvs)
            del_re = re.compile('^del')
            if len(del_re.findall(str(var.posedit.edit))) == 1:
                del_length = len(var.posedit.edit.ref)
                if var.type == 'g':
                    if (var.posedit.pos.end.base - var.posedit.pos.start.base + 1) == del_length:
                        return True
                    else:
                        return False
                if var.type in ['c', 'r', 'p']:
                    if ((var.posedit.pos.end.base + var.posedit.pos.end.offset) -
                            (var.posedit.pos.start.base + var.posedit.pos.start.offset) + 1) == del_length:
                        return True
                    else:
                        return False
            else:
                return False


class ExtrinsicValidation():
    """
    Attempts to determine if the HGVS name validates against external data sources
    """
    def __init__(self):
        self.hp = hgvs.parser.Parser()
        self.ivalid = IntrinsicValidation()

    def valid_ac(self, hgvs):
        if self.ivalid.valid_parse(hgvs):
            var = self.hp.parse_hgvs_variant(hgvs)
            r = requests.get('http://uta.locusdev.net/api/v0/transcripts/fetch_seq?ac={ac}&start=0&end=10'
                            .format(ac=var.ac))
            if r.status_code == 200:
                return True
            else:
                return False
        else:
            return False

    def valid_ref(self, hgvs):
        if self.ivalid.valid_parse(hgvs):
            var = self.hp.parse_hgvs_variant(hgvs)
            if len(var.posedit.edit.ref) >= 1:
                r = requests.get('http://uta.locusdev.net/api/v0/transcripts/fetch_seq?ac={ac}&start={start}&end={end}'
                                .format(ac=var.ac, start=var.posedit.pos.start.base, end=var.posedit.pos.end.base))
                if r.status_code == 200:
                    return True
                else:
                    return False



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