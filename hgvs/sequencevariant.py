# -*- coding: utf-8 -*-
"""represents simple sequence-based variants

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import recordtype
import re

import hgvs.variantmapper
from hgvs.utils.validationlevel import ValidationLevel


class SequenceVariant(recordtype.recordtype("SequenceVariant", ["ac", "type", "posedit"])):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """

    def format(self, conf=None):
        """Formatting the stringification of sequence variants

        :param conf: a dict comprises formatting options. None is to use global settings.

        See :class:`hgvs.config`.
        """
        if self.ac is not None:
            return "{ac}:{type}.{posedit}".format(ac=self.ac, type=self.type, posedit=self.posedit.format(conf))
        else:
            return "{type}.{posedit}".format(type=self.type, posedit=self.posedit.format(conf))

    __str__ = format

    def fill_ref(self, hdp):
        hm = hgvs.variantmapper.VariantMapper(hdp)
        type = self.posedit.edit.type
        if type in ["del", "delins", "identity", "dup", "inv"] and self.posedit.edit.ref_s is None:
            hm._replace_reference(self)
        if type == "identity" and isinstance(self.posedit.edit, hgvs.edit.NARefAlt):
            self.posedit.edit.alt = self.posedit.edit.ref
        return self

    def validate(self):
        if self.ac and self.type:
            (res, msg) = self._validate_ac_type_pair()
            if res != ValidationLevel.VALID:
                return (res, msg)
        return self.posedit.validate()

    def _validate_ac_type_pair(self):
        g_ac  = r'^(NC|NG|NT|NW|NZ|CM|LRG_\d+$)'
        cr_ac = r'^(NM|XM|ENST|LRG_\d+t\d+$)'
        n_ac  = r'^(NR|XR|ENST|LRG_\d+t\d+$)'
        p_ac  = r'^(NP|XP|ENSP|LRG_\d+p\d+$)'
        m_ac  = r'^(NC)'

        valid_pairs = {
            'g' : re.compile(g_ac),
            'c' : re.compile(cr_ac),
            'r' : re.compile(cr_ac),
            'n' : re.compile(n_ac),
            'p' : re.compile(p_ac),
            'm' : re.compile(m_ac)
        }

        invalid_pairs = {
            'g' : re.compile('|'.join((cr_ac, n_ac, p_ac))),
            'c' : re.compile('|'.join((g_ac, p_ac))),
            'r' : re.compile('|'.join((g_ac, p_ac))),
            'n' : re.compile('|'.join((g_ac, p_ac))),
            'p' : re.compile('|'.join((g_ac, cr_ac, n_ac))),
            'm' : re.compile('|'.join((cr_ac, n_ac, p_ac)))
        }

        if valid_pairs[self.type].match(self.ac):
            return (ValidationLevel.VALID, None)
        elif invalid_pairs[self.type].match(self.ac):
            return (ValidationLevel.ERROR,
                'Variant accession ({ac}) and type ({type}) do not match'.format(ac=self.ac, type=self.type))
        else:
            return (ValidationLevel.WARNING,
                'Variant accession ({ac}) and type ({type}) may not match'.format(ac=self.ac, type=self.type))


# <LICENSE>
# Copyright 2013-2015 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
