# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import recordtype

from hgvs.exceptions import HGVSError, HGVSUnsupportedOperationError


class PosEdit(recordtype.recordtype('PosEdit', [('pos', None), ('edit', None), ('uncertain', False)])):
    """
    represents a **simple** variant, consisting of a single position and edit pair
    """

    def __str__(self):
        rv = str(self.edit) if self.pos is None else '{self.pos}{self.edit}'.format(self=self)
        if self.uncertain:
            if self.edit in ['0', '']:
                rv = rv + '?'
            else:
                rv = '(' + rv + ')'
        return rv

    def _set_uncertain(self):
        """sets the uncertain flag to True; used primarily by the HGVS grammar

        :returns: self
        """
        self.uncertain = True
        return self


    def length_change(self, error_value=None):
        """returns net length change for this posedit

        For an interval length, ilen, determined from the position,
        calls _del_ins_lengths(ilen) on the Edit instance, which
        returns a (del_len, ins_len) tuple.  This method then returns
        ins_len - del_len as the net change.

        This Heimlich maneuver is necessary to accommodate the
        optionality of alleles. For example c.10_15del requires that
        we determine the del length from the interval, where as
        c.10_15delinsCAT requires information from the interval and
        from the edit.

        There are many circumstances in which the net length change
        cannot be determined or is ill-defined. In these cases, the
        result depends on the value of `error_value`. When
        `error_value` is None, an exception is raised; when not None,
        the exception is caught and the value is returned.

        """

        try:
            ilen = self.pos._length()
            (del_len, ins_len) = self.edit._del_ins_lengths(ilen)
        except HGVSUnsupportedOperationError:
            if error_value is not None:
                return error_value
            raise

        return ins_len - del_len


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
