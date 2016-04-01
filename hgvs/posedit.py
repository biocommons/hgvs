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


    def length_change(self, on_error_raise=True):
        """Returns the net length change for this posedit.

        The method for computing the net length change depends on the
        type of variant (dup, del, ins, etc).  The length_change
        method hides this complexity from callers.

        :param hgvs.posedit.PosEdit self: a PosEdit instance
        :param bool on_error_raise: whether to raise an exception on errors 

        :returns: A signed int for the net change in length.  Negative \
        values imply net deletions, 0 implies a balanced insertion and \
        deletion (e.g., SNV), and positive values imply a net \
        insertion.

        :raises HGVSUnsupportedOperationError: When determining the \
        length for this variant type is ill-defined or unsupported.

        There are many circumstances in which the net length change
        cannot be determined, is ill-defined, or is unsupported.  In
        these cases, the result depends on the value of
        `on_error_raise`: when `on_error_raise` is True, an exception
        is raised; when False, the exception is caught and `None` is
        returned.  Callers might wish to pass `on_error_raise=False`
        in list comprehensions to avoid dealing with exceptions.

        """

        try:
            ilen = self.pos._length()
            (del_len, ins_len) = self.edit._del_ins_lengths(ilen)
            return ins_len - del_len
        except HGVSUnsupportedOperationError:
            if on_error_raise:
                raise
            return None



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
