# -*- coding: utf-8 -*-
"""implements a (position,edit) tuple that represents a localized sequence change

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import attr

from hgvs.exceptions import HGVSUnsupportedOperationError
from hgvs.enums import ValidationLevel


@attr.s(slots=True, repr=False)
class PosEdit(object):
    """
    represents a **simple** variant, consisting of a single position and edit pair
    """
    pos = attr.ib(default=None)
    edit = attr.ib(default=None)
    uncertain = attr.ib(default=False)

    def format(self, conf=None):
        """Formatting the string of PosEdit
        """
        if self.pos is None:
            rv = str(self.edit.format(conf))
        else:
            rv = "{pos}{edit}".format(pos=self.pos.format(conf), edit=self.edit.format(conf))

        if self.uncertain:
            if self.edit in ["0", ""]:
                rv = rv + "?"
            else:
                rv = "(" + rv + ")"
        return rv

    __str__ = format

    def __repr__(self):
        return "{0}({1})".format(
            self.__class__.__name__, ", ".join(
                (a.name + "=" + str(getattr(self, a.name))) for a in self.__attrs_attrs__))

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

    def validate(self):
        if self.pos:
            (res, msg) = self.pos.validate()
            if res != ValidationLevel.VALID:
                return (res, msg)
            try:
                # Check ins length is 1
                if self.edit.type == "ins" and self.pos.end - self.pos.start != 1:
                    return (ValidationLevel.ERROR, "insertion length must be 1")
                # Check del length
                if self.edit.type in ["del", "delins"]:
                    ref_len = self.edit.ref_n
                    if ref_len is not None and ref_len != self.pos.end - self.pos.start + 1:
                        return (ValidationLevel.ERROR,
                                "Length implied by coordinates must equal sequence deletion length")
            except HGVSUnsupportedOperationError as err:
                return (ValidationLevel.WARNING, str(err))
        return (ValidationLevel.VALID, None)


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
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
