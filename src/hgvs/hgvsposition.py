# -*- coding: utf-8 -*-
"""Represent partial HGVS tags that refer to a position without alleles

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import attr


@attr.s(slots=True, repr=False)
class HGVSPosition(object):
    """
    HGVSPosition -- Represent partial HGVS tags that refer to a position without alleles

    :param str ac: sequence accession
    :param str type: type of sequence and coordinate
    :param str pos: sequence position
    :param str gene: gene symbol (may be None)

    """
    ac = attr.ib()
    type = attr.ib()
    pos = attr.ib()
    gene = attr.ib(default=None)

    def __str__(self):
        g = "" if not self.gene else "(" + self.gene + ")"
        return "{self.ac}{g}:{self.type}.{self.pos}".format(self=self, g=g)

    def __repr__(self):
        return "{0}({1})".format(
            self.__class__.__name__, ", ".join(
                (a.name + "=" + str(getattr(self, a.name))) for a in self.__attrs_attrs__))


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
