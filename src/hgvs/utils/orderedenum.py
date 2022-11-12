# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from enum import Enum


class OrderedEnum(Enum):
    def __ge__(self, other):
        assert self.__class__ is other.__class__, "OrderedEnum can only compare to OrderedEnum"
        return self.value >= other.value

    def __gt__(self, other):
        assert self.__class__ is other.__class__, "OrderedEnum can only compare to OrderedEnum"
        return self.value > other.value

    def __le__(self, other):
        assert self.__class__ is other.__class__, "OrderedEnum can only compare to OrderedEnum"
        return self.value <= other.value

    def __lt__(self, other):
        assert self.__class__ is other.__class__, "OrderedEnum can only compare to OrderedEnum"
        return self.value < other.value


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
