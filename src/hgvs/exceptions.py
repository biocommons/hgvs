# -*- coding: utf-8 -*-
"""defined exceptions used by the hgvs package

"""

from __future__ import absolute_import, division, print_function, unicode_literals


class HGVSError(Exception):
    pass


class HGVSDataNotAvailableError(HGVSError):
    pass


class HGVSInternalError(HGVSError):
    pass


class HGVSInvalidIntervalError(HGVSError):
    pass


class HGVSInvalidVariantError(HGVSError):
    """Exception raised when variant is inconsistent or invalid"""
    pass


class HGVSNormalizationError(HGVSError):
    pass


class HGVSParseError(HGVSError):
    pass


class HGVSUnsupportedOperationError(HGVSError):
    pass


class HGVSUsageError(HGVSError):
    """Exception raised when client/caller has made an invalid request"""
    pass


class HGVSVerifyFailedError(HGVSError):
    """Exception raised when the cached hdp result is not consistent with latest remote result"""
    pass


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
