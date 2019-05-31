# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import pickle


class PersistentDict(dict):
    ''' Persistent dictionary 
    '''

    def __init__(self, filename, flag='c', *args, **kwds):
        self.filename = filename
        self.flag = flag    # r=readonly, c=create,write,read
        try:
            with open(self.filename, 'rb') as f:
                self.update(pickle.load(f))    # noqa: B301,BAN-B301
        except IOError:
            if self.flag == 'r':
                raise IOError('Cannot open file ' + self.filename)
        dict.__init__(self, *args, **kwds)

    def sync(self):
        if self.flag == 'r':
            return
        with open(self.filename, 'wb') as f:
            pickle.dump(dict(self), f, -1)

    def close(self):
        self.sync()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


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
