# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import collections
import warnings


class deprecated(object):
    '''Decorator factory class which returns a decorator function that
    marks a function as deprecated. It will result in a warning being
    emitted when the function is used, once per invocation point and
    deprecated function. The doc string of the deprecated function is
    overwritten with a deprecation message.

    @deprecated(use_instead="c_to_p(...)')
    def hgvsc_to_hgvsp():

    results in warnings like this:

    /home/reece/projects/hgvs/hgvs/variantmapper.py:280: DeprecationWarning: Call to deprecated function hgvsc_to_hgvsp; use c_to_p(...) instead

    '''

    def __init__(self, use_instead=None):
        # only warn once per caller location and function
        self.seen = collections.Counter()
        self.use_instead = use_instead

    def __call__(self, func):
        msg = "Call to deprecated function {}".format(func.__name__)
        if self.use_instead:
            msg += '; use ' + self.use_instead + ' instead'

        def wrapper(*args, **kwargs):
            fingerprint = (func.__name__, func.func_code.co_filename, func.func_code.co_firstlineno)
            if fingerprint not in self.seen:
                warnings.warn_explicit(
                    msg,
                    category=DeprecationWarning,
                    filename=func.func_code.co_filename,
                    lineno=func.func_code.co_firstlineno + 1)
            self.seen.update([fingerprint])
            return func(*args, **kwargs)

        wrapper.__doc__ = "Deprecated"
        if self.use_instead:
            wrapper.__doc__ += '; use ' + self.use_instead + ' instead'
        return wrapper
