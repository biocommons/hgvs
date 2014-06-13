import collections
import functools
import warnings

# Make sure we're showing DeprecationWarnings
warnings.filterwarnings('default','',DeprecationWarning)

# only warn once per caller location and function
seen = collections.Counter()

def deprecated(func):
    '''This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.'''

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        fingerprint = (func.__name__,func.func_code.co_filename,func.func_code.co_firstlineno)
        if fingerprint not in seen:
            warnings.warn_explicit(
                "Call to deprecated function {}.".format(func.__name__),
                category=DeprecationWarning,
                filename=func.func_code.co_filename,
                lineno=func.func_code.co_firstlineno + 1)
        seen.update([fingerprint])
        return func(*args, **kwargs)
    return wrapper
