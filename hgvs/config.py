"""Defines a default set of config for use throughout the hgvs
package

Config are read from an ini-format file using configparser. A thin
wrapper on the ConfigParser instance provides *attribute* based
lookups (rather than key), and returns heuristically typed values
(e.g., "True" becomes True).

Although keys are settable, they are stringified on setting and
type-inferred on getting, which means that round-tripping works only
for str, int, and boolean.

>>> import hgvs.config

.. data:: hgvs.config.global_config

  global config, initialized once with defaults

"""

from configparser import ConfigParser, ExtendedInterpolation
from copy import copy
import logging
from pkg_resources import resource_stream
import re


logger = logging.getLogger(__name__)


class Config(object):
    def __init__(self, extended_interpolation=True):
        if extended_interpolation:
            cp = ConfigParser(interpolation=ExtendedInterpolation())
        else:
            cp = ConfigParser()
        cp.optionxform = _name_xform
        self._cp = cp

    def __dir__(self):
        return self._cp.keys()

    def __getattr__(self, k):
        return ConfigGroup(self._cp[k])

    __getitem__ = __getattr__

    def _read_file(self, flo):
        self._cp.read_file(flo)


class ConfigGroup(object):
    def __init__(self, section):
        self.__dict__['_section'] = section

    def __dir__(self):
        return self.__dict__['_section'].keys()

    def __getattr__(self, k):
        return _val_xform(self.__dict__['_section'][k])

    __getitem__ = __getattr__

    def __setattr__(self, k, v):
        logger.info(str(self.__class__.__name__) + '.__setattr__({k}, ...)'.format(k=k))
        self.__dict__['_section'][k] = str(v)

    __setitem__ = __setattr__


def _name_xform(o):
    """transform names to lowercase, without symbols (except underscore)
    Any chars other than alphanumeric are converted to an underscore
    """
    return re.sub('\W', '_', o.lower())


def _val_xform(v):
    if v == 'True':
        return True
    if v == 'False':
        return False
    try:
        return int(v)
    except ValueError:
        pass
    return v


_default_config = Config()
_default_config._read_file(resource_stream(__name__, '_data/defaults.ini'))

global_config = copy(_default_config)

