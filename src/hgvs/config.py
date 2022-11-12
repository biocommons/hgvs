# -*- coding: utf-8 -*-
"""The hgvs package uses a single, package-wide configuration instance
to control package behavior.  The hgvs.config module provides that
configuration instance, via the `hgvs.global_config` variable.

You should not import hgvs.config directly.

Config are read from an ini-format file.  `hgvs.config` implements a
thin wrapper on the ConfigParser instance in order to provide
*attribute* based lookups (rather than key). It also returns
heuristically typed values (e.g., "True" becomes True). 

Although keys are settable, they are stringified on setting and
type-inferred on getting, which means that round-tripping works only
for str, int, and boolean.

>>> import hgvs.config

.. data:: hgvs.config.global_config

   Package-wide ("global") configuration, initialized with package
   defaults.  Setting configuration in this object will change global
   behavior of the hgvs package.

   global_config, an instance of :ref:``hgvs.config.Config``, supports
   reading ini-like files that updates

"""

from __future__ import absolute_import, division, print_function, unicode_literals

from configparser import ConfigParser, ExtendedInterpolation
from copy import copy
import logging
from pkg_resources import resource_stream
import re

logger = logging.getLogger(__name__)


class Config(object):
    """provides an attribute-based lookup of configparser sections and
    settings.

    """

    def __init__(self, extended_interpolation=True):
        if extended_interpolation:
            cp = ConfigParser(interpolation=ExtendedInterpolation())
        else:
            cp = ConfigParser()
        cp.optionxform = _name_xform
        self._cp = cp

    def read_stream(self, flo):
        """read configuration from ini-formatted file-like object

        """
        self._cp.read_string(flo.read().decode('ascii'))

    def __copy__(self):
        new_config = Config.__new__(Config)
        new_config._cp = object.__getattribute__(self, '_cp')
        return new_config

    def __dir__(self):
        return list(self._cp.keys())

    def __getattr__(self, k):
        # Work around PyCharm bug https://youtrack.jetbrains.com/issue/PY-4213
        if k == "_cp":
            return
        try:
            return ConfigGroup(self._cp[k])
        except KeyError:
            raise AttributeError(k)

    __getitem__ = __getattr__


class ConfigGroup(object):
    def __init__(self, section):
        self.__dict__["_section"] = section

    def __dir__(self):
        return list(self.__dict__["_section"].keys())

    def __getattr__(self, k):
        return _val_xform(self.__dict__["_section"][k])

    __getitem__ = __getattr__

    def __setattr__(self, k, v):
        logger.info(str(self.__class__.__name__) + ".__setattr__({k}, ...)".format(k=k))
        self.__dict__["_section"][k] = str(v)

    __setitem__ = __setattr__


def _name_xform(o):
    """transform names to lowercase, without symbols (except underscore)
    Any chars other than alphanumeric are converted to an underscore
    """
    return re.sub(r"\W", "_", o.lower())


def _val_xform(v):
    if v == "True":
        return True
    if v == "False":
        return False
    if v == "None":
        return None
    try:
        return int(v)
    except ValueError:
        pass
    return v


_default_config = Config()
_default_config.read_stream(resource_stream(__name__, "_data/defaults.ini"))

global_config = copy(_default_config)
