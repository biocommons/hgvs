import inspect
import logging

import six


_logger = logging.getLogger(__name__)


class LazyWrapper:
    """Wrap a function call with a wrapper that delays loading until the
    object is used.  "Used" here means any attribute lookup on the
    object whatsoever.

    """

    def __init__(self, load):
        self._lazywrapper_load_function = load
        self._lazywrapper_value = None

    def __dir__(self):
        """expose methods of underlying object, not wrapper

        NOTE: Calling __dir__ will trigger loading if necessary.
        Therefore, code that introspects on the wrapper (e.g., a
        debugger) might trigger loading.

        """
        return self._get().__dir__()

    def __repr__(self):
        return "<{}({})>".format(
            self.__class__.__name__,
            self._lazywrapper_value or self._lazywrapper_load_function.__name__)

    def __getattr__(self, attr):
        """Python calls __getattr__ only when an attribute couldn't be found
        through normal mechanisms.  For LazyWrapper, this is our cue
        that we need to pass the request to the shadow object

        """
        shadow = self._get()
        return shadow.__getattribute__(attr)


    def _get(self):
        """load object if necessary, cache it, and return it"""
        if not self._lazywrapper_value:
            src_fn = inspect.getsourcefile(self._lazywrapper_load_function)
            sl = inspect.getsourcelines(self._lazywrapper_load_function)
            src, src_ln = " ".join(sl[0]).strip(), sl[1]
            _logger.debug("Lazy loading `{src}` ({src_fn}:{src_ln})".format(
                src=src, src_fn=src_fn, src_ln=src_ln))
            if not six.PY2:
                tb = []
                for sf in inspect.stack()[2:]:
                    if sf.filename.startswith("<"):
                        break
                    tb += ["{sf.filename}:{sf.lineno}".format(sf=sf)]
                _logger.debug("\n  ".join(["Call path:"] + tb))
            self._lazywrapper_value = self._lazywrapper_load_function()
        return self._lazywrapper_value
