"""Deprecated: import from ``hgvs.parsers.pyparsing_grammar`` instead."""

import warnings

from hgvs.parsers.pyparsing_grammar import HGVSGrammar

warnings.warn(
    "hgvs.grammar is deprecated; import HGVSGrammar from hgvs.parsers.pyparsing_grammar instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["HGVSGrammar"]
