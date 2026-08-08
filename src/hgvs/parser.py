"""Deprecated: import from ``hgvs.parsers`` instead.

``Parser`` moved to ``hgvs.parsers`` so that each grammar backend (OMeta/
parsley, pyparsing) lives in its own module with its own isolated
dependencies. This module re-exports the same names for backward
compatibility and will be removed in a future major release.
"""

import warnings

from hgvs.parsers import OMETA_GRAMMAR, PYPARSING_GRAMMAR, Parser

warnings.warn(
    "hgvs.parser is deprecated; import Parser from hgvs.parsers instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["OMETA_GRAMMAR", "PYPARSING_GRAMMAR", "Parser"]
