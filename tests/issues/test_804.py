"""https://github.com/biocommons/hgvs/issues/804

HGVS spec 21.1.3: an N-terminal extension is "Met1ext-5" -- aaterm must not
be embedded in the description (it supersedes the legacy "Met1extMet-5" form).
"""

import pytest

import hgvs.edit
import hgvs.location
import hgvs.posedit
from hgvs.parsers import Parser


def test_aaext_issue804():
    n_term = hgvs.posedit.PosEdit(
        pos=hgvs.location.Interval(
            start=hgvs.location.AAPosition(base=1, aa="M"),
            end=hgvs.location.AAPosition(base=1, aa="M"),
        ),
        edit=hgvs.edit.AAExt(ref="M", aaterm="LGM", length=-2),
    )
    assert str(n_term) == "Met1ext-2"

    # A C-terminal extension is "Ter110GlnextTer17": the aa replacing the
    # stop codon, then extTer, then the extension length.  A multi-residue
    # aaterm tail sequence must be reduced to those two slots.
    c_term = hgvs.posedit.PosEdit(
        pos=hgvs.location.Interval(
            start=hgvs.location.AAPosition(base=10, aa="*"),
            end=hgvs.location.AAPosition(base=10, aa="*"),
        ),
        edit=hgvs.edit.AAExt(ref="*", aaterm="LG*", length=2),
    )
    assert str(c_term) == "Ter10LeuextTer2"

    # spec-conformant constructions are unchanged
    assert str(hgvs.edit.AAExt(ref="", alt=None, aaterm="M", length=-5)) == "ext-5"
    assert str(hgvs.edit.AAExt(ref="", alt="Q", aaterm="*", length=17)) == "GlnextTer17"
    assert str(hgvs.edit.AAExt(ref="", alt="R", aaterm="*", length="?")) == "ArgextTer?"


# formatted pdot -> expected re-formatted pdot (identity when already
# spec-conformant; legacy forms normalize to the spec 21.1.3 shape)
_roundtrips = [
    ("Met1ext-2", "Met1ext-2"),
    ("Ter10LeuextTer2", "Ter10LeuextTer2"),
    ("Ter110GlnextTer17", "Ter110GlnextTer17"),
    ("Ter10extTer2", "Ter10extTer2"),
    # legacy single-residue N-terminal form: aaterm is dropped
    ("Met1extMet-5", "Met1ext-5"),
    # legacy multi-residue tails (as older formatters emitted):
    # reduced to the spec slots on format
    ("Met1extLeuGlyMet-2", "Met1ext-2"),
    ("Ter10GlnextLG*17", "Ter10GlnextTer17"),
    ("Ter10extLG*17", "Ter10LeuextTer17"),
]


@pytest.mark.parametrize("pdot,expected", _roundtrips)
def test_aaext_issue804_roundtrip(pdot, expected):
    # Full-circle: every pdot string the formatter emits must parse back
    # to a variant that re-formats identically.
    parser = Parser()
    var = parser.parse_hgvs_variant("NP_000000.0:p." + pdot)
    reformatted = str(var.posedit)
    assert reformatted == expected
    # second lap: the normalized form is a fixed point
    var2 = parser.parse_hgvs_variant("NP_000000.0:p." + reformatted)
    assert str(var2.posedit) == expected
