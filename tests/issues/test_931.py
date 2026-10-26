"""Insertions and duplications at the CDS/3'UTR boundary (#931).

The protein description depends on where the inserted bases land relative
to the stop codon:

* An ``ins`` or ``dup`` whose interval reaches into the 3'UTR (``c.N_*1ins``,
  ``c.N_*1dup``, ``c.*1dup``) puts the new bases after the stop codon. The
  coding sequence is untouched, so it is classified as a 3'UTR variant and
  reported as ``p.?`` like every other UTR variant (see #715).
* An ``ins`` or ``dup`` whose interval ends on the last base of the stop codon
  (``c.Ndup``, ``c.N-1_Nins``) still touches the stop codon. It is reported as
  stop-retained (``p.(TerNTer)``) so the stop-codon involvement is not missed
  (see #895 / #896).

Both rules are checked on a forward-strand transcript (NM_024529.5, stop TGA at
c.1594_1596) and a reverse-strand transcript (NM_004985.5 / KRAS, stop TAA at
c.565_567). Both transcripts are already in the test cache.
"""

import pytest

BOUNDARY_CASES = [
    # --- forward strand: NM_024529.5, CDS ends at c.1596 (TGA), c.*1 = A ---
    # inserted bases land after the stop codon -> 3'UTR variant -> p.?
    ("NM_024529.5:c.1596_*1insC", "NP_078805.3:p.?"),
    ("NM_024529.5:c.1596_*1insTAA", "NP_078805.3:p.?"),
    ("NM_024529.5:c.1596_*1dup", "NP_078805.3:p.?"),
    ("NM_024529.5:c.1595_*1dup", "NP_078805.3:p.?"),
    ("NM_024529.5:c.*1dup", "NP_078805.3:p.?"),
    ("NM_024529.5:c.*1_*2insC", "NP_078805.3:p.?"),
    # interval ends on the last stop base -> stop codon involved -> stop retained
    ("NM_024529.5:c.1596dup", "NP_078805.3:p.(Ter532Ter)"),
    ("NM_024529.5:c.1595_1596insA", "NP_078805.3:p.(Ter532Ter)"),
    # --- reverse strand: NM_004985.5 (KRAS), CDS ends at c.567 (TAA), c.*1 = A ---
    ("NM_004985.5:c.567_*1insC", "NP_004976.2:p.?"),
    ("NM_004985.5:c.567_*1insTAA", "NP_004976.2:p.?"),
    ("NM_004985.5:c.567_*1dup", "NP_004976.2:p.?"),
    ("NM_004985.5:c.566_*1dup", "NP_004976.2:p.?"),
    ("NM_004985.5:c.*1dup", "NP_004976.2:p.?"),
    ("NM_004985.5:c.*1_*2insC", "NP_004976.2:p.?"),
    ("NM_004985.5:c.567dup", "NP_004976.2:p.(Ter189Ter)"),
    ("NM_004985.5:c.566_567dup", "NP_004976.2:p.(Ter189Ter)"),
    ("NM_004985.5:c.565_567dup", "NP_004976.2:p.(Ter189Ter)"),
    ("NM_004985.5:c.566_567insA", "NP_004976.2:p.(Ter189Ter)"),
]


@pytest.mark.parametrize(("var_c_str", "var_p_str"), BOUNDARY_CASES)
def test_c_to_p_at_cds_utr_boundary(var_c_str, var_p_str, parser, am37):
    var_c = parser.parse(var_c_str)
    var_p = am37.c_to_p(var_c)
    assert str(var_p) == var_p_str


@pytest.mark.parametrize(
    ("var_c_str", "var_p_str"),
    [
        ("NM_024529.5:c.1596dup", "NP_078805.3:p.(*532*)"),
        ("NM_024529.5:c.1596_*1insC", "NP_078805.3:p.?"),
        ("NM_004985.5:c.567dup", "NP_004976.2:p.(*189*)"),
        ("NM_004985.5:c.567_*1insC", "NP_004976.2:p.?"),
    ],
)
def test_c_to_p_at_cds_utr_boundary_1_letter(var_c_str, var_p_str, parser, am37):
    var_c = parser.parse(var_c_str)
    var_p = am37.c_to_p(var_c)
    assert var_p.format(conf={"p_3_letter": False}) == var_p_str
