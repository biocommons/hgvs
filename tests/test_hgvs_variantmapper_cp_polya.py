"""c_to_p for transcripts whose stop codon is completed by the poly(A) tail.

RefSeq annotates some transcripts whose CDS ends in a partial stop codon (T or TA)
that is completed by polyadenylation, so the CDS length is not a multiple of 3.
NM_001424184.1 (TMEM247) has a 659 nt CDS ending in ...GTCATA; the final TA plus
the first A of the poly(A) tail form the TAA stop. Previously c_to_p raised
NotImplementedError for such transcripts (ported from openvar/vv_hgvs).
"""

import pytest

from hgvs.utils.reftranscriptdata import RefTranscriptData


@pytest.mark.parametrize(
    ("cds_seq", "expected"),
    [
        ("ATGT", "AA"),
        ("ATGTA", "A"),
        ("ATGTAA", None),  # complete
        ("ATGC", None),  # cannot be completed to a stop
        ("ATGCA", None),
        ("ATGTAAC", None),
    ],
)
def test_poly_a_completion(cds_seq, expected):
    assert RefTranscriptData._poly_a_completion(cds_seq) == expected


@pytest.mark.parametrize(
    ("var_c_str", "var_p_str"),
    [
        ("NM_001424184.1:c.3G>A", "NP_001411113.1:p.Met1?"),
        ("NM_001424184.1:c.658T>C", "NP_001411113.1:p.(Ter220GlnextTer?)"),
        ("NM_001424184.1:c.659A>G", "NP_001411113.1:p.(Ter220Ter)"),
        # the stop codon is lost and nothing follows it, so the C-terminus is unknown
        ("NM_001424184.1:c.657_659del", "NP_001411113.1:p.?"),
    ],
)
def test_c_to_p_poly_a_completed_stop(var_c_str, var_p_str, parser, vm):
    var_c = parser.parse(var_c_str)
    assert str(vm.c_to_p(var_c)) == var_p_str
