"""In-frame stop-loss variants that also change the residue before the stop (#941).

``NM_000098.2`` ends ``... AAA AGT TAA | CTT CTG GGC AGA TGA`` (Lys657 Ser658 Ter659,
then 3'UTR). ``c.1973_1975del`` removes ``GTT``: the Ser658 and Ter659 codons collapse
into ``AAA`` (Lys) and translation continues in frame through Leu Leu Gly Arg to a new
stop. The description names both stops, ``p.(Ser658_Ter659delinsLysLeuLeuGlyArgTer)``,
so that a consumer of the p. description can tell that the reference stop is lost.
The earlier ``p.(Ser658delinsLysLeuLeuGlyArg)`` describes the same protein but reads
like an ordinary in-frame delins.
"""

import pytest

STOP_LOSS_CASES = [
    # Ser658 changes and the stop is lost: both stops are named
    ("NM_000098.2:c.1973_1975del", "NP_000089.1:p.(Ser658_Ter659delinsLysLeuLeuGlyArgTer)"),
    ("NM_000098.2:c.1971_1976del", "NP_000089.1:p.(Ser658_Ter659delinsLeuLeuGlyArgTer)"),
    # only the stop codon changes: extension, as before
    ("NM_000098.2:c.1975_1977del", "NP_000089.1:p.(Ter659LeuextTer4)"),
    # true frameshift, as before
    ("NM_000098.2:c.1973_1974del", "NP_000089.1:p.(Ser658IlefsTer35)"),
    # in-frame deletion that leaves the stop codon alone, as before
    ("NM_000098.2:c.1969_1971del", "NP_000089.1:p.(Lys657del)"),
]


@pytest.mark.parametrize(("var_c_str", "var_p_str"), STOP_LOSS_CASES)
def test_c_to_p_in_frame_stop_loss(var_c_str, var_p_str, parser, vm):
    var_c = parser.parse(var_c_str)
    var_p = vm.c_to_p(var_c)
    assert str(var_p) == var_p_str


def test_c_to_p_in_frame_stop_loss_1_letter(parser, vm):
    var_p = vm.c_to_p(parser.parse("NM_000098.2:c.1973_1975del"))
    assert var_p.format(conf={"p_3_letter": False}) == "NP_000089.1:p.(S658_*659delinsKLLGR*)"
