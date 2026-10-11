import pytest

@pytest.mark.parametrize("var_c_str,var_p_str", [
    ('NM_000159.4:c.1074_1075insAGTTGAAGGACC', 'NP_000150.1:p.(Asp358_Gln359insSerTer)'),
    ('NM_000159.4:c.1074_1075insAGTTGAAGGAC', 'NP_000150.1:p.(Asp358_Gln359insSerTer)'),
    ('NM_004380.3:c.138_139insTCATCATGAGCTCCC', 'NP_004371.2:p.(Pro46_Asn47insSerSerTer)'),
    ('NM_000038.6:c.4189_4197delinsATATAAAAAA', 'NP_000029.2:p.(Glu1397IlefsTer2)'),
    ('NM_004380.2:c.139delinsTCATCATGAGCTG', 'NP_004371.2:p.(Asn47delinsSerSerTer)'),
    ('NM_004448.4:c.2339dup', 'NP_004439.2:p.(Tyr781IlefsTer93)'),
    ('NM_024529.5:c.1595_*2delinsCCCC', 'NP_078805.3:p.(Ter532SerextTer29)'),
    ('NM_004985.5:c.567_*1delinsTC', 'NP_004976.2:p.(Ter189TyrextTer8)'),
    ('NM_004985.5:c.567_*1insTAGAAA', 'NP_004976.2:p.?'),
    ('NM_000203.4:c.1960T>C', 'NP_000194.2:p.(Ter654ArgextTer?)'),
    ('NM_000203.4:c.1730del', 'NP_000194.2:p.(Cys577SerfsTer?)'),
    # Insertion point inside a codon (frame_offset 2)
    ('NM_000159.4:c.1073_1074insTAGTTGAAGGA', 'NP_000150.1:p.(Asp358_Gln359insSerTer)'),
    # Insertion splits codon 358 (frame_offset 1) and changes Asp to Val
    ('NM_000159.4:c.1072_1073insTTAGTTGAAGG', 'NP_000150.1:p.(Lys357_Asp358insValSerTer)'),
    # Same protein as above with a different untranslated base after the stop
    ('NM_000159.4:c.1072_1073insTTAGTTGAGGG', 'NP_000150.1:p.(Lys357_Asp358insValSerTer)'),
    ('NM_000159.4:c.1072_1073insTTAGTTGACGG', 'NP_000150.1:p.(Lys357_Asp358insValSerTer)'),
    # Stop ends one base before the end of the insertion (frame_offset 2)
    ('NM_000159.4:c.1073_1074insATAGG', 'NP_000150.1:p.(Lys357_Asp358insGluTer)'),
    # Codon-aligned dup that encodes a stop
    ('NM_000159.4:c.1061_1080dup', 'NP_000150.1:p.(Asp360_Lys361insAlaAlaTer)'),
    # Mid-codon dup (frame_offset 1) where the split Lys361 codon becomes Thr
    ('NM_000159.4:c.1062_1081dup', 'NP_000150.1:p.(Asp360_Lys361insThrAlaTer)'),
    # Stop rebuilt at the reference stop (NP_078805.3 is 532 residues with the stop)
    ('NM_024529.5:c.1593_1594insTGA', 'NP_078805.3:p.(Ter532Ter)'),
    ('NM_024529.5:c.1593_1594insTAA', 'NP_078805.3:p.(Ter532Ter)'),
    ('NM_024529.5:c.1593_1594insTAG', 'NP_078805.3:p.(Ter532Ter)'),
    # Stop rebuilt at the reference stop with untranslated bases after it
    ('NM_024529.5:c.1593_1594insTAGAAA', 'NP_078805.3:p.(Ter532Ter)'),
    # Insertion inside the stop codon rebuilds a stop in place (T|AAT|GA to TAA TGA)
    ('NM_024529.5:c.1594_1595insAAT', 'NP_078805.3:p.(Ter532Ter)'),
    # Dup of the stop codon, equivalent to c.1593_1594insTGA
    ('NM_024529.5:c.1594_1596dup', 'NP_078805.3:p.(Ter532Ter)'),
    # Stop retained by an insertion after the first stop base (already correct before this fix)
    ('NM_024529.5:c.1595_1596insAGT', 'NP_078805.3:p.(Ter532Ter)'),
    # Residue inserted before a rebuilt stop is an extension, not stop-retained
    ('NM_024529.5:c.1593_1594insGCTTGA', 'NP_078805.3:p.(Ter532AlaextTer2)'),
    # Insertion destroys the stop codon
    ('NM_024529.5:c.1595_1596insGAT', 'NP_078805.3:p.(Ter532TrpextTer30)'),
    # Stop ends exactly at the end of the insertion, giving a nonsense variant
    ('NM_000159.4:c.1073_1074insTTAG', 'NP_000150.1:p.(Gln359Ter)'),
])
def test_c_to_p(var_c_str, var_p_str, parser, am37):
    var_c = parser.parse(var_c_str)
    var_p = am37.c_to_p(var_c)
    assert str(var_p) == var_p_str


@pytest.mark.parametrize("var_c_str,var_p_str", [
    ('NM_000159.4:c.1074_1075insAGTTGAAGGAC', 'NP_000150.1:p.(D358_Q359insS*)')
])
def test_c_to_p_1_letter(var_c_str, var_p_str, parser, am37):
    var_c = parser.parse(var_c_str)
    var_p = am37.c_to_p(var_c)
    assert var_p.format(conf={"p_3_letter": False}) == var_p_str
