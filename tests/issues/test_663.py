import pytest

from bioutils.sequences import TranslationTable

testdata = [
    # snv
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.586G>A',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.10G>A',
            'var_p': None,
        },
        id='NC_012920.1:m.586G>A'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.583G>A',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.7G>A',
            'var_p': None,
        },
        id='NC_012920.1:m.583G>A'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.579T>C',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.3T>C',
            'var_p': None,
        },
        id='NC_012920.1:m.579T>C'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.578T>C',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.2T>C',
            'var_p': None,
        },
        id='NC_012920.1:m.578T>C'
    ),
    # inv
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.15326_15327inv',
            'tx_ac': 'NC_012920.1_14746_15887',
            'var_t': 'NC_012920.1_14746_15887:c.580_581inv',
            'var_p': 'YP_003024038.1:p.(Thr194Val)',
        },
        id='NC_012920.1:m.15326_15327inv'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.8993_8994inv',
            'tx_ac': 'NC_012920.1_08526_09207',
            'var_t': 'NC_012920.1_08526_09207:c.467_468inv',
            'var_p': 'YP_003024031.1:p.(Leu156Pro)',
        },
        id='NC_012920.1:m.8993_8994inv'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.7146_7147inv',
            'tx_ac': 'NC_012920.1_05903_07445',
            'var_t': 'NC_012920.1_05903_07445:c.1243_1244inv',
            'var_p': 'YP_003024028.1:p.(Thr415Val)',
        },
        id='NC_012920.1:m.7146_7147inv'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.3902_3908inv',
            'tx_ac': 'NC_012920.1_03306_04262',
            'var_t': 'NC_012920.1_03306_04262:c.596_602inv',
            'var_p': 'YP_003024026.1:p.(Asp199_Ala201delinsGlyLysVal)',
        },
        id='NC_012920.1:m.3902_3908inv'
    ),
    # ins
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.5537_5538insT',
            'tx_ac': 'NC_012920.1_05511_05579',
            'var_t': 'NC_012920.1_05511_05579:n.26_27insT',
            'var_p': None,
        },
        id='NC_012920.1:m.5537_5538insT'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.5536_5537insT',
            'tx_ac': 'NC_012920.1_05511_05579',
            'var_t': 'NC_012920.1_05511_05579:n.25_26insT',
            'var_p': None,
        },
        id='NC_012920.1:m.5536_5537insT'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.4369_4370insA',
            'tx_ac': 'NC_012920.1_04328_04400',
            'var_t': 'NC_012920.1_04328_04400:n.35dup',
            'var_p': None,
        },
        id='NC_012920.1:m.4369_4370insA'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.1619_1620insT',
            'tx_ac': 'NC_012920.1_01601_01670',
            'var_t': 'NC_012920.1_01601_01670:n.18_19insT',
            'var_p': None,
        },
        id='NC_012920.1:m.1619_1620insT'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.960_961insC',
            'tx_ac': 'NC_012920.1_00647_01601',
            'var_t': 'NC_012920.1_00647_01601:n.313dup',
            'var_p': None,
        },
        id='NC_012920.1:m.960_961insC'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.597_598insT',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.21_22insT',
            'var_p': None,
        },
        id='NC_012920.1:m.597_598insT'
    ),
    # eq
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.1438=',
            'tx_ac': 'NC_012920.1_00647_01601',
            'var_t': 'NC_012920.1_00647_01601:n.791=',
            'var_p': None,
        },
        id='NC_012920.1:m.1438='
    ),
    # delins
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.15452_15453delinsAC',
            'tx_ac': 'NC_012920.1_14746_15887',
            'var_t': 'NC_012920.1_14746_15887:c.706_707delinsAC',
            'var_p': 'YP_003024038.1:p.(Leu236Thr)',
        },
        id='NC_012920.1:m.15452_15453delinsAC'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.15257_15258delinsAG',
            'tx_ac': 'NC_012920.1_14746_15887',
            'var_t': 'NC_012920.1_14746_15887:c.511_512delinsAG',
            'var_p': 'YP_003024038.1:p.(Asp171Ser)',
        },
        id='NC_012920.1:m.15257_15258delinsAG'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.4768_4769delinsCG',
            'tx_ac': 'NC_012920.1_04469_05511',
            'var_t': 'NC_012920.1_04469_05511:c.299_300delinsCG',
            'var_p': 'YP_003024027.1:p.(Met100Thr)',
        },
        id='NC_012920.1:m.4768_4769delinsCG'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.3308delinsAC',
            'tx_ac': 'NC_012920.1_03306_04262',
            'var_t': 'NC_012920.1_03306_04262:c.2delinsAC',
            'var_p': 'YP_003024026.1:p.Met1?',
        },
        id='NC_012920.1:m.3308delinsAC'
    ),
    # del
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.1656del',
            'tx_ac': 'NC_012920.1_01601_01670',
            'var_t': 'NC_012920.1_01601_01670:n.55del',
            'var_p': None,
        },
        id='NC_012920.1:m.1656del'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.961_962del',
            'tx_ac': 'NC_012920.1_00647_01601',
            'var_t': 'NC_012920.1_00647_01601:n.314_315del',
            'var_p': None,
        },
        id='NC_012920.1:m.961_962del'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.960delC',
            'tx_ac': 'NC_012920.1_00647_01601',
            'var_t': 'NC_012920.1_00647_01601:n.313del',
            'var_p': None,
        },
        id='NC_012920.1:m.960delC'
    ),
    # dup
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.958_960dupCCC',
            'tx_ac': 'NC_012920.1_00647_01601',
            'var_t': 'NC_012920.1_00647_01601:n.311_313dup',
            'var_p': None,
        },
        id='NC_012920.1:m.958_960dupCCC'
    ),
    pytest.param(
        {
            'var_g': 'NC_012920.1:m.595dup',
            'tx_ac': 'NC_012920.1_00576_00647',
            'var_t': 'NC_012920.1_00576_00647:n.19dup',
            'var_p': None,
        },
        id='NC_012920.1:m.595dup'
    ),
]

@pytest.mark.parametrize('testcase', testdata)
def test_g_to_t(testcase, parser, am38):
    var_m = parser.parse(testcase['var_g'])
    var_t = am38.g_to_t(var_m, testcase['tx_ac'])
    assert str(var_t) == testcase['var_t']
    if var_t.type == 'c':
        var_p = am38.c_to_p(var_t, translation_table=TranslationTable.vertebrate_mitochondrial)
        assert str(var_p) == testcase['var_p']
