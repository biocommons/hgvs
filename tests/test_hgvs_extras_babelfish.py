def test_hgvs_to_vcf(parser, babelfish38):
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """
    
    def _h2v (h):
        return babelfish38.hgvs_to_vcf(parser.parse(h))

    # no-op
    assert None == _h2v("NC_000006.12:g.49949407=")

    # snv
    assert ('6', 49949406, 'AA', 'AT', 'sub') == _h2v("NC_000006.12:g.49949407A>T")

    # delins
    assert ('6', 49949412, 'AAA', 'ACC', 'delins') == _h2v("NC_000006.12:g.49949413_49949414delinsCC")
    
    # del, no shift
    assert ('6', 49949414, 'AT', 'A', 'del') == _h2v("NC_000006.12:g.49949415del")

    # del, w/ shift
    assert ('6', 49949409, 'GA', 'G', 'del') == _h2v("NC_000006.12:g.49949413del")
    assert ('6', 49949409, 'GA', 'G', 'del') == _h2v("NC_000006.12:g.49949414del")
    assert ('6', 49949409, 'GAA', 'G', 'del') == _h2v("NC_000006.12:g.49949413_49949414del")

    # ins, no shift
    assert ('6', 49949413, 'A', 'AC', 'ins') == _h2v("NC_000006.12:g.49949413_49949414insC")
    assert ('6', 49949414, 'A', 'ACC', 'ins') == _h2v("NC_000006.12:g.49949414_49949415insCC")

    # ins/dup, w/shift
    assert ('6', 49949409, 'G', 'GA', 'dup') == _h2v("NC_000006.12:g.49949413_49949414insA")
    assert ('6', 49949409, 'G', 'GA', 'dup') == _h2v("NC_000006.12:g.49949414_49949415insA")
    assert ('6', 49949409, 'G', 'GAA', 'dup') == _h2v("NC_000006.12:g.49949414_49949415insAA")
