def test_hgvs_to_vcf(parser, babelfish38):
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """
    
    def _h2v (h):
        return babelfish38.hgvs_to_vcf(parser.parse(h))

    # no-op
    assert _h2v("NC_000006.12:g.49949407=") == None

    # snv
    assert _h2v("NC_000006.12:g.49949407A>T") == ('6', 49949406, 'AA', 'AT', 'sub')

    # delins
    assert _h2v("NC_000006.12:g.49949413_49949414delinsCC") == ('6', 49949412, 'AAA', 'ACC', 'delins')
    
    # del, no shift
    assert _h2v("NC_000006.12:g.49949415del") == ('6', 49949414, 'AT', 'A', 'del')

    # del, w/ shift
    assert _h2v("NC_000006.12:g.49949413del") == ('6', 49949409, 'GA', 'G', 'del')
    assert _h2v("NC_000006.12:g.49949414del") == ('6', 49949409, 'GA', 'G', 'del')
    assert _h2v("NC_000006.12:g.49949413_49949414del") == ('6', 49949409, 'GAA', 'G', 'del')

    # ins, no shift
    assert _h2v("NC_000006.12:g.49949413_49949414insC") == ('6', 49949413, 'A', 'AC', 'ins')
    assert _h2v("NC_000006.12:g.49949414_49949415insCC") == ('6', 49949414, 'A', 'ACC', 'ins')

    # ins/dup, w/shift
    assert _h2v("NC_000006.12:g.49949413_49949414insA") == ('6', 49949409, 'G', 'GA', 'dup')
    assert _h2v("NC_000006.12:g.49949414_49949415insA") == ('6', 49949409, 'G', 'GA', 'dup')
    assert _h2v("NC_000006.12:g.49949414_49949415insAA") == ('6', 49949409, 'G', 'GAA', 'dup')
