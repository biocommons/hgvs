
Manuscript Example
==================

.. code:: python

    import hgvs.parser
    hp = hgvs.parser.Parser()
    var_c1 = hp.parse_hgvs_variant('NM_182763.2:c.688+403C>T')
    var_c1, var_c1.posedit.pos.start, str(var_c1)



.. parsed-literal::

    (SequenceVariant(ac=NM_182763.2, type=c, posedit=688+403C>T),
     BaseOffsetPosition(base=688, offset=403, datum=1, uncertain=False),
     'NM_182763.2:c.688+403C>T')



.. code:: python

    import hgvs.dataproviders.uta
    import hgvs.hgvsmapper
    uta = hgvs.dataproviders.uta.connect(db_url='postgresql://localhost/uta')
    hm = hgvs.hgvsmapper.HGVSMapper(uta, cache_transcripts=True)
.. code:: python

    var_p1 = hm.c_to_p(var_c1)
    var_p1, str(var_p1)



.. parsed-literal::

    (SequenceVariant(ac=NP_877495.1, type=p, posedit=?), 'NP_877495.1:p.?')



.. code:: python

    var_g = hm.c_to_g(var_c1,'NC_000001.10','splign')
    var_g, str(var_g)



.. parsed-literal::

    (SequenceVariant(ac=NC_000001.10, type=g, posedit=150550916G>A),
     'NC_000001.10:g.150550916G>A')



.. code:: python

    var_c2 = hm.g_to_c(var_g,'NM_001197320.1','splign')
    var_c2, str(var_c2)



.. parsed-literal::

    (SequenceVariant(ac=NM_001197320.1, type=c, posedit=281C>T),
     'NM_001197320.1:c.281C>T')



.. code:: python

    var_p2 = hm.c_to_p(var_c2)
    var_p2, str(var_p2)



.. parsed-literal::

    (SequenceVariant(ac=NP_001184249.1, type=p, posedit=(Ser94Phe)),
     'NP_001184249.1:p.(Ser94Phe)')



