
http://www.ncbi.nlm.nih.gov/projects/SNP/snp\_ref.cgi?rs=509749

.. code:: python

    import bdi.sources.uta0_pg
    import hgvs.hgvsmapper
    import hgvs.parser
.. code:: python

    bdi = bdi.sources.uta0_pg.UTA0()
    hgvsmapper = hgvs.hgvsmapper.HGVSMapper(bdi)
    hgvsparser = hgvs.parser.Parser()
.. code:: python

    var_c1 = hgvsparser.parse_hgvs_variant('NM_001261456.1:c.1762A>G')
    var_p1 = hgvsmapper.hgvsc_to_hgvsp(var_c1, None)
    var_c1, var_p1



.. parsed-literal::

    (Variant(ac=NM_001261456.1, type=c, posedit=1762A>G),
     Variant(ac=NP_001248385.1, type=p, posedit=Met588Val))



.. code:: python

    var_g = hgvsmapper.hgvsc_to_hgvsg(var_c1,'GRCh37.p10')
    var_g



.. parsed-literal::

    Variant(ac=NC_000001.10, type=g, posedit=160793560A>G)



.. code:: python

    txs = bdi.get_tx_for_gene('LY9')
    len(txs)
    [ tx['ac'] for tx in txs ] 



.. parsed-literal::

    ['NM_002348.3', 'NM_001261456.1', 'NM_001261457.1', 'NM_001033667.2']



.. code:: python

    var_c2 = hgvsmapper.hgvsg_to_hgvsc(var_g,'NM_001261457.1')
    var_p2 = hgvsmapper.hgvsc_to_hgvsp(var_c2, None)
    var_c2, var_p2



.. parsed-literal::

    (Variant(ac=NM_001261457.1, type=c, posedit=1534A>G),
     Variant(ac=NP_001248386.1, type=p, posedit=Met512Val))



.. code:: python

    