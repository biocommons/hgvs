
Manual liftover of NM\_001261456.1:c.1762A>G (rs509749) to NM\_001261457.1 via GRCh37
=====================================================================================

                http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=509749
                
.. code:: python

    import hgvs.dataproviders.uta
    import hgvs.variantmapper
    import hgvs.parser
.. code:: python

    hdp = hgvs.dataproviders.uta.connect()
    variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    hgvsparser = hgvs.parser.Parser()
.. code:: python

    var_c1 = hgvsparser.parse_hgvs_variant('NM_001261456.1:c.1762A>G')
    var_p1 = variantmapper.c_to_p(var_c1, None)
    var_c1, var_p1



.. parsed-literal::

    (SequenceVariant(ac=NM_001261456.1, type=c, posedit=1762A>G),
     SequenceVariant(ac=MD5_e999a940ca422ec8cab9bc3cc64e0d7d, type=p, posedit=(Met588Val)))



.. code:: python

    var_g = variantmapper.c_to_g(var_c1,'NC_000001.10')
    var_g



.. parsed-literal::

    SequenceVariant(ac=NC_000001.10, type=g, posedit=160793560A>G)



.. code:: python

    txs = hdp.get_tx_for_gene('LY9')
    txs



.. parsed-literal::

    [['LY9', 30, 1998, 'ENST00000263285', 'NC_000001.10', 'genebuild'],
     ['LY9', 1, 583, 'ENST00000368039', 'NC_000001.10', 'genebuild'],
     ['LY9', 0, 1648, 'ENST00000392203', 'NC_000001.10', 'genebuild'],
     ['LY9', 0, 1833, 'ENST00000368037', 'NC_000001.10', 'genebuild'],
     ['LY9', 211, 1024, 'ENST00000368035', 'NC_000001.10', 'genebuild'],
     ['LY9', 50, 1616, 'ENST00000341032', 'NC_000001.10', 'genebuild'],
     ['LY9', 170, 1751, 'ENST00000368041', 'NC_000001.10', 'genebuild'],
     ['LY9', 1094, 1907, 'ENST00000368040', 'NC_000001.10', 'genebuild'],
     ['LY9', 114, 2040, 'NM_001261456.1', 'AC_000133.1', 'splign'],
     ['LY9', 114, 2040, 'NM_001261456.1', 'NC_000001.10', 'blat'],
     ['LY9', 114, 2040, 'NM_001261456.1', 'NC_000001.10', 'splign'],
     ['LY9', 114, 2040, 'NM_001261456.1', 'NC_018912.2', 'splign'],
     ['LY9', 114, 696, 'NM_001033667.2', 'AC_000133.1', 'splign'],
     ['LY9', 114, 696, 'NM_001033667.2', 'NC_000001.10', 'blat'],
     ['LY9', 114, 696, 'NM_001033667.2', 'NC_000001.10', 'splign'],
     ['LY9', 114, 696, 'NM_001033667.2', 'NC_018912.2', 'splign'],
     ['LY9', 114, 2082, 'NM_002348.3', 'AC_000133.1', 'splign'],
     ['LY9', 114, 2082, 'NM_002348.3', 'NC_000001.10', 'blat'],
     ['LY9', 114, 2082, 'NM_002348.3', 'NC_000001.10', 'splign'],
     ['LY9', 114, 2082, 'NM_002348.3', 'NC_018912.2', 'splign'],
     ['LY9', 114, 1812, 'NM_001261457.1', 'AC_000133.1', 'splign'],
     ['LY9', 114, 1812, 'NM_001261457.1', 'NC_000001.10', 'blat'],
     ['LY9', 114, 1812, 'NM_001261457.1', 'NC_000001.10', 'splign'],
     ['LY9', 114, 1812, 'NM_001261457.1', 'NC_018912.2', 'splign']]



.. code:: python

    var_c2 = variantmapper.g_to_c(var_g,'NM_001261457.1',alt_aln_method='splign')
    var_p2 = variantmapper.c_to_p(var_c2, None)
    var_c2, var_p2



.. parsed-literal::

    (SequenceVariant(ac=NM_001261457.1, type=c, posedit=1534A>G),
     SequenceVariant(ac=MD5_921ebefe79bff479f4bfa17e133fc084, type=p, posedit=(Met512Val)))



