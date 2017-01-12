
Automated liftover of NM\_001261456.1:c.1762A>G (rs509749) to NM\_001261457.1 via GRCh37
========================================================================================

Automatically project variant from one transcript to another via
common reference.

http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=509749
                
.. code:: python

    import hgvs.parser
    hgvsparser = hgvs.parser.Parser()
    var_c1 = hgvsparser.parse_hgvs_variant('NM_001261456.1:c.1762A>G')

.. code:: python

    import hgvs.dataproviders.uta
    hdp = hgvs.dataproviders.uta.connect()

.. code:: python

    import hgvs.projector
    pj = hgvs.projector.Projector(hdp=hdp,
                                  alt_ac='NC_000001.10',
                                  src_ac=var_c1.ac,
                                  dst_ac='NM_001261457.1')

.. code:: python

    pj.project_variant_forward(var_c1)



.. parsed-literal::

    SequenceVariant(ac=NM_001261457.1, type=c, posedit=1534A>G)



