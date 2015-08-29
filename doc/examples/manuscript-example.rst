
Manuscript Example
==================

.. code:: python

    import hgvs
    hgvs.__version__



.. parsed-literal::

    '0.3dev-283858cb6466'



Parse an HGVS string into a Python structure
--------------------------------------------

.. code:: python

    import hgvs.parser 
    hp = hgvs.parser.Parser()
    var_c1 = hp.parse_hgvs_variant('NM_182763.2:c.688+403C>T')
    var_c1, var_c1.posedit.pos.start



.. parsed-literal::

    (SequenceVariant(ac=NM_182763.2, type=c, posedit=688+403C>T),
     BaseOffsetPosition(base=688, offset=403, datum=1, uncertain=False))



Open the UTA public data source for mapping and validation
----------------------------------------------------------

.. code:: python

    import hgvs.dataproviders.uta
    hdp = hgvs.dataproviders.uta.connect()
Project transcript variant NM\_182763.2:c.688+403C>T to GRCh37 primary assembly using splign alignments
-------------------------------------------------------------------------------------------------------

.. code:: python

    import hgvs.variantmapper
    vm = hgvs.variantmapper.EasyVariantMapper(
        hdp, primary_assembly='GRCh37', alt_aln_method='splign')
    var_g = vm.c_to_g(var_c1)
    var_g



.. parsed-literal::

    SequenceVariant(ac=NC_000001.10, type=g, posedit=150550916G>A)



Project genomic variant to a new transcript
-------------------------------------------

.. code:: python

    vm.relevant_transcripts(var_g)



.. parsed-literal::

    ['NM_182763.2', 'NM_021960.4', 'NM_001197320.1']



.. code:: python

    var_c2 = vm.g_to_c(var_g,'NM_001197320.1')
    var_c2



.. parsed-literal::

    SequenceVariant(ac=NM_001197320.1, type=c, posedit=281C>T)



Infer protein changes for these transcript variants
---------------------------------------------------

.. code:: python

    var_p1 = vm.c_to_p(var_c1)
    var_p2 = vm.c_to_p(var_c2)
    var_p1, var_p2



.. parsed-literal::

    (SequenceVariant(ac=NP_877495.1, type=p, posedit=?),
     SequenceVariant(ac=NP_001184249.1, type=p, posedit=(Ser94Phe)))



Format the results by "stringification"
---------------------------------------

.. code:: python

    print("""mapped {var_c1} ({var_p1})
        to {var_c2} ({var_p2})
       via {var_g}""".format(
            var_c1=var_c1, var_p1=var_p1,
            var_c2=var_c2, var_p2=var_p2,
            var_g=var_g))

.. parsed-literal::

    mapped NM_182763.2:c.688+403C>T (NP_877495.1:p.?)
        to NM_001197320.1:c.281C>T (NP_001184249.1:p.(Ser94Phe))
       via NC_000001.10:g.150550916G>A


Validate a variant
------------------

.. code:: python

    import hgvs.validator
    import hgvs.exceptions
    vr = hgvs.validator.Validator(hdp=hdp)
    try:
        vr.validate( hp.parse_hgvs_variant('NM_001197320.1:c.281C>T') )
        vr.validate( hp.parse_hgvs_variant('NM_001197320.1:c.281A>T') )
    except hgvs.exceptions.HGVSError as e:
        print(e)

.. parsed-literal::

    NM_001197320.1:c.281A>T: Variant reference does not agree with reference sequence

