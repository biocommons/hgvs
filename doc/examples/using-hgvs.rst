
Using hgvs
==========

This notebook demonstrates major features of the hgvs package.

.. code:: ipython2

    import hgvs
    hgvs.__version__




.. parsed-literal::

    '0.5.0a6.dev3+nf998c16a46b3.d20161012'



Variant I/O
-----------

Initialize the parser
~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    # You only need to do this once per process
    import hgvs.parser
    hp = hgvsparser = hgvs.parser.Parser()

Parse a simple variant
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    v = hp.parse_hgvs_variant("NC_000007.13:g.21726874G>A")

.. code:: ipython2

    v




.. parsed-literal::

    SequenceVariant(ac=NC_000007.13, type=g, posedit=21726874G>A)



.. code:: ipython2

    v.ac, v.type




.. parsed-literal::

    ('NC_000007.13', 'g')



.. code:: ipython2

    v.posedit




.. parsed-literal::

    PosEdit(pos=21726874, edit=G>A, uncertain=False)



.. code:: ipython2

    v.posedit.pos




.. parsed-literal::

    Interval(start=21726874, end=21726874, uncertain=False)



.. code:: ipython2

    v.posedit.pos.start




.. parsed-literal::

    SimplePosition(base=21726874, uncertain=False)



Parsing complex variants
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    v = hp.parse_hgvs_variant("NM_003777.3:c.13552_*36del57")

.. code:: ipython2

    v.posedit.pos.start, v.posedit.pos.end




.. parsed-literal::

    (BaseOffsetPosition(base=13552, offset=0, datum=1, uncertain=False),
     BaseOffsetPosition(base=36, offset=0, datum=2, uncertain=False))



.. code:: ipython2

    v.posedit.edit




.. parsed-literal::

    NARefAlt(ref=57, alt=None, uncertain=False)



Formatting variants
~~~~~~~~~~~~~~~~~~~

All objects may be formatted simply by “stringifying” or printing them
using ``str``, ``print()``, or ``"".format()``.

.. code:: ipython2

    str(v)




.. parsed-literal::

    'NM_003777.3:c.13552_*36del57'



.. code:: ipython2

    print(v)


.. parsed-literal::

    NM_003777.3:c.13552_*36del57


.. code:: ipython2

    "{v} spans the CDS end".format(v=v)




.. parsed-literal::

    'NM_003777.3:c.13552_*36del57 spans the CDS end'



Projecting variants between sequences
-------------------------------------

Set up a dataprovider
~~~~~~~~~~~~~~~~~~~~~

Mapping variants requires exon structures, alignments, CDS bounds, and
raw sequence. These are provided by a ``hgvs.dataprovider`` instance.
The only dataprovider provided with hgvs uses UTA. You may write your
own by subsclassing hgvs.dataproviders.interface.

.. code:: ipython2

    import hgvs.dataproviders.uta
    hdp = hgvs.dataproviders.uta.connect()

Initialize mapper classes
~~~~~~~~~~~~~~~~~~~~~~~~~

The VariantMapper class projects variants between two sequence
accessions using alignments from a specified source. In order to use it,
you must know that two sequences are aligned. VariantMapper isn’t
demonstrated here.

AssemblyMapper builds on VariantMapper and handles identifying
appropriate sequences. It is configured for a particular genome
assembly.

.. code:: ipython2

    import hgvs.variantmapper
    #vm = variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    am37 = easyvariantmapper = hgvs.variantmapper.AssemblyMapper(hdp, assembly_name='GRCh37')
    am38 = easyvariantmapper = hgvs.variantmapper.AssemblyMapper(hdp, assembly_name='GRCh38')

c_to_g
~~~~~~

This is the easiest case because there is typically only one alignment
between a transcript and the genome. (Exceptions exist for
pseudoautosomal regions.)

.. code:: ipython2

    var_c = hp.parse_hgvs_variant("NM_015120.4:c.35G>C")
    var_g = am37.c_to_g(var_c)
    var_g

.. code:: ipython2

    am38.c_to_g(var_c)

g_to_c
~~~~~~

In order to project a genomic variant onto a transcript, you must tell
the AssemblyMapper which transcript to use.

.. code:: ipython2

    am37.relevant_transcripts(var_g)




.. parsed-literal::

    ['NM_015120.4']



.. code:: ipython2

    am37.g_to_c(var_g, "NM_015120.4")




.. parsed-literal::

    SequenceVariant(ac=NM_015120.4, type=c, posedit=35T>C)



c_to_p
~~~~~~

.. code:: ipython2

    var_p = am37.c_to_p(var_c)
    str(var_p)




.. parsed-literal::

    'NP_055935.4:p.(Leu12Pro)'



.. code:: ipython2

    var_p.posedit.uncertain = False
    str(var_p)




.. parsed-literal::

    'NP_055935.4:p.Leu12Pro'



Projecting in the presence of a genome-transcript gap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As of Oct 2016, 1033 RefSeq transcripts in 433 genes have gapped
alignments. These gaps require special handlingin order to maintain the
correspondence of positions in an alignment. hgvs uses the precomputed
alignments in UTA to correctly project variants in exons containing
gapped alignments.

This example demonstrates projecting variants in the presence of a gap
in the alignment of NM_015120.4 (ALMS1) with GRCh37 chromosome 2. (The
alignment with GRCh38 is similarly gapped.) Specifically, the adjacent
genomic positions 73613031 and 73613032 correspond to the non-adjacent
CDS positions 35 and 39.

::

    NM_015120.4  c         15 >                           >       58
    NM_015120.4  n        126 > CCGGGCGAGCTGGAGGAGGAGGAG  >      169
                                |||||||||||   ||||||||||  21=3I20= 
    NC_000002.11 g   73613021 > CCGGGCGAGCT---GGAGGAGGAG  > 73613041
    NC_000002.11 g   73613021 < GGCCCGCTCGA---CCTCCTCCTC  < 73613041                                                  

.. code:: ipython2

    str(am37.c_to_g(hp.parse_hgvs_variant("NM_015120.4:c.35G>C")))




.. parsed-literal::

    'NC_000002.11:g.73613031T>C'



.. code:: ipython2

    str(am37.c_to_g(hp.parse_hgvs_variant("NM_015120.4:c.39G>C")))




.. parsed-literal::

    'NC_000002.11:g.73613032G>C'



Normalizing variants
--------------------

In hgvs, normalization means shifting variants 3’ (as requried by the
HGVS nomenclature) as well as rewriting variants. The variant
“NM_001166478.1:c.30_31insT” is in a poly-T run (on the transcript). It
should be shifted 3’ and is better written as dup, as shown below:

::

                                            *                       NC_000006.11:g.49917127dupA
      NC_000006.11 g   49917117 > AGAAAGAAAAATAAAACAAAG  > 49917137 
      NC_000006.11 g   49917117 < TCTTTCTTTTTATTTTGTTTC  < 49917137 
                                  |||||||||||||||||||||  21= 
    NM_001166478.1 n         41 < TCTTTCTTTTTATTTTGTTTC  <       21 NM_001166478.1:n.35dupT
    NM_001166478.1 c         41 <                        <       21 NM_001166478.1:c.30_31insT

.. code:: ipython2

    import hgvs.normalizer
    hn = hgvs.normalizer.Normalizer(hdp)

.. code:: ipython2

    v = hp.parse_hgvs_variant("NM_001166478.1:c.30_31insT")
    str(hn.normalize(v))




.. parsed-literal::

    'NM_001166478.1:c.35dupT'



A more complex normalization example
------------------------------------

This example is based on https://github.com/biocommons/hgvs/issues/382/.

::

      NC_000001.11 g   27552104 > CTTCACACGCATCCTGACCTTG > 27552125
      NC_000001.11 g   27552104 < GAAGTGTGCGTAGGACTGGAAC < 27552125
                                  |||||||||||||||||||||| 22= 
    NM_001029882.3 n        843 < GAAGTGTGCGTAGGACTGGAAC <      822 
    NM_001029882.3 c         12 <                        <      -10 
                                            ^^  
                                            NM_001029882.3:c.1_2del
                                            NM_001029882.3:n.832_833delAT
                                            NC_000001.11:g.27552114_27552115delAT

.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.1A>G"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552115T>C)



.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.2T>G"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552114A>C)



.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.1_2del"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552114_27552115delAT)



The genomic coordinates for the SNVs at c.1 and c.2 match those for the
del at c.1_2. Good!

Now, notice what happens with c.1_3del, c.1_4del, and c.1_5del:

.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.1_3del"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552114_27552116delATC)



.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.1_4del"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552112_27552115delGCAT)



.. code:: ipython2

    am38.c_to_g(hp.parse_hgvs_variant("NM_001029882.3:c.1_5del"))




.. parsed-literal::

    SequenceVariant(ac=NC_000001.11, type=g, posedit=27552112_27552116delGCATC)



Explanation:

On the transcript, c.1_2delAT deletes AT from …AGGATGCG…, resulting in
…AGGGCG…. There’s no ambiguity about what sequence was actually deleted.

c.1_3delATG deletes ATG, resulting in …AGGCG…. Note that you could also
get this result by deleting GAT. This is an example of an indel that is
subject to normalization and hgvs does this.

c.1_4delATGC and 1_5delATGCG have similar behaviors.

Normalization is always 3’ with respect to the reference sequence. So,
after projecting from a - strand transcript to the genome, normalization
will go in the opposite direction to the transcript. It will have
roughly the same effect as being 5’ shifted on the transcript (but
revcomp’d).

For more precise control, see the ``normalize`` and
``replace_reference`` options of ``AssemblyMapper``.

Validating variants
-------------------

``hgvs.validator.Validator`` is a composite of two classes,
``hgvs.validator.IntrinsicValidator`` and
``hgvs.validator.ExtrinsicValidator``. Intrinsic validation evaluates a
given variant for *internal* consistency, such as requiring that
insertions specify adjacent positions. Extrinsic validation evaluates a
variant using external data, such as ensuring that the reference
nucleotide in the variant matches that implied by the reference sequence
and position. Validation returns ``True`` if successful, and raises an
exception otherwise.

.. code:: ipython2

    import hgvs.validator
    hv = hgvs.validator.Validator(hdp)

.. code:: ipython2

    hv.validate(hp.parse_hgvs_variant("NM_001166478.1:c.30_31insT"))




.. parsed-literal::

    True



.. code:: ipython2

    from hgvs.exceptions import HGVSError
    
    try:
        hv.validate(hp.parse_hgvs_variant("NM_001166478.1:c.30_32insT"))
    except HGVSError as e:
        print(e)


.. parsed-literal::

    insertion length must be 1


