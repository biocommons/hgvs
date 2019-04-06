
Creating a SequenceVariant from scratch
---------------------------------------

0. Overview
~~~~~~~~~~~

A SequenceVariant consists of an accession (a string), a sequence type
(a string), and a PosEdit, like this:

var = hgvs.sequencevariant.SequenceVariant(ac=‘NM_01234.5’, type=‘c’,
posedit=…)

Unsurprisingly, a PosEdit consists of separate position and Edit
objects. A position is generally an Interval, which in turn is comprised
of SimplePosition or BaseOffsetPosition objects. An edit is a subclass
of Edit, which includes classes like NARefAlt for substitutions,
deletions, and insertions) and Dup (for duplications).

Importantly, each of the objects we’re building has a rule in the
parser, which means that you have the tools to serialize and deserialize
(parse) each of the components that we’re about to construct.

1. Make an Interval to define a position of the edit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import hgvs.location
    import hgvs.posedit

.. code:: ipython3

    start = hgvs.location.BaseOffsetPosition(base=200,offset=-6,datum=hgvs.location.Datum.CDS_START)
    start, str(start)




.. parsed-literal::

    (BaseOffsetPosition(base=200, offset=-6, datum=Datum.CDS_START, uncertain=False),
     '200-6')



.. code:: ipython3

    end = hgvs.location.BaseOffsetPosition(base=22,datum=hgvs.location.Datum.CDS_END)
    end, str(end)




.. parsed-literal::

    (BaseOffsetPosition(base=22, offset=0, datum=Datum.CDS_END, uncertain=False),
     '*22')



.. code:: ipython3

    iv = hgvs.location.Interval(start=start,end=end)
    iv, str(iv)




.. parsed-literal::

    (Interval(start=200-6, end=*22, uncertain=False), '200-6_*22')



2. Make an edit object
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import hgvs.edit, hgvs.posedit

.. code:: ipython3

    edit = hgvs.edit.NARefAlt(ref='A',alt='T')
    edit, str(edit)




.. parsed-literal::

    (NARefAlt(ref='A', alt='T', uncertain=False), 'A>T')



.. code:: ipython3

    posedit = hgvs.posedit.PosEdit(pos=iv,edit=edit)
    posedit, str(posedit)




.. parsed-literal::

    (PosEdit(pos=200-6_*22, edit=A>T, uncertain=False), '200-6_*22A>T')



3. Make the variant
~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import hgvs.sequencevariant

.. code:: ipython3

    var = hgvs.sequencevariant.SequenceVariant(ac='NM_01234.5', type='c', posedit=posedit)
    var, str(var)




.. parsed-literal::

    (SequenceVariant(ac=NM_01234.5, type=c, posedit=200-6_*22A>T),
     'NM_01234.5:c.200-6_*22A>T')



**Important: It is possible to bogus variants with the hgvs package. For
example, the above interval is incompatible with a SNV. See
hgvs.validator.Validator for validation options.**

4. Update your variant
~~~~~~~~~~~~~~~~~~~~~~

The stringification happens on-the-fly. That means that you can update
components of the variant and see the effects immediately.

.. code:: ipython3

    import copy

.. code:: ipython3

    var2 = copy.deepcopy(var)
    var2.posedit.pos.start.base=456
    str(var2)




.. parsed-literal::

    'NM_01234.5:c.456-6_*22A>T'



.. code:: ipython3

    var2 = copy.deepcopy(var)
    var2.posedit.edit.alt='CT'
    str(var2)




.. parsed-literal::

    'NM_01234.5:c.200-6_*22delinsCT'



.. code:: ipython3

    var2 = copy.deepcopy(var)
    str(var2)




.. parsed-literal::

    'NM_01234.5:c.200-6_*22A>T'



