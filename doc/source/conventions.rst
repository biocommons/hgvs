Conventions
-----------

.. note:: This section is being written.

Variables
~~~~~~~~~

The following code variable conventions are used for most of the ``hgvs``
code base.  They should be considered aspirations rather than reality or
policy.  Understanding these conventions will help uses and developers
understand the code.

.. note:: A note on variable suffixes
  If a particular variant type is expected, a suffix is often added to
  variable names. |eg| ``var_c`` in a function argument list signifies
  that a SequenceVariant object with type='c' is expected.

:hgvs*: a string representing an HGVS variant name.  

:var*: a :class:`hgvs.variant.SequenceVariant` object

:pos: 

:posedit: 

:hgvs_position:


.. |eg| replace:: *e.g.,*
