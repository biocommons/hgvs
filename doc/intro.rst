Introduction
!!!!!!!!!!!!

Genome, transcript, and protein sequence variants are typically
reported using the `mutation nomenclature ("mutnomen") recommendations
<http://www.hgvs.org/mutnomen/>`_ provided by the `Human Genome
Variation Society (HGVS) <http://www.hgvs.org/>`_ (`Taschner and den
Dunnen, 2011 <http://www.ncbi.nlm.nih.gov/pubmed/21309030>`_).  Most
variants are deceptively simple looking, such as
NM_021960.4:c.740C>T. In reality, the mutnomen standard provides for
much more complex concepts and representations.

As high-throughput sequencing becomes commonplace in the investigation
and diagnosis of disease, it is essential that communicating variants
from sequencing projects to the scientific community and from
diagnostic laboratories to health care providers is easy and
accurate. The HGVS mutation nomenclature recommendations⁠ are generally
accepted for the communication of sequence variation: they are widely
endorsed by professional organizations, mandated by numerous journals,
and the prevalent representation used by databases and interactive
scientific software tools. The guidelines – originally devised to
standardize the representation of variants discovered before the
advent of high-throughput sequencing – are now approved by the HGVS
and continue to evolve under the auspices of the Human Variome
Project. Unfortunately, the complexity of biological phenomena and the
breadth of the mutnomen standard makes it difficult to implement the
standard in software, which in turn makes using the standard in
high-throughput analyses difficult.

This package, ``hgvs``, is an easy-to-use Python library for parsing,
representing, formatting, and mapping variants between genome, transcript,
and protein sequences.  The current implementation handles most (but not
all) of the mutnomen standard for precisely defined sequence variants.
The intent is to centralize the subset of HGVS variant manipulation that
is routinely used in modern, high-throughput sequencing analysis.


.. _features:

Features of the hgvs Package
@@@@@@@@@@@@@@@@@@@@@@@@@@@@

* **Convenient object representation**. Manipulate variants
  conceptually rather than by modifying text strings. Classes model
  HGVS concepts such as :class:`Interval <hgvs.location.Interval>`,
  intronic offsets (in :class:`BaseOffsetPosition
  <hgvs.location.BaseOffsetPosition>`), uncertainty, and types of
  variation (:mod:`hgvs.edit`).
* **A grammar-based parser**. `hgvs` uses :doc:`a formal grammar
  <hgvs_railroad>` to parse HGVS variants rather than string
  partitioning or regular expression pattern matching.  This makes
  parsing easier to understand, extend, and validate.
* **Simple variant formatting**. Object representations of variants
  may be turned into HGVS strings simply by printing or "stringifying"
  them.
* **Robust variant mapping**. The package includes tools to map variants between
  genome, transcript, and protein sequences (:class:`VariantMapper
  <hgvs.variantmapper.VariantMapper>` and to perform liftover between
  two transcript via a common reference (:class:`Projector
  <hgvs.projector.Projector>`).  The hgvs mapper is specifically
  designed to reliably handl of regions reference-transcript indel
  discrepancy that are not covered by other tools.
* **Additional variant validation**. The package includes tools to
  validate variants, separate from syntactic validation provided by
  the grammar.
* **Extensible data sources**. Mapping and sequence data come from
  `UTA <https://bitbucket.org/biocommons/uta/>`_ by default, but the
  package includes a well-defined service interface that enables
  alternative data sources.
* **Extensive automated tests**. We run extensive automated tests
  consisting of all supported variant types on many genes for every
  single commit to the source code repository. Test results are
  displayed publicly and immediately.


.. _limitations:

Current limitations of the hgvs Package
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

.. note::

   All issues are public. For a full set of bugs, feature requests,
   and tasks, see `hgvs issues
   <https://bitbucket.org/biocommons/hgvs/issues?status=new&status=open>`__.

* **Compound, complex, and mosaic variants will be available in the next release (0.5)**.
  These have been implemented, and are awaiting further testing before
  merging.  See :issue:`104`.

* **Some HGVS recommendations are intentionally absent.**. We chose to focus on the subset
  of the HGVS recommendations that are relevant for high-throughput
  sequencing. Some features, such as translocations, are not currently
  on our roadmap.





Related tools
@@@@@@@@@@@@@

* `Mutalyzer <http://www.humgen.nl/mutalyzer.html>`_ provides a web
  interface to variant validation and mapping.
* `Counsyl hgvs package <https://github.com/counsyl/hgvs>`_ provides
  functionality conceptually similar to that of the Invitae hgvs
  package.


Getting Help
@@@@@@@@@@@@

There are several ways to get help with the ``hgvs`` package.

The `hgvs-discuss mailing list
<https://groups.google.com/forum/#!forum/hgvs-discuss>`_ is the preferred
way to reach the ``hgvs`` package authors.  Please file bugs and feature
requests on the `hgvs issue tracker
<https://bitbucket.org/biocommons/hgvs/issues?status=new&status=open>`_.

If you have questions about the `HGVS Mutation Nomenclature Recommendations
<http://www.hgvs.org/mutnomen/>`_, consider posting your questions to the
`HGVS Facebook page <https://www.facebook.com/HGVSmutnomen>`_.


Links
@@@@@

* `HGVS Mutation Nomenclature Recommendations <http://www.hgvs.org/mutnomen/>`_
* `Human Genome Variation Society (HGVS) <http://www.hgvs.org/>`_
* `Parsley <https://pypi.python.org/pypi/Parsley>`_
* `Universal Transcript Archive (UTA) <https://bitbucket.org/biocommons/uta/>`_


References
@@@@@@@@@@

Describing structural changes by extending HGVS sequence variation nomenclature.
  | Taschner, P. E. M., & den Dunnen, J. T.
  | Human mutation, 32(5), 507–11. (2011).
  | http://www.ncbi.nlm.nih.gov/pubmed/21309030

A formalized description of the standard human variant nomenclature in Extended Backus-Naur Form.
  | Laros, J. F. J., Blavier, A., den Dunnen, J. T., & Taschner, P. E. M.
  | BMC bioinformatics, 12 Suppl 4(Suppl 4), S5. (2011). 
  | http://www.ncbi.nlm.nih.gov/pubmed/21992071





.. _`Parsley`: https://pypi.python.org/pypi/Parsley
.. _`HGVS`: http://www.hgvs.org/
.. _`HGVS Recommendations`: http://hgvs.org/mutnomen/
.. _PEG: http://en.wikipedia.org/wiki/Parsing_expression_grammar
