Introduction
!!!!!!!!!!!!

Genome, transcript, and protein sequence variants are typically
reported using the `variation nomenclature ("varnomen")
recommendations <http://varnomen.hgvs.org/>`_ provided by the `Human
Genome Variation Society (HGVS) <http://www.hgvs.org/>`_ (`Taschner
and den Dunnen, 2011 <http://www.ncbi.nlm.nih.gov/pubmed/21309030>`_).
Most variants are deceptively simple looking, such as
NM_021960.4:c.740C>T. In reality, the varnomen standard provides for
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
breadth of the varnomen standard makes it difficult to implement the
standard in software, which in turn makes using the standard in
high-throughput analyses difficult.

This package, `hgvs`, is an easy-to-use Python library for parsing,
representing, formatting, and mapping variants between genome, transcript,
and protein sequences.  The current implementation handles most (but not
all) of the varnomen standard for precisely defined sequence variants.
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
  `UTA`_ by default, but the package includes a well-defined service
  interface that enables alternative data sources.
* **Extensive automated tests**. We run extensive automated tests
  consisting of all supported variant types on many genes for every
  single commit to the source code repository. Test results are
  displayed publicly and immediately.


.. note:: **Some HGVS recommendations are intentionally absent.** This
   package is primarily concerned with the subset of the `VarNomen`_
   recommendations that are relevant for high-throughput
   sequencing. See `issues`_ for a full set of bugs and feature
   requests.




Related tools
@@@@@@@@@@@@@

* `Mutalyzer <http://www.humgen.nl/mutalyzer.html>`_ provides a web
  interface to variant validation and mapping.
* `Counsyl hgvs package <https://github.com/counsyl/hgvs>`_ provides
  functionality conceptually similar to that of the Invitae hgvs
  package.


Support
@@@@@@@

See the section :doc:`getting_help` for information about connecting
with the community, asking questions, and :ref:`filing bug reports
correctly<bug-reports>`.


Links
@@@@@

* `Variation Nomenclature Recommendations <http://varnomen.hgvs.org/>`_
* `Human Genome Variation Society (HGVS) <http://www.hgvs.org/>`_
* Parsley_, an Python wrapper for the OMeta Parser Expression Grammar (PEG_)
* `Universal Transcript Archive (UTA) <https://github.com/biocommons/uta/>`_


References
@@@@@@@@@@

hgvs: A Python package for manipulating sequence variants using HGVS nomenclature: 2018 Update.
  | Wang M, Callenberg KM, Dalgleish R, Fedtsov A, Fox N, Freeman PJ, Jacobs KB, Kaleta P, McMurry AJ, Prlić A, Rajaraman V, Hart RK
  | Human Mutation. 2018
  | https://www.ncbi.nlm.nih.gov/pubmed/30129167

A Python package for parsing, validating, mapping and formatting sequence variants using HGVS nomenclature.
  | Hart RK, Rico R, Hare E, Garcia J, Westbrook J, Fusaro VA
  | Bioinformatics. 31(2):268-70 (2014).
  | https://www.ncbi.nlm.nih.gov/pubmed/25273102

Describing structural changes by extending HGVS sequence variation nomenclature.
  | Taschner, P. E. M., & den Dunnen, J. T.
  | Human Mutation, 32(5), 507–11. (2011).
  | http://www.ncbi.nlm.nih.gov/pubmed/21309030

A formalized description of the standard human variant nomenclature in Extended Backus-Naur Form.
  | Laros, J. F. J., Blavier, A., den Dunnen, J. T., & Taschner, P. E. M.
  | BMC Bioinformatics, 12 Suppl 4(Suppl 4), S5. (2011). 
  | http://www.ncbi.nlm.nih.gov/pubmed/21992071

