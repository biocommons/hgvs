Grammar
@@@@@@@

Grammar Overview
################

.. note:: This section is being written.

  Provide an overview of the grammar rules
  Also consider a document link to the grammar itself

Grammar Implementations
#######################

Two grammar implementations are shipped, selected with the ``grammar_fn``
argument to :class:`hgvs.parser.Parser`:

* **OMeta/parsley** (``grammar_fn="__ometa__"``) is the default. The grammar is
  written in ``src/hgvs/_data/hgvs.pymeta`` and compiled offline by
  ``sbin/generate_parser.py`` into ``src/hgvs/generated/hgvs_grammar.py``.
* **pyparsing** (``grammar_fn="__pyparsing__"``) is a newer implementation in
  ``src/hgvs/grammar.py``, where the grammar is ordinary Python rather than a
  separate file plus a code-generation step. It is expected to become the
  default in a future major release, after which OMeta/parsley will be
  deprecated and removed. See `issue #505
  <https://github.com/biocommons/hgvs/issues/505>`_.

The two are intended to be interchangeable and are tested for equivalence. The
one deliberate difference is that the OMeta grammar silently skips leading
spaces and tabs, while the pyparsing grammar rejects them; HGVS strings contain
no whitespace, so the stricter reading is intended.

.. include:: hgvs_railroad.rst


.. toctree::
   :hidden:

   hgvs_railroad
