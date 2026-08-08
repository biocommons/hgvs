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
argument to :class:`hgvs.parsers.Parser`, which dispatches (via ``__new__``) to
a backend-specific subclass:

* **OMeta/parsley** (``grammar_fn="__ometa__"``), the default, is implemented
  by ``hgvs.parsers.parsley_parser.ParsleyParser``. The grammar is written in
  ``src/hgvs/parsers/_data/hgvs.pymeta`` and compiled offline by
  ``sbin/generate_parser.py`` into
  ``src/hgvs/parsers/generated/hgvs_grammar.py``.
* **pyparsing** (``grammar_fn="__pyparsing__"``) is a newer implementation,
  ``hgvs.parsers.pyparsing_parser.PyParsingParser``, backed by the grammar in
  ``src/hgvs/parsers/pyparsing_grammar.py``, where the grammar is ordinary
  Python rather than a separate file plus a code-generation step. It is
  expected to become the default in a future major release, after which
  OMeta/parsley will be deprecated and removed. See `issue #505
  <https://github.com/biocommons/hgvs/issues/505>`_.

Each backend's dependencies (``parsley``/``ometa`` or ``pyparsing``) are
imported only when that backend is selected.

The two are intended to be interchangeable and are tested for equivalence. The
one deliberate difference is that the OMeta grammar silently skips leading
spaces and tabs, while the pyparsing grammar rejects them; HGVS strings contain
no whitespace, so the stricter reading is intended.

.. include:: hgvs_railroad.rst


.. toctree::
   :hidden:

   hgvs_railroad
