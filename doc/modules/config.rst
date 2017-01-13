Configuration
.............

:mod:`hgvs.config`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: hgvs.config


The defaults are::

  [mapping]
  alt_aln_method = splign
  assembly = GRCh38
  in_par_assume = X
  inferred_p_is_uncertain = True
  normalize = True
  replace_reference = True
  
  [formatting]
  p_3_letter = True
  p_term_asterisk = False
  
  [normalizer]
  cross_boundaries = False
  shuffle_direction = 3
  validate = True
  window_size = 20
  
  [lru_cache]
  maxsize = 100

