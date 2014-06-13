.. _index:

====
hgvs
====

.. include:: ../description.txt

|build_status| | `Source <https://bitbucket.org/invitae/hgvs>`_ | `Documentation <http://pythonhosted.org/hgvs/>`_ | `Discuss <https://groups.google.com/forum/#!forum/hgvs-discuss>`_ | `Issues <https://bitbucket.org/invitae/hgvs/issues?status=new&status=open>`_

Features
--------

* :doc:`A formal grammar <hgvs_railroad>` for HGVS variant names
* :doc:`Classes <modules>` that model HGVS concepts such as :class:`Interval
  <hgvs.location.Interval>`, intronic offsets (in
  :class:`BaseOffsetPosition <hgvs.location.BaseOffsetPosition>`), frameshifts, uncertain
  positions, and types of variation (:mod:`hgvs.edit`)
* Formatters that generate HGVS strings from internal representations
* Tools to map variants between genome, transcript, and protein sequences
  (:class:`VariantMapper <hgvs.variantmapper.VariantMapper>` and :class:`Projector
  <hgvs.projector.Projector>`)
* Reliable handling of regions reference-transcript discrepancy (requires
  `UTA <https://bitbucket.org/invitae/uta/>`_)
* Tools to validate variants (coming soon)
* Support for alternative sources of reference and transcript mapping
  information via a defined abstract interface
* Extensive automated tests


Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   getting_started

   installation
   examples

   reference

   getting_help
   license


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |build_status| image:: https://drone.io/bitbucket.org/invitae/hgvs/status.png
  :target: https://drone.io/bitbucket.org/invitae/hgvs
  :align: middle

