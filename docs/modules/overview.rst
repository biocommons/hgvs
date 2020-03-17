Module Overview
...............

+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| Module                                  |Classes                                              |Description                              |
+=========================================+=====================================================+=========================================+
|                                                                                                                                         |
| *Variant Object Representation*                                                                                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.edit`                        | | :class:`hgvs.edit.AAExt`                          |:mod:`hgvs.edit` classes implement       |
|                                         | | :class:`hgvs.edit.AAFs`                           |various kinds of sequence edits. For     |
|                                         | | :class:`hgvs.edit.AARefAlt`                       |nucleic acids, these edits are           |
|                                         | | :class:`hgvs.edit.AASub`                          |independent of location; amino acids     |
|                                         | | :class:`hgvs.edit.Dup`                            |edits currently contain the location.    |
|                                         | | :class:`hgvs.edit.Edit`                           |                                         |
|                                         | | :class:`hgvs.edit.NACopy`                         |                                         |
|                                         | | :class:`hgvs.edit.NADupN`                         |                                         |
|                                         | | :class:`hgvs.edit.NARefAlt`                       |                                         |
|                                         | | :class:`hgvs.edit.Repeat`                         |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.hgvsposition`                | | :class:`hgvs.hgvsposition.HGVSPosition`           |A non-standard representation of a       |
|                                         |                                                     |sequence location without an edit. For   |
|                                         |                                                     |example, NM_012345.6:c.72+5_73-2.        |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.location`                    | | :class:`hgvs.location.AAPosition`                 |Various kinds of locations. Interval is a|
|                                         | | :class:`hgvs.location.BaseOffsetPosition`         |span from ``start`` to ``end``; the      |
|                                         | | :class:`hgvs.location.Interval`                   |others are points in a sequence.         |
|                                         | | :class:`hgvs.location.SimplePosition`             |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.posedit`                     | | :class:`hgvs.posedit.PosEdit`                     |A position+edit (really, an interval and |
|                                         |                                                     |edit).                                   |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.variant`                     | | :class:`hgvs.variant.SequenceVariant`             |A sequence variant of any type (g, c, m, |
|                                         |                                                     |r, n, p). A SequenceVariant is returned  |
|                                         |                                                     |by :class:`hgvs.parser.Parser`, and it is|
|                                         |                                                     |the input and output type for            |
|                                         |                                                     |:class:`hgvs.variantmapper.VariantMapper`|
|                                         |                                                     |operations.                              |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
|                                                                                                                                         |
| *Parsing and Formatting*                                                                                                                |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.parser`                      | | :class:`hgvs.parser.Parser`                       |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
|                                                                                                                                         |
| *Coordinate, Interval, and Variant Mapping/Transformation*                                                                              |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.projector`                   | | :class:`hgvs.projector.Projector`                 |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.alignmentmapper`             | | :class:`hgvs.alignmentmapper.AlignmentMapper`     |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.variantmapper`               | | :class:`hgvs.variantmapper.VariantMapper`         |                                         |
|                                         | | :class:`hgvs.assemblymapper.AssemblyMapper`       |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
|                                                                                                                                         |
| *Variant Normalization and Validation*                                                                                                  |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.normalizer`                  | | :class:`hgvs.normalizer.Normalizer`               |                                         |
|                                         |                                                     |                                         |
|                                         |                                                     |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.validator`                   | | :class:`hgvs.validator.Validator`                 |                                         |
|                                         | | :class:`hgvs.validator.IntrinsicValidator`        |                                         |
|                                         | | :class:`hgvs.validator.ExtrinsicValidator`        |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
|                                                                                                                                         |
| *External Data Providers*                                                                                                               |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.dataproviders.interface`     | | :class:`hgvs.dataproviders.interface.Interface`   |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
| :mod:`hgvs.dataproviders.uta`           | | :class:`hgvs.dataproviders.uta.UTABase`           |                                         |
+-----------------------------------------+-----------------------------------------------------+-----------------------------------------+
