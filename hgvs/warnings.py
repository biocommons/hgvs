# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

"""HGVS Warnings

The `hgvs` package issues warnings under specific circumstances during
mapping and validation. 

Types of Warnings
=================

Warnings have three basic classes: 
- transcript warnings
- coordinate mapping warnings
- sequence warnings

These classes of warnings, and specific causes, are represented as
subclasses of Python's Warning class. Details are given below.


Issuing and Capturing Warnings
==============================

An HGVSWarning (or subclass) should be issued like this:

  import warnings
  warnings.warn( HGVSWarning( message, optional_data_dict ) )

HGVSWarning may be any of the HGVSWarning subclasses.

Warnings (in general) may be captured using a context manager:

  with warnings.catch_warnings(record=True) as warning_messages:
    warnings.simplefilter('always')
    do_something()
    print("{} messages found".format(len(warning_messages)))

"""

import contextlib
import warnings

#* cross-source exon structure ambiguity (e.g., ncbi v. ucsc)
#* intra-source exon structure ambiguity (e.g., ucsc paralogous alignments)
#* alignment indel count > x events (or x events/kb)
#* alignment within indel
#* alignment adjacent to indel


# TODO: implement hgvs-specific context manager
#class HGVSWarningContextManager(object):
#    def __init__(self):
#        pass
#    def __enter__(self):
#        with warnings.catch_warnings(record=True) as self.warnings_messages:
#            warnings.simplefilter('always')
#    def __exit__(self, type, value, traceback):
        
# TODO: implement warning severity map


class HGVSWarning(Warning):
    """Base class for warnings emitted by the `hgvs` package.
    """

    def __init__(self,message,data={}):
        super(HGVSWarning,self).__init__(message,data)
        self.type = self.__class__.__name__
        self.message = message
        self.data = data

    def __str__(self):
        return '{self.type}: {self.message} ({self.data})'.format(self=self)



class TranscriptWarning(HGVSWarning):
    pass

class AmbiguousExonStructure(TranscriptWarning):
    pass

class ComplexAlignment(TranscriptWarning):
    pass

class ReferenceSNVDiscrepancy(TranscriptWarning):
    pass

class ReferenceIndelDiscrepancy(TranscriptWarning):
    pass



class MappingWarning(HGVSWarning):
    pass

class MappingOverlapsSNV(MappingWarning):
    pass

class MappingOverlapsIndel(MappingWarning):
    pass



class SequenceWarning(HGVSWarning):
    pass
