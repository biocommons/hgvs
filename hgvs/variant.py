import warnings

import recordtype

class SequenceVariant( recordtype.recordtype('Variant', ['ac','type','posedit']) ):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """
    
    @property
    def seqref(self):
        warnings.warn('seqref is deprecated; use ac instead',DeprecationWarning,stacklevel=2)
        return self.ac

    def __str__(self):
        return '{self.ac}:{self.type}.{self.posedit}'.format(self=self)
