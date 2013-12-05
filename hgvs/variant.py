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



if __name__ == '__main__':
# expect output like this:
# (default2.7)snafu$ PYTHONPATH=. python ./hgvs/variant.py
# a
# ./hgvs/variant.py:27: DeprecationWarning: seqref is deprecated; use ac instead
#   print(v.seqref)
# a
# a
# a
# (no error on second call to t() because warning has already been reported)

    import hgvs.variant
    def t():
        v = hgvs.variant.SequenceVariant('a','b','c')
        print(v.ac)
        print(v.seqref)

    warnings.simplefilter("default")
    t()
    t()
    
