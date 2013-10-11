import recordtype

class Variant( recordtype.recordtype('Variant', ['seqref','type','posedit']) ):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an hgvs.location.CDSInterval (for example) are both intended uses
    """
    
    def __str__(self):
        return '{self.seqref}:{self.type}.{self.posedit}'.format(self=self)
