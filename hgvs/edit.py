import recordtype

from hgvs.exceptions import HGVSError

class Edit(object):
    pass

class DelIns( Edit, recordtype.recordtype('DelIns', ['pre','post'], default=None) ):
    def __str__(self):
        if self.pre is None and self.post is None:
            raise HGVSError('DelIns: pre and post sequences are both empty')
        if self.pre is not None and self.post is not None:
            if self.pre == self.post:
                return '='
            if len(self.pre) == 1 and len(self.post) == 1:
                return '{self.pre}>{self.post}'.format(self=self)
            return 'del{self.pre}ins{self.post}'.format(self=self)
        if self.pre is not None and self.post is None:
            return 'del{self.pre}'.format(self=self)
        if self.pre is None and self.post is not None:
            return 'ins{self.post}'.format(self=self)

class Dup( Edit, recordtype.recordtype('Dup', ['seq'], default=None) ):
    def __str__(self):
        return 'dup' + (self.seq or '')

class Repeat( Edit, recordtype.recordtype('Repeat', ['seq','min','max'], default=None) ):
    def __str__(self):
        if self.min > self.max:
            raise HGVSError('Repeat min count must be less than or equal to max count')
        return '{self.seq}({self.min}_{self.max})'.format(self=self)
