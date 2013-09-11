import recordtype

from hgvs.exceptions import HGVSError

class Edit(object):
    @property
    def allele1(self):
        return self.__str__()

    @property
    def allele2(self):
        return 'ref'

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

    # hacking these to return the required 'allele1' and 'allele2' fields
    @property
    def allele1(self):
        if self.pre is None and self.post is None:
            raise HGVSError('DelIns: pre and post sequences are both empty')
        if self.pre is not None and self.post is not None:
            if len(self.pre) == 1 and len(self.post) == 1:
                return self.pre
            return 'del{self.pre}ins{self.post}'.format(self=self)
        if self.pre is not None and self.post is None:
            return 'del{self.pre}'.format(self=self)
        if self.pre is None and self.post is not None:
            return 'ins{self.post}'.format(self=self)

    @property
    def allele2(self):
        if self.pre is None and self.post is None:
            raise HGVSError('DelIns: pre and post sequences are both empty')
        if self.pre is not None and self.post is not None:
            if len(self.pre) == 1 and len(self.post) == 1:
                return self.post
            return 'ref'
        if self.pre is not None and self.post is None:
            return 'ref'
        if self.pre is None and self.post is not None:
            return 'ref'

class Dup( Edit, recordtype.recordtype('Dup', ['seq'], default=None) ):
    def __str__(self):
        return 'dup' + (self.seq or '')


class Repeat( Edit, recordtype.recordtype('Repeat', ['seq','min','max'], default=None) ):
    def __str__(self):
        if self.min > self.max:
            raise HGVSError('Repeat min count must be less than or equal to max count')
        return '{self.seq}({self.min}_{self.max})'.format(self=self)

