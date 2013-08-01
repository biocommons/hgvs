import recordtype

class Edit(object):
    pass

class DelIns( Edit, recordtype.recordtype('DelIns', ['pre','post'], default=None) ):
    def __str__(self):
        assert self.pre or self.post, 'DelIns: pre and post sequences are both empty'
        if self.pre and self.post:
            if self.pre == self.post:
                return '='
            if len(self.pre) == 1 and len(self.post) == 1:
                return '{self.pre}>{self.post}'.format(self=self)
            return 'del{self.pre}ins{self.post}'.format(self=self)
        if self.pre and not self.post:
            return 'del{self.pre}'.format(self=self)
        if not self.pre and self.post:
            return 'ins{self.post}'.format(self=self)

class Dup( Edit, recordtype.recordtype('Dup', ['seq'], default=None) ):
    def __str__(self):
        return 'dup' + (self.seq or '')

class Repeat( Edit, recordtype.recordtype('Repeat', ['seq','min','max'], default=None) ):
    def __str__(self):
        assert self.min <= self.max, 'Repeat min count must be less than or equal to max count'
        return '{self.seq}({self.min}_{self.max})'.format(self=self)
