import os, re, subprocess

class HgId(object):
    __slots__ = ['root', 'branch', 'id', 'tags']

    def __init__(self,root,id,branch,tags):
        self.root = root
        self.branch = branch
        self.id = id
        self.tags = tags if isinstance(tags,list) else tags.split()
        
    @staticmethod
    def init_from_string(root,id_string):
        """
        create HgId instance from output of `hg id -bit`
        The string should look like:
        95aa60212dd4 default foo goo
        where foo and goo are tags on this changeset
        """
        args = id_string.split(' ',2)
        return HgId(root,*args)

    @staticmethod
    def init_from_hg(path):
        """
        create HgId instance by calling `hg -R <root> id -bit`, where
        root is determined as the nearest ancestor directory with a .hg directory
        """
        def find_hg_root(path):
            "return hg root for a given file/dir path (e.g., package.__file__)"
            root = os.path.abspath(path)
            while True:
                if os.path.exists( os.path.join( root, '.hg' ) ):
                    return root
                next_root = os.path.dirname(root)
                if root == next_root:
                    break
                root = next_root
            return None
        
        root = find_hg_root(path)
        if root is None:
            raise RuntimeError('.hg directory not found')
        hg_id_out = subprocess.check_output([
            'hg', '-R', root, 'id', '-bit', 
            '--config', 'trusted.users=*'])    # prevents Not trusting file... message
        return HgId.init_from_string(root,hg_id_out.rstrip())

    @property
    def tag0(self):
        return self.tags[0]

    @property
    def semver(self):
        """Return Semantic Version (http://semver.org/) if the tag
        has that format, else None"""
        matches = [ re.match('\d+\.\d+\.\d+.*',t) for t in self.tags ]
        return matches[0].group(0) if any(matches) else None

    @property
    def version(self):
        """Return semver if not None, otherwise the first 7 chars of the
        changeset it (same as used by BitBucket"""
        return self.semver or self.id[:7]

    def __str__(self):
        return '{self.tag0} ({self.id} {self.branch} {tags_str})'.format(
            self = self, tags_str = ' '.join(self.tags))
            

if __name__ == '__main__':
    import pprint
    hg_id = HgId.init_from_hg(__file__)
    print( 'repr:', repr(hg_id) )
    print( 'str:', str(hg_id) )
