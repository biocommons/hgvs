"""
populate __version__ and hg_id

Two cases:
1) Try to fetch SCM info by running 'hg id'.  This case applies
   mostly to development code.
2) If that fails, look for _release.py.
   _release.py is created during packaging, but is NOT added in hg.
   Therefore, this case applies mostly to deployed code that
   doesn't have a .hg directory somewhere along the path to __file__.
"""

from hgvs.utils.hgid import HgId

try:
    
    hg_id = HgId.init_from_hg(__file__)
    __version__ = hg_id.version

except:

    try:

        from _release import hg_id
        hg_id = HgId.init_from_string( hg_id )
        __version__ = hg_id.version

    except ImportError:

        hg_id = None
        __version__ = None

