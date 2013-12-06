"""
MultiFastaDB -- present a collection of indexed fasta files as a single source

from multifastafile import MultiFastaDB
mff = MultiFastaDB([ '/a/file.fasta', '/a/dir/of/fastas/' ])
mff.fetch('NM_01234.5',60,70)
mff['NM_01234.5'][60:70] (equivalent to the above)
mff.where_is('NM_01234.5')

NOTE: This is temporarily in the HGVS package. It will likely move to a
separate package of bioinformatics tools.
"""


import collections
import itertools
import logging
import os

import pysam

class MultiFastaDB(object):
    """Class that opens a set of indexed fasta files and provides ordered
    lookup of a given accession across all of them.  The intent is to
    simplify accessing a virtual database of sequences that is distributed
    across multiple files."""

    class SequenceProxy(object):

        """Represents a sequence fetch but defers the actual fetch until a
        slice is applied, allowing convenient and sexy slice syntax like this:
        frsf = MultiFastaDB('/my/seqs')
        seq = frsf['NM_01234.5'][4:20]
        """

        def __init__(self,frsf,ac):
            self.frsf = frsf
            self.ac = ac

        def __getslice__(self,start_i,end_i):
            return self.frsf.fetch(self.ac,start_i,end_i)

        def __getitem__(self,i):
            return self[i,i+1]

        def __str__(self):
            return self.frsf.fetch(self.ac)

    def __init__(self,sources=[],suffixes=['fa','fasta','fa.gz','fasta.gz']):
        self._fastas = None
        self._logger = logging
        self.sources = sources
        self.suffixes = suffixes
        self.open_sources()

    def open_sources(self):
        """Opens or reopens fasta sources (directories or files) provided
        when the instance was created."""
        self._fastas = collections.OrderedDict()
        for source in self.sources:
            fa_paths = []
            if os.path.isfile(source):
                fa_paths = [ source ]
            elif os.path.isdir(source):
                # sort so that searches have defined ordering
                fa_paths = sorted([ os.path.join(source,de)
                                    for de in os.listdir(source)
                                    if any(de.endswith(sfx) for sfx in self.suffixes) ])
            else:
                raise RuntimeError(source + ": expected a directory or file for fasta sources")
            for fa_path in fa_paths:
                self._logger.debug("opening "+fa_path)
                fai_path = fa_path+'.fai'
                if (os.path.exists(fai_path)
                    and os.stat(fa_path).st_mtime > os.stat(fai_path).st_mtime):
                    self._logger.warn(fai_path + " is out-of-date (older than fasta file)")
                self._fastas[fa_path] = pysam.Fastafile(fa_path)

    def where_is(self,ac):
        """return list of all (filename,pysam.Fastafile) pairs in which
        accession occurs"""
        return [ (fp,fh)
                 for fp,fh in self._fastas.iteritems()
                 if ac in fh ]

    def fetch(self,ac,start_i=None,end_i=None):
        """return a sequence, or subsequence if start_i and end_i are provided"""
        for fah in self._fastas.values():
            seq = fah.fetch(ac,start_i,end_i)
            if seq != '':
                return seq
        return None

    @property
    def references(self):
        return list( itertools.chain.from_iterable([
            fa.references for fa in self._fastas.values() ]))
        
    @property
    def lengths(self):
        return list( itertools.chain.from_iterable([
            fa.lengths for fa in self._fastas.values() ]))


    def __contains__(self,ac):
        return any([ ac in fh for fh in self._fastas.itervalues() ])

    def __getitem__(self,ac):
        return self.SequenceProxy(self,ac)
