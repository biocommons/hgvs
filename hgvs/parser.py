import os

import parsley

import hgvs.edit
import hgvs.hgvsposition
import hgvs.location
import hgvs.posedit
import hgvs.variant

from hgvs.exceptions import *

class Parser(object):
    default_grammar_fn = os.path.join(os.path.dirname(__file__),'grammar.txt')

    def __init__(self,grammar_fn=default_grammar_fn):
        self._grammar_fn = grammar_fn
        self._grammar = parsley.makeGrammar(open(grammar_fn,'r').read(),{
            'CDSPosition': hgvs.location.CDSPosition,
            'Interval': hgvs.location.Interval,
            'Position': hgvs.location.Position,

            'DelIns': hgvs.edit.DelIns,
            'Dup': hgvs.edit.Dup,
            'Repeat': hgvs.edit.Repeat,
            
            'PosEdit': hgvs.posedit.PosEdit,

            'HGVSPosition': hgvs.hgvsposition.HGVSPosition,
            'Variant': hgvs.variant.Variant,

            'hgvs': hgvs,
            })

        self.rules = [ m.replace('rule_','')
                       for m in dir(self._grammar._grammarClass)
                       if m.startswith('rule_') ]

    def __getattr__(self,name):
        # http://docs.python.org/2/reference/datamodel.html#object.__getattr__
        """
        attempt to call name as a grammar rule
        
        This is challenging because the call structure is self._grammar(s).name(),
        where s is the string to parse, which isn't provided to __getattr__.
        Notice that we need to call name() on an object that doesn't exist yet.

        This function returns a function that takes a string and returns the parsing result.
        """
        return lambda s: self._grammar(s).__getattr__(name)()
        

    def parse(self,variant):
        return self._grammar(variant).hgvs_variant()
        #try:
        #except Exception as e:
        #    raise HGVSParseError('{variant}: parsing raised {type} ({e.message})'.format(
        #        variant = variant, type = type(e), e = e))

        
