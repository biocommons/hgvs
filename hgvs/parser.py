import os

import parsley

import hgvs.location

grammar_fn = os.path.join(os.path.dirname(__file__),'grammar.txt')
grammar = parsley.makeGrammar(open(grammar_fn,'r').read(),{
    'Position': hgvs.location.Position,
    'CDSPosition': hgvs.location.CDSPosition,
    'Interval': hgvs.location.Interval,
    'CDSInterval': hgvs.location.CDSInterval,
    })

