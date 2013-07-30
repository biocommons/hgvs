import parsley
import os


grammar_fn = os.path.join(os.path.dirname(__file__),'grammar.txt')
grammar = parsley.makeGrammar(open(grammar_fn,'r').read(),{})
