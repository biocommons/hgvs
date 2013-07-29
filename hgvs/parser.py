import parsley
import os


class Parser(object):
    grammer_fn = os.path.join(os.path.dirname(__file__),'grammar.txt')
    parsley.makeGrammar(open('hgvs/hgvs.parsley').read(),{})
