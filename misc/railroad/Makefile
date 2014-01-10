.SUFFIXES:
.PHONY: FORCE

RDG:=${HOME}/opt/railroad-diagram-generator/build/railroad-diagram-generator

default: hgvs-syntax-diagram.html

%.ebnf: ../../hgvs/data/%.pymeta FORCE
	./pymeta-to-ebnf $< >$@

%.ebnf_grammar/all.html: %.ebnf
	${RDG} $<

%-syntax-diagram.html: hgvs.ebnf_grammar/all.html
	cp $< $@
