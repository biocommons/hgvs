.SUFFIXES:
.PHONY: FORCE

RDG:=${HOME}/opt/railroad-diagram-generator/build/railroad-diagram-generator

default: hgvs-syntax-diagram.html

.PRECIOUS: %.ebnf
%.ebnf: ../../hgvs/data/%.pymeta FORCE
	./pymeta-to-ebnf $< >$@

.PRECIOUS: %.ebnf_grammar/all.html
%.ebnf_grammar/all.html: %.ebnf
	${RDG} $<

%-syntax-diagram.html: hgvs.ebnf_grammar/all.html
	./redact-html $< >$@


.PHONY: clean cleaner cleanest
clean:
	/bin/rm -f *~
cleaner: clean
cleanest: cleaner
	/bin/rm -fr hgvs.ebnf hgvs.ebnf_grammar hgvs-syntax-diagram.html
