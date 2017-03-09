.SUFFIXES:
.PHONY: FORCE

SHELL:=/bin/bash

RDG:=${HOME}/opt/railroad-diagram-generator/build/railroad-diagram-generator

default: hgvs-syntax-diagram.html hgvs_railroad.rst

.PRECIOUS: %.pymeta
%.pymeta:: ../../hgvs/_data/%.pymeta
	ln -s $<

.PRECIOUS: %.ebnf
%.ebnf: %.pymeta
	./pymeta-to-ebnf $< >$@

.PRECIOUS: %.ebnf_grammar/all.html
%.ebnf_grammar/all.html: %.ebnf
	${RDG} $<

%-syntax-diagram.html: hgvs.ebnf_grammar/all.html
	./redact-html $< >$@

hgvs_railroad.rst: %.rst: hgvs.ebnf_grammar/all.html
	/bin/mkdir -p $*
	/bin/cp -l ${<D}/*.svg $*/
	./make-hgvs-railroad-doc --hgvs-version "$$(hg id -bit)" hgvs_railroad/*svg >$@

doc-install:
	rm -fr ../../doc/source/hgvs_railroad*
	make hgvs_railroad.rst && /bin/mv -b hgvs_railroad* ../../doc/source/


.PHONY: clean cleaner cleanest
clean:
	/bin/rm -f *~
cleaner: clean
cleanest: cleaner
	/bin/rm -fr hgvs.ebnf hgvs.ebnf_grammar hgvs-syntax-diagram.html hgvs_railroad
