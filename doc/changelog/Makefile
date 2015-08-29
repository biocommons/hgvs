# Note: You must have biocommons.pytemplate/bin in your PATH

SHELL:=/bin/bash

VTAGS:=$(shell hg tags | cut -f1 -d' ' | grep -E '^[0-9]+\.' | fgrep -v 0.0.0)

default all: $(addsuffix .rst,${VTAGS})
	make {0.0,0.1,0.2,0.3,0.4}.rst

.PRECIOUS: %.clog
%.clog:
	clogger $* | perl -pe 's/CORE-(\d+)/#\1/' >$@ >$@.tmp
	mv $@.tmp $@

%.rst: %.clog
	clogger-fmt \
		-I '`#{issue_id} <https://bitbucket.org/biocommons/hgvs/issues/{issue_id}/>`_' \
		-C '`{cset} <https://bitbucket.org/biocommons/hgvs/commits/{cset}>`_' \
		<$< >$@.tmp
	mv $@.tmp $@

0.0.rst 0.1.rst 0.2.rst 0.3.rst 0.4.rst 0.5.rst: %.rst: $(wildcard %.[0-9]*.rst)
	clogger-group-XY $* >$@.tmp
	mv $@.tmp $@


.PHONY: clean cleaner cleanest

clean:
	/bin/rm -f *~

cleaner: clean
	/bin/rm -f *.rst

cleanest: cleaner
	/bin/rm -f *.clog
