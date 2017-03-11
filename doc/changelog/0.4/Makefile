.PHONY: FORCE
.DELETE_ON_ERROR:

# N.B. this scripts requires tools that are not publicly available
# yet. Eventually, clogger will be released.
# Sorry.
PATH:=/home/reece/projects/reece/clogger/bin:${PATH}
SHELL:=/bin/bash -o pipefail

-include .tags.mk
.tags.mk: .tags
	perl -ane 'BEGIN {$$\=" "; print "VTAGS="}; print $$F[0] if m/^\d+\.\d+\./' $< >$@.tmp
	mv $@.tmp $@
.tags:
	hg tags >$@.tmp
	mv $@.tmp $@

XY:=0.4
XY_RST:=$(addsuffix .rst,${XY})
XY_VTAGS:=$(filter ${XY}%,${VTAGS})
XY_VTAGS_RST:=$(addsuffix .rst,${XY_VTAGS})

default:
	@echo "no $@ target"; exit 1

all: ${XY_VTAGS_RST} ${XY_RST}


${XY_RST}: %.rst:
	clogger-group-XY $* >$@.tmp
	mv $@.tmp $@


.PRECIOUS: %.clog
%.clog:
	clogger $* | perl -pe 's/CORE-(\d+)/#\1/' >$@.tmp
	mv $@.tmp $@
next.clog:
	clogger . | perl -pe 's/CORE-(\d+)/#\1/' >$@.tmp
	mv $@.tmp $@

# TODO: use git-mapfile to map commits
%.rst: %.clog
	clogger-fmt \
		-I '`#{issue_id} <https://github.com/biocommons/hgvs/issues/{issue_id}/>`_' \
		-C '`BROKEN: {cset} <https://github.com/biocommons/hgvs/commit/{cset}>`_' \
		-M '`milestone {mstone} <https://github.com/biocommons/hgvs/milestones/{mstone}>`_' \
		<$< >$@.tmp
	mv $@.tmp $@



.PHONY: clean cleaner cleanest

clean:
	/bin/rm -f *~

cleaner: clean
	#/bin/rm -f *.rst

cleanest: cleaner
	/bin/rm -f *.clog
