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

XY:=1.2
XY_RST:=$(addsuffix .rst,${XY})
XY_VTAGS:=$(filter ${XY}%,${VTAGS})
XY_VTAGS_RST:=$(addsuffix .rst,${XY_VTAGS})

default:
	@echo "no $@ target"; exit 1

all: ${XY_VTAGS_RST} ${XY_RST}


index.rst ${XY_RST}: %.rst:
	clogger-group-XY $* >$@.tmp
	mv $@.tmp $@


hg-git-remap.pl: ../../../.hg/git-mapfile
	perl -ane 'printf("s/%.12s/%.7s/g;\n", $$F[1], $$F[0])' <$< | sort >$@


next.clog:: hg-git-remap.pl
	clogger . | perl -p0  hg-git-remap.pl >$@

.PRECIOUS: %.clog
%.clog:: # hg-git-remap.pl
	clogger $* | perl -p0  hg-git-remap.pl >$@


# TODO: use git-mapfile to map commits
%.rst: %.clog
	clogger-fmt \
		-I '`#{issue_id} <https://github.com/biocommons/hgvs/issues/{issue_id}/>`_' \
		-C '`{cset} <https://github.com/biocommons/hgvs/commit/{cset}>`_' \
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
