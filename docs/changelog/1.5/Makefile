.PHONY: FORCE
.DELETE_ON_ERROR:

# N.B. this scripts requires tools that are not publicly available
# yet. Eventually, clogger will be released.
# Sorry.
PATH:=/home/reece/projects/reece/clogger/bin:${PATH}
SHELL:=/bin/bash -o pipefail

default:
	@echo "no $@ target"; exit 1

next.clog::
	biocommons-changelog . >$@

%.rst: %.clog
	clogger-fmt \
		-I '`#{issue_id} <https://github.com/biocommons/hgvs/issues/{issue_id}/>`_' \
		-C '`{cset} <https://github.com/biocommons/hgvs/commit/{cset}>`_' \
		<$< >$@.tmp
	mv $@.tmp $@



.PHONY: clean cleaner cleanest

clean:
	/bin/rm -f *~

cleaner: clean
	#/bin/rm -f *.rst

cleanest: cleaner
	/bin/rm -f *.clog
