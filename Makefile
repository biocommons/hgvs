.SUFFIXES :
.PRECIOUS :
.PHONY : FORCE
.DELETE_ON_ERROR:

SHELL:=/bin/bash -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))
export PYTHONPATH=lib/python


############################################################################
#= BASIC USAGE
default: help

#=> help -- display this help message
help:
	@sbin/extract-makefile-documentation "${SELF}"


############################################################################
#= UTILITY FUNCTIONS

#=> lint -- run lint
# TBD

#=> test -- run tests
test:
	PYTHONPATH=lib/python python setup.py nosetests -v --with-xunit

#=> docs -- make sphinx docs
docs: build_sphinx

#=> develop, build_sphinx, sdist, upload_sphinx
develop build build_sphinx install sdist upload_sphinx: %:
	python setup.py $*

#=> upload-<tag>
upload-%:
	hg up -r $*
	python setup.py sdist upload

invitae-upload-%:
	hg up -r $*
	python setup.py sdist upload -r invitae


############################################################################
#= CLEANUP
.PHONY: clean cleaner cleanest pristine
#=> clean: clean up editor backups, etc.
clean:
	find . -name \*~ -print0 | xargs -0r /bin/rm
#=> cleaner: above, and remove generated files
cleaner: clean
	find . -name \*.pyc -print0 | xargs -0r /bin/rm -f
	/bin/rm -fr distribute-* *.egg *.egg-info
	make -C doc clean
#=> cleanest: above, and remove the virtualenv, .orig, and .bak files
cleanest: cleaner
	find . \( -name \*.orig -o -name \*.bak \) -print0 | xargs -0r /bin/rm -v
	/bin/rm -fr build ve dist bdist
#=> pristine: above, and delete anything unknown to mercurial
pristine: cleanest
	hg st -un0 | xargs -0r echo /bin/rm -fv
