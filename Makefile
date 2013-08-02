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
#= INSTALLATION/SETUP

#=> setup -- prepare python and perl environment for prequisites
#=>> This is optional; the only requirement is that packages are discoverable
#=>> in PYTHONPATH and PERL5LIB
setup: setup-python

#=> setup-python: create a virtualenv with base packages
# NOTE: setup-python only makes the virtualenv. You must actvate it
# yourself (source ve/bin/activate)
setup-python: ve
	source ve/bin/activate; python setup.py develop
ve: virtualenv.py
	python $< --distribute ve
virtualenv.py:
	curl https://raw.github.com/pypa/virtualenv/master/virtualenv.py >$@


############################################################################
#= UTILITY FUNCTIONS

#=> lint -- run lint


#=> test -- run tests
test:
	PYTHONPATH=lib/python python setup.py nosetests -v --with-xunit

#=> docs -- make sphinx docs
docs: build_sphinx

#=> develop, build_sphinx, sdist, upload_sphinx
develop build_sphinx install sdist upload_sphinx: %:
	python setup.py $*

#=> upload-<tag>
upload-%:
	hg up -r $*
	python setup.py sdist upload



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
	/bin/rm -fr ve dist bdist
#=> pristine: above, and delete anything unknown to mercurial
pristine: cleanest
	hg st -un0 | xargs -0r echo /bin/rm -fv
