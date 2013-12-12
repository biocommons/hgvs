.SUFFIXES :
.PRECIOUS :
.PHONY : FORCE
.DELETE_ON_ERROR:

# for ve
PYTHON_VERSION=2.7
SYSTEM_PYTHON:=/usr/bin/python${PYTHON_VERSION}
VE_VERSION=1.10.1
VE_DIR=ve

#SHELL:=/bin/bash -o pipefail
SELF:=$(firstword $(MAKEFILE_LIST))
export PYTHONPATH=lib/python

PYPI_SERVICE:=-r invitae

ifdef LOCAL_UTA
export UTA_DB_URL=postgresql://localhost/
endif

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


# virtualenv
VE_PY_DIR:=virtualenv-${VE_VERSION}
VE_PY:=${VE_PY_DIR}/virtualenv.py

build_ve:
	curl -sO  https://pypi.python.org/packages/source/v/virtualenv/virtualenv-${VE_VERSION}.tar.gz
	tar xvfz virtualenv-${VE_VERSION}.tar.gz
	rm -f virtualenv-${VE_VERSION}.tar.gz
	${SYSTEM_PYTHON} ${VE_PY} ${VE_DIR} 2>$@.log | tee "$@.err"

#=> test -- run tests
test:
	PYTHONPATH=lib/python python setup.py nosetests -v --with-xunit

#=> docs -- make sphinx docs
docs: build_sphinx

#=> develop, build_sphinx, sdist, upload_sphinx
develop bdist bdist_egg build build_sphinx install sdist upload_sphinx: %:
	python setup.py $*

setup_dev:
	pip install -r requirements_dev.txt
	python setup.py develop

#=> upload-<tag>
upload-%:
	[ -z "$$(hg st -admnr)" ] || { echo "Directory contains changes; aborting." 1>&2; hg st -admr; exit 1; }
	R=$$(hg id -t); hg up -r $*; python setup.py sdist upload ${PYPI_SERVICE}; hg up -r $$R


############################################################################
#= CLEANUP
.PHONY: clean cleaner cleanest pristine
#=> clean: clean up editor backups, etc.
clean:
	find . \( -name \*~ -o -name \*.fai \) -print0 | xargs -0r /bin/rm
#=> cleaner: above, and remove generated files
cleaner: clean
	find . -name \*.pyc -print0 | xargs -0r /bin/rm -f
	/bin/rm -fr distribute-* *.egg *.egg-info nosetests.xml
	make -C doc clean
	/bin/rm -rf virtualenv-${VE_VERSION}* build_ve.*
#=> cleanest: above, and remove the virtualenv, .orig, and .bak files
cleanest: cleaner
	find . \( -name \*.orig -o -name \*.bak \) -print0 | xargs -0r /bin/rm -v
	/bin/rm -fr build ve dist bdist
	/bin/rm -fr ${VE_DIR}
#=> pristine: above, and delete anything unknown to mercurial
pristine: cleanest
	hg st -un0 | xargs -0r echo /bin/rm -fv
