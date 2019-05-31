# Makefile to convert example notebooks into rst files to be included
# in docs.
#
# To update documentation:
# 1. Add/edit notebooks in this directory (hgvs/examples/)
# 2. Type `make install` to convert and copy to the documentation
#    examples directory (hgvs/doc/examples)
# 3. In the top level directory (hgvs/), type `make docs`. Verify that
#    the examples are rendered as intended.

.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

SHELL:=/bin/bash

DOC_EXAMPLES_DIR:=../doc/examples

nbs=$(subst .ipynb,,$(wildcard *.ipynb))
rst_nbs:=$(addsuffix .rst,${nbs})

default: install

install: $(addprefix ${DOC_EXAMPLES_DIR}/,${rst_nbs})

${DOC_EXAMPLES_DIR}/%.rst: %.ipynb
	ipython nbconvert --to rst --output=$@ $<

.PHONY: clean cleaner cleanest
clean:
	/bin/rm -f *~
cleaner: clean
cleanest: cleaner
	/bin/rm -fr *_files *.rst
