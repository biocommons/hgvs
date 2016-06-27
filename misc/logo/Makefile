.SUFFIXES:
.PHONY: FORCE
.DELETE_ON_ERROR:


all: hgvs-logo.png hgvs-logo.ico favicon.ico

%.png: %.svg
	convert $< $@

%.ico: %.svg
	convert -geometry 32x32 $< $@

favicon.ico: hgvs-logo.ico
	ln -s $< $@
