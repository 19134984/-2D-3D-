###########################################################################
## OpenLB user guide (Quarto)
##
## The user guide is rendered twice from the same sources: once as a PDF
## book (_quarto-book.yml) and once as a website (_quarto.yml). Quarto
## applies one profile per invocation, so the two renders cannot be merged
## into a single command. The order matters, because the website only shows
## the PDF download button in its navbar when the PDF is already on disk --
## see ug_pdf_download.sh. Plain `make` therefore renders the PDF first and
## the website second.
##
## The targets are phony: Quarto tracks freshness itself, and the names of
## the generated files are set in the Quarto configs rather than here.

QUARTO ?= quarto

# The release is kept in OLB_ROOT/VERSION; the guide shows the series only,
# so 1.9r0 becomes 1.9.
OLB_RELEASE := $(shell cat ../VERSION)
OLB_SERIES  := $(firstword $(subst r, ,$(OLB_RELEASE)))

# Both renders share the .quarto state directory, so they must not run at
# the same time. With .NOTPARALLEL the prerequisites of `all` are also made
# in the order they are listed, even under `make -j`.
.NOTPARALLEL:

.PHONY: all pdf website prepare clean

all: pdf website

# Quarto builds its list of inputs -- for the book, the chapters and the
# example include list they pull in -- before it runs its own pre-render
# script, so on a fresh checkout the generated example files have to be in
# place beforehand. link_examples.sh is idempotent and Quarto runs it again
# as its pre-render step.
prepare:
	./link_examples.sh
	./example_table.py

# doc/ug_pdf/<name>.pdf
pdf: prepare
	$(QUARTO) render --profile book --to pdf -M "subtitle=Associated with Release $(OLB_SERIES) of the Code"

# doc/ug/index.html -- picks up a PDF rendered beforehand
website: prepare
	$(QUARTO) render

clean:
	rm -rf ug ug_pdf .quarto
	rm -f include/pdf_download.html examples/_includes.qmd examples/_table.qmd examples/_table_pdf.qmd
	find examples -type l -name '*.qmd' -delete
