PAPER=report
SLIDES=slides

all: $(PAPER).pdf $(SLIDES).pdf

$(PAPER).pdf: $(PAPER).tex $(PAPER).bib
	pdflatex -interaction=nonstopmode $(PAPER) || true
	bibtex $(PAPER) || true
	pdflatex -interaction=nonstopmode $(PAPER) || true
	pdflatex -interaction=nonstopmode $(PAPER)

$(SLIDES).pdf: $(SLIDES).tex $(PAPER).bib
	-pdflatex -interaction=nonstopmode $(SLIDES)
	-bibtex $(SLIDES)
	-pdflatex -interaction=nonstopmode $(SLIDES)
	-pdflatex -interaction=nonstopmode $(SLIDES)

slides: $(SLIDES).pdf

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.fdb_latexmk *.fls *.synctex.gz *.nav *.snm *.toc *.vrb

distclean: clean
	rm -f $(PAPER).pdf $(SLIDES).pdf

.PHONY: all slides clean distclean
