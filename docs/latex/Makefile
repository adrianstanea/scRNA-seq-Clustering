# =============================================================================
# Makefile - LaTeX Build Automation
# scRNA-seq Clustering Comparative Study
# =============================================================================

MAIN = main
LATEX = pdflatex
BIBTEX = bibtex

# Default target
all: $(MAIN).pdf

# Main PDF compilation
$(MAIN).pdf: $(MAIN).tex preamble.tex sections/*.tex bibliography/references.bib
	$(LATEX) $(MAIN)
	$(BIBTEX) $(MAIN)
	$(LATEX) $(MAIN)
	$(LATEX) $(MAIN)

# Quick compile (no bibliography update)
quick:
	$(LATEX) $(MAIN)

# Clean auxiliary files
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz

# Clean everything including PDF
distclean: clean
	rm -f $(MAIN).pdf

# View PDF (Linux)
view: $(MAIN).pdf
	xdg-open $(MAIN).pdf &

# Continuous compilation (requires latexmk)
watch:
	latexmk -pdf -pvc $(MAIN)

.PHONY: all quick clean distclean view watch
