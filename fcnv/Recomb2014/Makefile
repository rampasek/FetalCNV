all: main.pdf

clean:
	-rm main.aux main.bbl main.blg main.log main.log main.lot main.out \
	main.dvi main.ps main.pdf *~

main.ps: main.dvi
	dvips main -o main.ps

main.pdf: main.tex *.tex *.bib
	pdflatex main
	pdflatex main
	-bibtex main
	pdflatex main
	pdflatex main

main.dvi: main.tex *.tex *.bib
	latex main
	latex main
	-bibtex main
	latex main
	latex main

.PHONY: all bib
