all:
	pdflatex report.tex
clean:
	rm *.aux
	rm *.bbl
	rm *.bcf
	rm *.blg
	rm *.log
	rm *.run.xml

