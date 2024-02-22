#!/bin/bash

for texfile in `ls *.tex` ; do
    echo $texfile
    cat $texfile | sed 's#{figures/#{#g' > toto
    mv toto $texfile

done

cat <<EOF > Makefile
all : 	paper.pdf

paper.pdf : *.tex *.png Makefile
	pdflatex -halt-on-error paper
	pdflatex -halt-on-error paper
	pdflatex -halt-on-error paper

clean:
	@rm -f *.log *.lof *.aux *.out paper.pdf *.dvi *.blg *.mtc *.mtc0 *.toc *.maf *.idx
EOF
