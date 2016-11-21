#!/bin/bash

directoryPath=$1

cd $directoryPath

pdflatex Report.tex

bibtex Report.aux

pdflatex Report.tex

pdflatex Report.tex

rm Report.aux Report.bbl Report.blg Report.log Report.out