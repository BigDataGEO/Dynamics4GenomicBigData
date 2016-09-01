#!/bin/bash

directoryPath=$1

cd $directoryPath

pdflatex Paper.tex

bibtex Paper.aux

pdflatex Paper.tex

pdflatex Paper.tex

rm Paper.aux Paper.bbl Paper.blg Paper.log Paper.out