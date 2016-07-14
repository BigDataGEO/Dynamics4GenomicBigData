#!/bin/bash

directoryPath=$1

cd $directoryPath

pdflatex Paper_dum.tex

bibtex Paper_dum.aux

pdflatex Paper_dum.tex

pdflatex Paper_dum.tex

rm Paper_dum.aux Paper_dum.bbl Paper_dum.blg Paper_dum.log Paper_dum.out