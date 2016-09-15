#!/bin/bash

directoryPath=$1

cd $directoryPath

pdflatex paper.tex

bibtex paper.aux

pdflatex paper.tex

pdflatex paper.tex

rm paper.aux paper.bbl paper.blg paper.log paper.out