
####################################################################################
##
## Objective: obtain summary of in/out-regulated modules, and plot the GRN
## Programmer: Jaejoon Song
## Date: July 26, 2016
##
####################################################################################


args <- commandArgs(trailingOnly = TRUE)

# install.packages("qgraph")
library(qgraph)

## Reads the matrix
easMat <- read.csv(args[1],header=FALSE)

## Just assigning module number to the matrix
mNumber <- 1:dim(easMat)[1]
moduleNumber <- paste("M",mNumber,sep="")
colnames(easMat) <- moduleNumber
rownames(easMat) <- moduleNumber

## Plots the GRN
pdf(args[2])
mygraph<- qgraph(t(easMat),details=FALSE,shape="square",repulsion=.7, usePCH=TRUE,fade=FALSE,curveScale=TRUE, layout = "spring", node.height=3.5,label.prop=1.1,asize=1)
title("Gene regulatory network",line=2.5)
dev.off()
