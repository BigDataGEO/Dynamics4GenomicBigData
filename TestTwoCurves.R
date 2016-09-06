library(fdrtool)

TestTwoCurves <- function(x, y){
  n <- length(x)
  xy.cor <- cor(x, y, method = "spearman")
  p.value.cor <- pcor0(xy.cor, kappa = n - 1, lower.tail = FALSE)
  p.value.median <- wilcox.test(x, y, method = "two.sided", paired = TRUE)$p.value
  return(min(p.adjust(c(p.value.cor, p.value.median), method="bonferroni")))
}

args <- commandArgs(TRUE);

filename = args[1];
 
curves=read.csv(filename, header = FALSE, sep = ",");

colnames(curves) = c('x', 'y');
 
x=curves$x;
 
y=curves$y;

TestTwoCurves(x,y);
