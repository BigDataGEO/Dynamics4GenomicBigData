########################################################################################################################
#input
#x: mean curve 1
#y: mean curve 2
library(fdrtool)
TestTwoCurves <- function(x, y){
  n <- length(x)
  xy.cor <- cor(x, y, method = "spearman")
  p.value.cor <- pcor0(xy.cor, kappa = n - 1, lower.tail = TRUE)
  p.value.median <- wilcox.test(x, y, method = "two.sided", paired = TRUE)$p.value
  print(c(p.value.cor, p.value.median))
  return(min(p.adjust(c(p.value.cor, p.value.median), method="bonferroni")))
}

#x and y are two mean curves for one module. The function TestTwoCurves will return a p-value of testing their association.
#Example:
#x <- rnorm(8)
#y <- rnorm(8)
#TestTwoCurves(x,y)

args <- commandArgs(TRUE);

filename = args[1];

# time_points_filename = args[2];
 
curves=read.csv(filename, header = FALSE, sep = ",");

colnames(curves) = c('x', 'y');
 
x=curves$x;
 
y=curves$y;

# t=as.numeric(unlist(read.csv(time_points_filename, header = FALSE, sep = ",")));

TestTwoCurves(x, y);
