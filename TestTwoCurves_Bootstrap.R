########################################################################################################################
#Bootstrap
#input
#x: mean curve 1
#y: mean curve 2
#t: associated time of x or y
#Rep: replicates
TestTwoCurves <- function(x, y, t, Rep = 10000){
  n <- length(x)
  xy <- c(x,y)
  time.coeff <- c(0, diff(t, 1)/2) + c(diff(t, 1)/2, 0)
  test.stat <- abs(x - y) %*% time.coeff
  test.stat.simu <- rep(0, Rep)
  for(i in 1:Rep){
    x.tmp <- sample(xy, size = n, replace = TRUE)
    y.tmp <- sample(xy, size = n, replace = TRUE)
    test.stat.simu[i] <- abs(x.tmp - y.tmp) %*% time.coeff
  }
  p.value <- length(which(test.stat.simu > as.numeric(test.stat))) / Rep
  return(p.value)
}

#Example
#x <- rnorm(8, mean = 0)
#y <- rnorm(8, mean = 3)
#p <- TestTwoCurves(x, y, t = 1:8, Rep = 10000) 
#p

args <- commandArgs(TRUE);

filename = args[1];

time_points_filename = args[2];
 
curves=read.csv(filename, header = FALSE, sep = ",");

colnames(curves) = c('x', 'y');
 
x=curves$x;
 
y=curves$y;

t=as.numeric(unlist(read.csv(time_points_filename, header = FALSE, sep = ",")));

TestTwoCurves(x, y, t);