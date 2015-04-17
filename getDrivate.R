##########################################
###estimate the derivate of function f####
##########################################
#-----------------------------------------
#Description:  A method depending interpolation 
#              to estimate derivate of function 
#----------------------------------------------
#Input
#x is one demension covariate
#y is the value of f(x)
#-----------------------------------------------
#Return 
#list: A list containing the covariate x and 
#      the derivate of f(x) at points x.
#-----------------------------------------------
#Author:         Guosheng Xu
#Insititution:   Zhejinag University of Fiance & Economics
#Date:           26/03/2015
#-----------------------------------------------
###########################
###the function getDeriv###
###########################
library(lokern)
library(pracma)
getDeriv <- function(x, y){
  fit1 <- glkerns(x, y)
  f <- gradient(fit1$est, fit1$x.out)
  fit2 <- glkerns(fit1$x.out, f, x.out = x)
  return (list(x = fit2$x.out, deriv = fit2$est))  
}


######################
###a toy simulation###
######################
#----------------------
#Model: y = sin(x)
#----------------------
set.seed(150326)
x = rnorm(500)
y = sin(x)
retVal = getDeriv(x,y)
temp = order(x)
plot(x[temp], y[temp], type = "l", lwd = 2, col = "red")
lines(x[temp], retVal[[2]], lwd = 2, col = "blue")
abline(h = 0, lty = 2, lwd = 2)


