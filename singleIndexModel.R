#######################################################
#######Single Index Model Quantile Regression##########
#######################################################
#-----------------------------------------------------
#Description 
#Main function : qrSIM, which is depend on lprq function. 
#other function: lprq, which dervied from the source code
#                of function lprq in quantReg package. 
#-----------------------------------------------------
#------------------------------------------------------
#Input
#y                n*1 response vector
#xx               n*p covariates matrix
#theta0           p*1 the initial value of parameter theta0
#                 (also be beta sometimes)
#tau              scalar quantile level
#maxit            maxmum times for iteration
#crit             ctiterion for convergence of theta
#------------------------------------------------------
#Return
#thetaNew         the estimation of the theta
#-----------------------------------------------------------
#Author           Guosheng Xu, 
#Insititution     Zhejinag University of Fiance & Economics
#------------------------------------------------------------
library(pracma)
library(KernSmooth)
library(quantreg)  
library(KernSmooth)
lprq = function (x, y, h, tau = 0.5,x0)
{
  fv <- x0
  dv <- x0    
  z <- x - x0
  wx <- dnorm(z/h)
  r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
  fv <- r$coef[1]
  dv <- r$coef[2]
  list(x0 = x0, fv = fv, dv = dv)
}
qrSIM <- function(y,xx, tau, theta0, maxit, crit){
  thetaNew = theta0*sign(theta0[1])/sqrt(sum(theta0^2))
  n = nrow(xx)
  d = ncol(xx)
  a = matrix(rep(0,2*n), nrow = n)
  iter = 0
  thetaOld = 2*thetaNew
  beta = 0
  betaNew = beta
  while((iter<maxit) && (sum((thetaNew-thetaOld)^2)>crit)){
    thetaOld = thetaNew
    iter = iter + 1
    t = xx%*%thetaOld
    ystar = y
    ordert = sort(t)
    ordery = sort(y)
    hm = dpill(ordert[20:(n-20)], ordery[20:(n-20)]) # this value can be adjusted
    hp = hm*(tau*(1-tau)/(dnorm(qnorm(tau)))^2)^.2
    h = hm
    for (i in 1:n) {
      fit = lprq(t, y, h, tau, t[i])
      a[i, 1] = fit$fv
      a[i, 2] = fit$dv
    }
    xnew = zeros(n^2, d)
    amean = a[, 1]
    ynew = repmat(y,n,1) - kron(amean, ones(n, 1))
    for (j in 1:n) {
      for (i in 1:n) {
        xnew[(j-1)*n+i, ] = a[j, 2]*(xx[i, ] - xx[j, ])
      }
    }
    xg = xnew[ , 1:d]%*%thetaOld
    xgh = xg/h
    wt = dnorm(xgh)*abs(xgh<3)
    ind = which(wt != 0)
    wt = wt[ind]
    ynew = ynew[ind]
    xnew = xnew[ind, ]
    fit <- rq(ynew ~ xnew, tau, weights = wt)
    thetaNew = fit$coef[2:(d+1)]
    thetaNew = sign(thetaNew[1])*thetaNew/sqrt(sum(thetaNew^2))
  }
  return (thetaNew)
}
#######################################################
#######A simulation  Y = exp(x*beta)+random term######
#######################################################
xxx <- matrix(rnorm(2400), ncol = 8)
betaa <- c(1/sqrt(2), 0, 1/sqrt(2), 0, 0, 0, 0, 0)
yyy <- xxx%*%betaa
yyy <- exp(yyy)  + rnorm(300)
maxit=20;crit=1e-6;
ans <- qrSIM(yyy, xxx, 0.5, betaa, maxit, crit)
