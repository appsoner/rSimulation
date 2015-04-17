##########S_u^x##########
suxElement <- function(x1, x2, h, n){
  return (1/n/h*exp(-(x1-x2)^2/2/h^2)/sqrt(2*pi))
}
sux <- function(x, u, h){
  suxMatrix = NULL
  xx <- x[ , u]
  vector <- rep(0, nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:nrow(x)){
      vector[j] <- suxElement(xx[i], xx[j], h, nrow(x))
      #print(vector)
    }
    suxMatrix = rbind(suxMatrix, vector)
  }
  return (suxMatrix)
}
############################################
demo <- matrix(1:12, nrow = 3)
sux(demo, 1, 0.2)
###########s_u^w##########################
distance <- function(x1, x2){
  return (mean(abs(x1-x2)))
}
suwElement <- function(x1, x2, h, n){
  return (1/n/h^length(x1)*exp(-distance(x1,x2)/2/h^2)/sqrt(2*pi))
}
suw <- function(x, u, h){
  suwMatrix <- NULL
  xx <- x[ ,-u]
  vector <- rep(0, nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:nrow(x)){
      vector[j] <- suwElement(xx[i, ], xx[j, ], h, nrow(x))
      #print(xx[i, ])
    }
    suwMatrix <- rbind(suwMatrix, vector)
  }
  return (suwMatrix)
}
#######################################
demo <- matrix(1:12, nrow = 3)
suw(demo, 2, 0.2)
############################################
####################S###################
distance <- function(x1, x2){
  return (mean(abs(x1-x2)))
}
sElements <- function(x1, x2, h, n){
  return (1/n/h^length(x1)*exp(-distance(x1,x2)/2/h^2)/sqrt(2*pi))
}
s <- function(x, h){
  sMatrix <- NULL
  vector <- rep(0, nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:nrow(x)){
      vector[j] <- sElements(x[i, ], x[j, ], h, nrow(x))
    }
    sMatrix <- rbind(sMatrix, vector)
  }
  return (sMatrix)
}
#########################################
demo <- matrix(1:12, nrow = 3)
s(demo, 0.2)
###########################################
###########q_hat###########################
library(quantreg)
qlprq <- function(x1, x2, x3, y, h, tau = .5)
{  
 # if(x2 == NULL)
#    x2 <- 0
  fv <- x1
  for(i in 1:length(x1)) {
    z1 <- x1 - x1[i]
    z2 <- x2 - x2[i]
    z3 <- x3 - x3[i]
    dist <- sqrt(z1^2+z2^2+z3^2)
    wx <- dnorm(dist/h)
    r <- rq(y~z1+z2+z3, weights=wx, tau=tau, ci=FALSE)
    fv[i] <- r$coef[1.]
  }
  return (fv)
}
qlprq2 <- function(x1, x2, x3, y, h, tau = .5)
{  
  # if(x2 == NULL)
  #    x2 <- 0
    fv = x1
    for(i in 1:length(x1)){
      z1 <- x1 - x1[i]
      z2 <- x2 - x2[i]
      z3 <- x3 - x3[i]
      dist <- apply(cbind(abs(z1), abs(z2), abs(z3)), 1, mean)
      wx <- dnorm(dist/h)
      r <- rq(y~x1+x2+x3, weight = wx, tau=tau, ci=FALSE)
      print(r$coefficients)
      fv[i] <- cbind(1,x1[i], x2[i], x3[i])%*%r$coefficients
    }
  return (fv)
}
#qlprq3 <- function(x1, y, h, tau = .5, m = 50)
#{  
#  fv <- x1
#  for(i in 1:length(x1)) {
#    z1 <- x1 - x1[i]
#    dist <- sqrt(z1^2)
#    wx <- dnorm(dist/h)
#    r <- rq(y~z1, weights=wx, tau=tau, ci=FALSE)
#    fv[i] <- r$coef[1.]
#  }
#  list(x1 = x1, fv = fv)
#}
###############################################
#x1 <- rnorm(10)
#x2 <- runif(10)
#x3 <- rnorm(10)
#y <- x1 + x2 + x3
#quantile(x2, probs=0.5)
#solu <- qlprq3(x2,x2,1)
#solu <- qlprq2(x1, x2, x3, y, 1)
#solu$fv
#mean(solu$fv)
###############################################
###########q_{u,\alpha}^\hat###################
qAlphaTime <- function(x, y, u, h){
  e <- rep(1, nrow(x))
  part <- qlprq(x[ ,1], x[, 2], x[ ,3] , y, h)*(suw(x, u, h)%*%e)/(s(x,h)%*%e)
  return (sux(x, u, h)%*%part)
}
#################################################
#x <- matrix(rnorm(18), ncol = 3)
#y <- 1:6
#ans <- qAlphaTime(x, y, 1, 1)
#################################################
#################c_/alpha##########################
cAlpha <- function(x, y, h){
  return (mean(qlprq(x[ ,1], x[ , 2], x[ ,3], y, h)))
}
##################################################
#x <- matrix(rnorm(18), ncol = 3)
#y <- 1:6
#cAlpha(x, y, 1)
##################################################
#############q_{u,\alpha}^hat#####################
qAlphaHat <- function(x, y, u, h){
  return (qAlphaTime(x, y, u, h) - cAlpha(x, y, h))
}
######################################################
#x <- matrix(rnorm(12), ncol = 2)
#y <- 1:6
#qAlphaHat(x, y, 1, 1)
########################Oracle estimatot###############################
qOracle <- function(x, y, u, h, tau = .5){
  qTime = 0
  set <- c(1,2,3)[-u]
  for(i in set){
    qTime <- qTime +  qAlphaTime(x, y, i, h)  
  }
  cA <- cAlpha(x, y, h)
  y <- y + cA - qTime
  return (qlprq(x[, 1], x[ , 2], x[ , 3], y, h))
}
#####################################################################
#x <- matrix(rnorm(18), ncol = 3)
#y <- 1:6
#ans <- qOracle(x, y, 1, 1)
#########################generate random variables#####################
num = 100
x1Vector <- NULL;
xMatrix <- NULL;
tVector <- NULL;
yVector <- NULL;
alphaVector <- NULL;
for(i in 1:num){
  beta <- c(1/sqrt(2), 0, 1/sqrt(2), 0, 0, 0, 0, 0)
  alpha <- rnorm(1)
  for(j in 1:5){
    xline <- rnorm(8)
    xMatrix <- rbind(xMatrix, xline)
    t <- runif(1)
    tVector <- rbind(tVector, t)
    x1Line <- xline%*%beta
    x1Vector <- rbind(x1Vector, x1Line)
    y <- exp(x1Line) + sin(2*pi*t) + alpha + rnorm(1)
    yVector <- rbind(yVector, y)
    alphaVector <- rbind(alphaVector, alpha)
  }  
}
##################Simulation##############################


g <- exp(x1Vector)
f <- sin(2*pi*tVector)
x <- cbind(x1Vector, tVector, alphaVector)
y <- yVector
h = sd(x1Vector)+sd(tVector)+sd(alphaVector)
h = h/3
gS <- qOracle(x, y, 2, h, 0.5)
mean(abs(f - gS))
###########################
gS <- qAlphaHat(x, y, 0.5, h)
mean(abs(g - gS))
