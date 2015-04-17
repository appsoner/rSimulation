######################estimate efficients Beta#####################
rho <- function(tau, u){
  ifelse(u<0, u*(tau-1), u*tau)
}
beta0 <- c(1/sqrt(2), 0 ,1/sqrt(2) ,0 , 0, 0, 0, 0)
h = 1
#################Kernel########################################
wKernel <- function(x1, x2, h, n){
  return (1/n/h*exp(-(x1-x2)^2/2/h^2)/sqrt(2*pi))
}
######################W_ij,ls###################################
wIjLs <- function(i,j,l,s, h, n, beta){
  sum = 0
  wH2 <- wKernel(dataX[j, , i]%*%beta, dataX[s, ,l]%*%beta, h, n)
  for(l1 in 1:n){
    for(s1 in 1:5){
      sum = sum + wKernel(dataX[s1, , l1]%*%beta, dataX[s, , l]%*%beta, h, n)
    }
  }
  return (wH2/sum)
}
##############step 1 function#######################################
step1 <- function(gg){
  sum = 0
  for(i in 1:dim(dataX)[3]){
    for(j in 1:5){
          u <- dataY[j, , i]-dataG[s, , l]-
            gg*(dataX[j, , i]%*%beta0-dataX[s, , l]%*%beta0)-
            dataF[j, , i] - dataAlpha[j, , i]
          sum = sum + rho(0.5, u)*wIjLs(i,j,l,s,h,dim(dataX)[3], beta0)
    }
  }
  return (sum)
}
#############################################################
####################################################
l = 1
s = 2
optim(1, step1, method = "Brent", lower = 0, upper = 12)
#################################################################
matrixMul <- function(data, beta){
  sum = 0
  for(i in 1:length(data)){
    sum = sum + data[i]*beta[i]
  }
  return (sum)
}
######################step2###################################
count = 1
step2 <- function(beta){
  sum = 0
  for(i in 1:dim(dataX)[3]){
    for(j in 1:5){
      for(l in 1:dim(dataX)[3]){
        for(s in 1:5){
          u <- dataY[j, , i]-dataG[s, , l]-
            dataDrivate[s, , l]*(dataX[j, , i]%*%beta0-matrixMul(dataX[s, , l], beta))-
            dataF[j, , i] - dataAlpha[j, , i]
          sum = sum + rho(0.5, u)*wIjLs(i,j,l,s,h,dim(dataX)[3], beta)
         # print(sum)
       #  print(count)
         count = count + 1
        }
      }
    }
  }
  print(sum)
  print(count)
  print(beta)
  return (sum) 
  
}
################################################################
dataDrivate <- dataG
optim(beta0, step2)
##################Simulation#####################################
beta0 <- c(1/sqrt(3), 0, 1/sqrt(3), 0, 0, 0, 0, 0)
std = 1
while(std > 0.01){
  gg <- NULL
  for(l in 1:dim(dataX)[3]){
    for(s in 1:5){
      g <- optim(1, step1, method = "Brent", lower = 0, upper = 12)$par
      gg <- rbind(gg, g)
      print(gg)
    }
  }
  dataDrivate <- array(rep(0,nrow(xMatrix)), dim = c(5, 1, nrow(xMatrix)/5))
  for(i in seq(1, nrow(xData), 5)){
    dataDrivate[, , (i+4)/5] <- gg[i:i+4]
  }
  
  betaRet <- optim(beta0, step2)$par
  std <- sum((beta0 - betaRet)^2)
  beta0 <- betaRet
}
