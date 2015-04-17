rho <- function(tau, u){
  ifelse(u<0, u*(tau-1), u*tau)
}
wIjLs <- function(i,j,dataX, beta,h){
  temp1 <- dnorm((dataX[i,]%*%beta - dataX[j,]%*%beta)/h)
  temp2 <- 0
  for(ic in 1:nrow(dataX)){
    temp2 = temp2 + dnorm((dataX[ic,]%*%beta - dataX[j,]%*%beta)/h)
  }
  return (temp1/temp2)
}
matrixMul <- function(data, beta){
  sum = 0
  for(i in 1:length(data)){
    sum = sum + data[i]*beta[i]
  }
  return (sum)
}
count <- 0
step2 <- function(beta){
  sum = 0
  for(i in 1:nrow(dataX)){
    for(j in 1:nrow(dataX)){
      
          u <- dataY[i]-dataG[j]-
            dataDrivate[j]*(matrixMul(dataX[i,]-dataX[j,], beta))-
            dataF[i] - dataAlpha[i]
         sum = sum + rho(0.5, u)*wIjLs(i,j,dataX,beta,(nrow(dataX))^(-1/5))
         # print(sum)
          #  print(count)
        
      
    }
  }
  count = count + 1
  print(sum)
  print(count)
  #  print(beta)
  return (sum) 
  
}

beta0 <- c(1/sqrt(2), 0, 1/sqrt(2), 0, 0, 0, 0, 0)
std = 1
while(std > 0.01){  
  betaRet <- optim(beta0, step2)$par
  std <- sum((beta0 - betaRet)^2)
  beta0 <- betaRet
}