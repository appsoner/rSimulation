num = 5
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
xData <- cbind(x1Vector, tVector, alphaVector)
xChengData <- cbind(x1Vector, tVector, alphaVector)
dataX <- array(rep(0,nrow(xMatrix)), dim = c(5, 8, nrow(xMatrix)/5))
dataY <- array(rep(0,nrow(xMatrix)), dim = c(5, 1, nrow(xMatrix)/5))
dataG <- array(rep(0,nrow(xMatrix)), dim = c(5, 1, nrow(xMatrix)/5))
dataF <- array(rep(0,nrow(xMatrix)), dim = c(5, 1, nrow(xMatrix)/5))
dataAlpha <- array(rep(0,nrow(xMatrix)), dim = c(5, 1, nrow(xMatrix)/5))

g <- exp(x1Vector)
f <- sin(2*pi*tVector)
for(i in seq(1, nrow(xData), 5)){
  dataX[ , , (i+4)/5] <- xMatrix[i:(i+4), ]
  dataY[ , , (i+4)/5] <- yVector[i:(i+4)]
  dataG[ , , (i+4)/5] <- g[i:(i+4)]
  dataF[ , , (i+4)/5] <- f[i:(i+4)]
  dataAlpha[ , , (i+4)/5] <- alphaVector[i:(i+4)]
  #print(xData[i:(i+4), ])
 # print("################")
}
