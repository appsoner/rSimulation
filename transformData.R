x = matrix(rnorm(600), nrow = 200)
beta = c(1,2,3)
y = x%*%beta + rnorm(100)
ft <- rq(y~x[,1] + x[, 2] + x[,3], tau = 0.5)
ya = cbind(1,x)%*%ft$coefficients
mean(abs(ya - x%*%beta))
h = 2*200^(-1/5)
fv = qlprq2(x[,1], x[,2], x[,3], y, h, 0.5)
mean(abs(fv - x%*%beta)^2)
library(rqPen)
c(1, x[1, ])%*%c(0.008511246, 0.813443862, 0.197095447, 0.867598401)
x[1,]%*%beta
