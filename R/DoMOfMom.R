data <- read.csv(file = 'data/offspring info.csv')
data <- data[!is.na(data[,1]),]

xbar <- mean(data$Calves.sired)
ssqr <- var(data$Calves.sired)

if(xbar >= ssqr){
  phihat <- 0
  thetahat <- xbar
}else{
  phihat <- (ssqr - xbar) / (xbar^2 + ssqr - xbar)
  thetahat <- xbar + (ssqr / xbar) - 1
}
phihat
thetahat
