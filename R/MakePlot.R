##### Windows #####
library(doSNOW)
library(foreach)
source('DoMatrix.R')
source('DoExact.R')
cl<-makeCluster(7, type="SOCK")
on.exit(stopCluster(cl))
opts <- list(preschedule = FALSE)
registerDoSNOW(cl)
##### Max/Unix #####
library(doMC)
library(foreach)
source('DoMatrix.R')
source('DoExact.R')
#cl<-makeCluster(7, type="SOCK")
#on.exit(stopCluster(cl))
opts <- list(preschedule = FALSE)
registerDoMC(7)

ngen <- 40

# names(par)<-c(pop, theta.a, phi.a, theta.i, phi.i, theta.n, phi.n, theta.b, phi.b, h, d.m)
par <-        c(613,6,      .69,    1,      .5,     1,      .5,     100,     0,    .05,.57)
# par <-          c(1500,)
aprx <- vector(length = ngen)
pop <- vector(length = ngen)
for(n in 1:ngen){
  temp <- DoMatrixApprox(par = par,
                         n = n)
  aprx[n] <- temp[1]
  pop[n] <- sum(temp[2:3])
}
pp <- pop


# par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)
par <-   c(.69,   .5,    1,      6,       .5,     1,      .57,  0,    100,    .05)
pop = 613 
tol = 10^-3 
num = 100

results <- array(dim = c(ngen, 3))
colnames(results) <- c('freq', 'time','bins')
results <- foreach(n = 1:ngen, .options.multicore=opts, .combine = 'rbind') %dopar% {
  start <- Sys.time()
  cat(n,'\n')
  exact <- DoExact(n = n, 
                  par = par,
                  pop = pop, 
                  tol = tol, 
                  num = num)
  bins <- exact[1]
  exact <- exact[2]
  end <- Sys.time()
  times <- end - start
  res <- c(exact, times,bins)
  res
}

exact <- results[,1]
times <- results[,2]
bins <- results[,3]

plot(aprx, type = 'l', xlab = 'Generation', ylab = 'frequency', lwd = 2,
     main = 'Frequency Over Time')
legend(x = 'topleft', lwd = 2, legend = c('Approx','Exact'), col = c('black','red'),
       bty = 'n')
lines(exact, col = 'red', lwd = 2)


plot(pp, type = 'l', xlab = 'Generation', ylab = 'Population Size', lwd = 2,
     main = 'Population Size Over Time')


plot(times, type = 'l', xlab = 'Generation', ylab = 'Runtime', lwd = 2,
     main = 'Runtime Over Time')
#lines(bins, col = 'red', lwd = 2)

plot(times~exact, type = 'l', lwd = 2, main = 'Frequency vs Runtime')
#lines(bins~exact, col = 'red', lwd = 2)

#plot(times~bins, main = 'Bins vs Runtime')
