
source('DoMatrix.R')
source('ExactTableApproach.R')


ngen <- 80

# names(par)<-c(pop, theta.a, phi.a, theta.i, phi.i, theta.n, phi.n, theta.b, phi.b, h, d.m)
par <-        c(713,6,      .69,    1,      .5,     1,      .5,     100,     0,    .05,.57)
# par <-          c(1500,)
aprx <- vector(length = ngen)
pop <- vector(length = ngen)
for(n in 1:ngen){
  temp <- DoMatrixApprox(par = par,
                         n = n)
  aprx[n] <- temp[1]
  pop[n] <- sum(temp[2:4])
}
pp <- pop


# par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)
par <-   c(.69,   .5,    1,      6,       .5,     1,      .57,  0,    100,    .05)
pop = 713 
tol = 10^-3 
num = 100


results <- array(dim = c(ngen,3))
colnames(results) <- c('freq', 'time','bins')
rownames(results) <- c(paste('gen =', 1:ngen))
tables <- MakeSeedTable(ngen = ngen, par = par, pop = pop)
for(n in ngen:1){
  start <- Sys.time()
  cat(n,'\n')
  exact <- DoExact(n = n, 
                   par = par,
                   pop = pop, 
                   tol = tol, 
                   num = num,
                   tables = tables)
  bins <- exact[[1]]
  freq <- exact[[2]]
  tables <- exact[[3]]
  end <- Sys.time()
  times <- end - start
  results[n,] <- c(freq,times,bins)
}







plot(aprx, type = 'l', xlab = 'Generation', ylab = 'frequency', lwd = 2,
     main = 'Frequency Over Time')
legend(x = 'topleft', lwd = 2, legend = c('Approx','Exact'), col = c('black','red'),
       bty = 'n')
lines(results[,1], col = 'red', lwd = 2)

secderiv.array <- FindDerivs(tables$dPsi.array)

sorted.dPsi.array <- SortTable(tables$dPsi.array)
svals <- sorted.dPsi.array[[2]]
sorted.dPsi.array <- sorted.dPsi.array[[1]]

plot(svals,sorted.dPsi.array[ngen,], type = 'l')
lines(svals, secderiv.array[ngen,], lwd = .5, col = 'blue')

var.part.1 <- CompositeSimpsonsWrapper(Psi.array = SortTable(tables$dPsi.array)[[1]], 
                                       deriv.array = sorted.dPsi.array, 
                                       secderiv.array, 
                                       svals)

var <- var.part.1 - (results[,1])^2

upper95 <- results[,1] + 2 * sqrt(var/pp)
lower95 <- results[,1] - 2 * sqrt(var/pp)


plot(results[,1], lwd = 2, type = 'l', xlab = 'Generation', ylab = 'Frequency', main = 'Frequency over Generations',
     ylim = c(0,1))
legend( x = 'topleft', lwd = 2, legend = c('Expectation','95% CI'), col = c('black','red'), bty = 'n')
lines(upper95, lwd = .5, col = 'red')
lines(lower95, lwd = .5, col = 'red')
abline(h = .2)






# plot(pp, type = 'l', xlab = 'Generation', ylab = 'Population Size', lwd = 2,
#      main = 'Population Size Over Time')
# 
# 
# plot(results[,2], type = 'l', xlab = 'Generation', ylab = 'Runtime', lwd = 2,
#      main = 'Runtime Over Time')

