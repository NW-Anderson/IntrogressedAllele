library(profvis)
source('DoMatrix.R')
source('GetPGF.R')
source('GetExp.R')

ngen <- 5
# names(par)<-c(pop, theta.a, phi.a, theta.i, phi.i, theta.n, phi.n, theta.b, phi.b, h, d.m)
par <-        c(1000,20,      .95,    1,      .1,     1,      .1,     12,      .5,    .05,.5)
aprx <- vector(length = ngen)
pop <- vector(length = ngen)
for(n in 1:ngen){
  temp <- DoMatrixApprox(par = par,
                         n = n)
  aprx[n] <- temp[1]
  pop[n] <- sum(temp[2:3])
}


# par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)
par <-   c(.95,   .1,    1,       20,      .1,    1,       .5,  .5,    12,      .05)
pop = 1000 
tol = 10^-3 
num = 100

# profvis({
#   GetExp(n = 20, 
#          par = par,
#          pop = pop, 
#          tol = tol, 
#          num = num)
# })

exact <- vector(length = ngen)
for(n in 1:ngen){
  cat(n,'\n')
  exact[n] <- GetExp(n = n, 
                  par = par,
                  pop = pop, 
                  tol = tol, 
                  num = num)
  
}

plot(aprx, type = 'l')
lines(exact, col = 'red')
