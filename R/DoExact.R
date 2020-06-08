
# par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)
par <-   c(.95,   .1,    1,       20,      .1,    1,       .5,  .5,    12,      .05)
pop = 1000 
tol = 10^-5 
num = 100
n=5

aprx <- vector(length = 200)
for(n in 1:200){
  temp <- GetExp(n = n, 
                 par = par,
                 pop = pop, 
                 tol = 10^-5, 
                 num = 100)
}