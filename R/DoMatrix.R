# par <-        c(1000,.5,      .5,    .5,      .5,    .5,      .5,    .5,      .5,    .5,.5)
# names(par) <- c(pop, theta.a, phi.a, theta.i, phi.i, theta.n, phi.n, theta.b, phi.b, h, d.m)
# n <- 8
DoMatrixApprox <- function(par, n){
  pop <- par[1]
  theta.a <- par[2]
  phi.a <- par[3] 
  theta.i <- par[4] 
  phi.i <- par[5] 
  theta.n <- par[6] 
  phi.n <- par[7] 
  theta.b <- par[8] 
  phi.b <- par[9] 
  h <- par[10] 
  d.m <- par[11]
  
  a <- theta.a * theta.i * (1 - phi.a) * (1 - phi.i)
  b <- theta.a * theta.n * (1 - phi.a) * (1 - phi.n) * h
  c <- theta.a * theta.n * (1 - phi.a) * (1 - phi.n) * (1 - h)
  d <- 1- d.m
  m <- theta.b *(1 - phi.b)
  
  z <- vector(length = n)
  z[1] <- c * pop
  if(n ==1){
    i <- a + b * pop
  }
  if(n >= 2){
    if(c == 1){
      for(k in 2:n){
        z[k] <- pop + d * m * (k - 1)
      }
    }else{
      for(k in 2:n){
        z[k] <- c^k * (pop + d * m /(c-1)) - ( d * m /(c-1))
      }
    }
    
    sm <- (a^(n-1) * pop)
    for(k in 1:(n-1)){
      sm <- sm + (a^(n-1-k) * z[k])
    }
    i <- a^n + b * sm
  }
  prop <- i / (i + z[n] + m)
  res <- c(prop, i, z[n], m)
  names(res) <- c('frequency','introgressed', 'native', 'migrant')
  return(res)
}