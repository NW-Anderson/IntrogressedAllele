# par <-        c(1000,.5,      .5,    .5,      .5,    .5,      .5,    .5,      .5,    .5,.5)
# names(par) <- c(pop, theta.a, phi.a, theta.i, phi.i, theta.n, phi.n, theta.b, phi.b, h, d.m)
DoMatrixApprox <- function(par, n){
  a <- theta.a * theta.i * (1 - phi.a) * (1 - phi.i)
  b <- theta.a * theta.n * (1 - phi.a) * (1 - phi.n) * h
  c <- theta.a * theta.n * (1 - phi.a) * (1 - phi.n) * (1 - h)
  d <- 1- d.m
  m <- theta.b *(1 - phi.b)
  if(c == 1){
    z <- pop + d * m * (n - 1)
    
    sm <- 0
    for(k in 0:(n-1)){
      sm <- sm + (a^(n-1-k) * (pop + d * m * (k - 1)))
    }
    i <- a^n + b * sm
  }else{
    z <- c^n * (pop + d * m /(c-1)) - ( d * m /(c-1))
    
    sm <- 0
    for(k in 0:(n-1)){
      sm <- sm + (a^(n-1-k) * ( c^k * (pop + d * m /(c-1)) - ( d * m /(c-1))))
    }
    i <- a^n + b * sm
  }
  prop <- i / (i + z + m)
  res <- c(prop, i, z, m)
  names(res) <- c('frequency','introgressed', 'native', 'migrant')
  return(prop)
}