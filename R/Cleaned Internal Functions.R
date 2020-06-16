# #          1      2      3        4        5      6        7    8      9        10
# # par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)
# par <-   c(.95,   .1,    1,       20,      .1,    1,       .5,  .5,    12,      .05)
# pop = 1000 
# tol = 10^-3 
# num = 100
# n <- 1

Psi <- function(s, n, pop, par){
  return(varPhi_I(s, n, par) * varPhi_N(s,n,par)^pop * psi_M(s, n, par))
}

dPsi <- function(s, n, pop, par){
  term1 <- dvarPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop) * psi_M(s, n, par)
  term2 <- varPhi_I(s, n, par)* pop*(varPhi_N(s, n, par)^(pop-1))* dvarPhi_N(s, n, par)*psi_M(s, n, par)
  term3 <- varPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop)*dpsi_M(s, n, par)
  return(term1 + term2 + term3)
} 

varPhi_I<- function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  if(n==1){
    return(phi_A + (1-phi_A)*exp(theta_A*(phi_I + (1-phi_I)*exp(theta_I*(s-1))-1))) 
  }else{
    return(phi_A + (1-phi_A)*exp(theta_A*(phi_I + (1-phi_I)*exp(theta_I*(varPhi_I(s, n-1, par)-1))-1)))
  }
}

dvarPhi_I<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  if(n==1){
    return(theta_I*theta_A*(1-phi_A)*(1-phi_I)*exp(theta_A*(phi_I+(1-phi_I)*exp(theta_I*(s-1))-1)+
                                                     theta_I*(s-1)))
  }else{
    return(theta_I*theta_A*(1-phi_A)*(1-phi_I)*dvarPhi_I(s,n-1, par)*
             exp(theta_A*(phi_I+(1-phi_I)*exp(theta_I*(varPhi_I(s, n-1, par)-1))-1)+
                   theta_I*(varPhi_I(s, n-1, par)-1)))
  }
}

varPhi_N<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  if(n==1){
    return(phi_A+(1-phi_A)*exp(theta_A*(phi_N+(1-phi_N)*exp(theta_N*(s-1))-1)))
  }else{
    return(phi_A + (1-phi_A)*exp(theta_A*(h*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_I(s, n-1, par)-1)))
                                           +(1-h)*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_N(s, n-1, par)-1)))-1)))
  }
}

dvarPhi_N<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  if(n==1){
    return((1-phi_A)*theta_A*((1-phi_N)*theta_N*h*exp(theta_N*(s-1)))*exp(theta_A*(phi_N+(1-phi_N)*exp(theta_N*(s-1))-1)))
  }else{
    return((1-phi_A)*theta_A*(h*(1-phi_N)*theta_N*dvarPhi_I(s, n-1, par)*exp(theta_N*(varPhi_I(s, n-1, par)-1))+
                                (1-h)*(1-phi_N)*theta_N*dvarPhi_N(s, n-1, par)*exp(theta_N*(varPhi_N(s, n-1, par)-1)))
           *exp(theta_A*(h*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_I(s,n-1,par)-1)))+
                           (1-h)*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_N(s, n-1, par)-1)))-1)))
  }
}

psi_M<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  
  if(n==1){
    return(phi_B+(1-phi_B)*exp(theta_B*(s-1)))
  }
  if(n==2){
    return((phi_B+(1-phi_B)*exp(theta_B*(s-1)))*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*s-1))))
  }
  if(n>=3){
    return((phi_B+(1-phi_B)*exp(theta_B*(s-1)))*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*s-1)))*psi_prod_helper(s, n, par))
  }
}

psi_prod_helper<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  
  i <- 0
  prod = 1
  while(i<=(n-2)){
    i=i+1
    prod = prod*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*varPhi_N(s,i,par)-1)))
  }
  return(prod)
}

dpsi_M<-function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  
  if(n <= 2){
    return(0)
  }
  if(n >= 3){
    return((phi_B+(1-phi_B)*exp(theta_B*(s-1)))*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*s-1)))*psi_sum(s, n, par))
  }
}

psi_sum <- function(s, n, par){
  phi_A <- par[1]
  phi_I <- par[2]
  theta_I <- par[3]
  theta_A <- par[4]
  phi_N <- par[5]
  theta_N <- par[6]
  d_M <- par[7]
  phi_B <- par[8]
  theta_B <- par[9]
  h <- par[10]
  
  i <- 0
  sum = 0
  while(i<=(n-2)){
    i = i+1
    num <- ((1-phi_B)*theta_B*(1-d_M)*dvarPhi_N(s, i, par)*exp(theta_B*(d_M+(1-d_M)*varPhi_N(s, i, par)-1)))*psi_prod_helper(s, n, par)
    den <- (phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*varPhi_N(s,i,par)-1)))
    sum <- sum + (num / den)
  }
  return(sum)
}
