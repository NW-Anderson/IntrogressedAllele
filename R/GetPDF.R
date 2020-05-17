GetPDF <- function(s, n, pop, par){#n = generation number
  term1 <- dvarPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop) * psi_M(s, n, par)
  term2 <- varPhi_I(s, n, par)* pop*(varPhi_N(s, n, par)^(pop-1))* dvarPhi_N(s, n, par)*psi_M(s, n, par)
  term3 <- varPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop)*dpsi_M(s, n, par)
  return(term1 + term2 + term3)
} 

dvarPhi_I<-function(s, n, par){##par[1] = d_I, par[2] = phi_I, par[3] = theta_I 
  if(n==0){
    return((1-par[1])*(par[2]+(par[3]+1)*(1-par[2])*exp(par[3]*(s-1))))
  }else{
    return((1-par[1])*dvarPhi_I(s,n-1, par)*(par[2]+par[3]*varPhi_I(s, n-1, par)+1)*(1-par[2])*exp(par[3]*varPhi_I(s, n-1, par)-1))
  }
}
  
varPhi_I<- function(s, n, par){
  if(n==0){
    return(par[1] + (1-par[1])*(par[2]+(1-par[2]*exp(par[3]*(s-1))))) 
  }else{
    return(par[1] + (1-par[1])*(varPhi_I(s, n-1, par))*(par[2] + (1-par[2])*exp(par[3]*(varPhi_I(s, n-1, par)-1))))
  }
}

dvarPhi_N<-function(s, n, par){##par[4] = d_N, par[5] = phi_N, par[6] = theta_N, par[10] = h
  if(n==0){
    return((1-par[4])*((par[5]+(1-par[5])*exp(-par[6])*(1+(exp(par[6]*par[10]*s)-1)*(exp(par[6]*(1-par[10])*s)-1)))+
                         s*((1-par[5])*exp(-par[6])*(par[6]*par[10]*exp(par[6]*par[10]*s)*(exp(par[6]*(1-par[10])*s)-1)
                                             +(exp(par[6]*par[10]*s)-1)*par[6]*(1-par[10])*exp(par[6]*(1-par[10])*s)))))
  }else{
    return((1-par[4])*dvarPhi_N(s, n-1, par)*(par[5]+(1-par[5])*exp(-par[6])*(1+(exp(par[6]*par[10]*varPhi_I(s, n-1, par))-1)*(exp(par[6]*(1-par[10])*varPhi_N(s, n-1, par)-1)))+
                         varPhi_N(s, n-1, par)*((1-par[5])*exp(-par[6])*(par[6]*par[10]*dvarPhi_I(s, n-1, par)*exp(par[6]*par[10]*varPhi_N(s, n-1, par))*(exp(par[6]*(1-par[10])*varPhi_N(s, n-1,par))-1))
                                             +(exp(par[6]*par[10]*varPhi_N(s, n-1,par))-1)*par[6]*(1-par[10])*dvarPhi_N(s, n-1,par)*exp(par[6]*(1-par[10])*varPhi_N(s, n-1,par)))))
  }
}

varPhi_N<-function(s, n, par){
  if(n==0){
    return(par[4] +((1-par[4])*s*(par[5]+(1-par[5])*exp(-par[6])*(1+(exp(par[6]*par[10]*s)-1)*(exp(par[6]*(1-par[10])*s)-1)))))
  }else{
    return(par[4] +((1-par[4])*varPhi_N(s, n-1, par)*(par[5]+(1-par[5])*exp(-par[6])*(1+(exp(par[6]*par[10]*varPhi_I(s, n-1, par))-1)*(exp(par[6]*(1-par[10])*varPhi_N(s, n-1, par))-1)))))
  }
}

psi_M<-function(s, n, par){###par[7] = d_M, par[8] = phi_B, par[9] = theta_B
  return((par[8] + (1-par[8])*exp(par[9]*(s-1)))*(par[8] + (1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)))*psi_prod(s, n, par))
}

dpsi_M<-function(s, n, par){
  return((par[8]+(1-par[8])*exp(par[9]*(s-1))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)*psi_sum(s, n, par)))))
}

psi_sum <- function(s, n, par){
  i = 1
  sum = 0
  while(i<=(n-2)){
    i = i+1
    sum = sum + ((1-par[8])*par[9]*(1-par[7])*dvarPhi_N(s, i, par)*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s, i, par)-1)))*psi_prod_helper(s, n, par)
    sum = sum - ((1-par[8])*par[9]*(1-par[7])*dvarPhi_N(s, i, par)*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s, i, par)-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s,i,par)-1)))
  }
  return(sum)
}

psi_prod<-function(s, n, par){
  if(n == 0){
    return(1)
  }else if(n==1){
    return(par[8]+(1-par[8])*exp(par[9]*(s-1)))
  }else if(n==2){
    return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1))))
  }else{
    return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)))*psi_prod_helper(s,n,par))
  }
}

psi_prod_helper<-function(s, n, par){
  i = 1
  prod = 1
  while(i<=(n-2)){
    i=i+1
    prod = prod*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s,i,par)-1)))
  }
  return(prod)
}
