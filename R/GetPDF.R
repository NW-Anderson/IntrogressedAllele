GetPDF <- function(n, pop, par){#n = generation number
  term1 <- dvarPhi_I(n, par)*(varPhi_N(n, par)^pop) * psi_M(n, par)
  term2 <- varPhi_I(n, par)* pop*(varPhi_N(n, par)^(pop-1))* dvarPhi_N(n, par)*psi_M(n, par)
  term3 <- varPhi_I(n, par)*(varPhi_N(n, par)^pop)*dpsi_M(n, par)
  return(term1 + term2 + term3)
} 

dvarPhi_I<-function(n, par){##par[1] = d_I, par[2] = phi_I, par[3] = theta_I 
  if(n==0){
    return((1-par[1])*(par[2]+(par[3]+1)*(1-par[2])*exp(par[3]*(0))))
  }else{
    return((1-par[1])*dvarPhi_I(n-1, par)*(par[2]+par[3]*varPhi_I(n-1, par)+1)*(1-par[2])*exp(par[3]*varPhi_I(n-1, par)-1))
  }
}
  
varPhi_I<- function(n, par){
  if(n==0){
    return(par[1] + (1-par[1])*(par[2]+(1-par[2]*exp(par[3]*0)))) 
  }else{
    return(par[1] + (1-par[1])*(varPhi_I(n-1, par))*(par[2] + (1-par[2])*exp(par[3]*(varPhi_I(n-1, par)-1))))
  }
}

dvarPhi_N<-function(n, par){##par[4] = d_N, par[5] = phi_N, par[6] = theta_N, par[10] = h
  if(n==0){
    return((1-par[4])*(par[5]+(1-par[5])*exp(par[6])*(1+(exp(par[6]*par[10]*varPhi_I(n-1, par))-1)*(exp(par[6]*(1-par[10])*varPhi_N(n-1, par)-1)))+
                         varPhi_N(n-1, par)*((1-par[5])*exp(-par[6])*(par[6]*par[10]*dvarPhi_I(n-1, par)*exp(par[6]*par[10]*varPhi_N(n-1, par))*(exp(par[6]*(1-par[10])*varPhi_N(n-1,par))-1))
                                             +(exp(par[6]*par[10]*varPhi_N(n-1,par))-1)*par[6]*(1-par[10])*dvarPhi_N(n-1,par)*exp(par[6]*(1-par[10])*varPhi_N(n-1,par)))))
  }
}
