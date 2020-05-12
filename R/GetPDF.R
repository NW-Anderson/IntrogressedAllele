GetPDF <- function(n, pop, par){#n = generation number
  term1 <- dvarPhi_I(n, par)*(varPhi_N(n, par)^pop) * psi_M(n, par)
  term2 <- varPhi_I(n, par)* pop*(varPhi_N(n, par)^(pop-1))* dvarPhi_N(n, par)*psi_M(n, par)
  term3 <- varPhi_I(n, par)*(varPhi_N(n, par)^pop)*dpsi_M(n, par)
  return(term1 + term2 + term3)
} 

dvarPhi_I<-function(n, par){##par[1] = d_I, par[2] = phi_I, par[3] = theta_I 
  if(n==0){
    return((1-par[1])*(par[2]+(par[3]*n+1)*(1-par[2])*exp(par[3]*(n-1))))
  }
  return((1-par[1])*dvarPhi_I(n-1, par)*(par[2]+par[3]*varPhi_I(n-1, par)+1)*(1-par[2])*exp(par[3]*varPhi_I(n-1, par)-1))
}
  