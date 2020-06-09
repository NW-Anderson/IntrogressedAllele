#          1      2      3        4        5      6        7    8      9        10
# par <- c(phi_A, phi_I, theta_I, theta_A, phi_N, theta_N, d_M, phi_B, theta_B, h)

Psi <- function(s, n, pop, par){
  return(varPhi_I(s, n, par) * varPhi_N(s,n,par)^pop * psi_M(s, n, par))
}
GetPGF <- function(s, n, pop, par){#n = generation number
  term1 <- dvarPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop) * psi_M(s, n, par)
  term2 <- varPhi_I(s, n, par)* pop*(varPhi_N(s, n, par)^(pop-1))* dvarPhi_N(s, n, par)*psi_M(s, n, par)
  term3 <- varPhi_I(s, n, par)*(varPhi_N(s, n, par)^pop)*dpsi_M(s, n, par)
  return(term1 + term2 + term3)
} 

dvarPhi_I<-function(s, n, par){##par[1] = phi_A, par[2] = phi_I, par[3] = theta_I, par[4] = theta_A
  if(n==1){
    return(par[3]*par[4]*(1-par[1])*(1-par[2])*exp(par[4]*(par[2]+(1-par[2])*exp(par[3]*(s-1))-1)+par[3]*(s-1)))
  }else{
    return(par[3]*par[4]*(1-par[1])*(1-par[2])*dvarPhi_I(s,n-1, par)*exp(par[4]*(par[2]+(1-par[2])*exp(par[3]*(varPhi_I(s, n-1, par)-1))-1)+par[3]*(varPhi_I(s, n-1, par)-1)))
  }
}
  
varPhi_I<- function(s, n, par){
  if(n==1){
    return(par[1] + (1-par[1])*exp(par[4]*(par[2] + (1-par[2])*exp(par[3]*(s-1))-1))) 
  }else{
    return(par[1] + (1-par[1])*exp(par[4]*(par[2] + (1-par[2])*exp(par[3]*(varPhi_I(s, n-1, par)-1))-1)))
  }
}

dvarPhi_N<-function(s, n, par){##par[4] = theta_A, par[5] = phi_N, par[6] = theta_N, par[10] = h
  if(n==1){
    return((1-par[1])*par[4]*((1-par[5])*par[6]*exp(par[6]*(s-1)))*exp(par[4]*(par[5]+(1-par[5])*exp(par[6]*(s-1))-1)))
  }else{
    return((1-par[1])*par[4]*(par[10]*(1-par[5])*par[6]*dvarPhi_I(s, n-1, par)*exp(par[6]*(varPhi_I(s, n-1, par)-1))+(1-par[10])*(1-par[5])*par[6]*dvarPhi_N(s, n-1, par)*exp(par[6]*(varPhi_N(s, n-1, par)-1)))
                         *exp(par[4]*(par[10]*(par[5]+(1-par[5])*exp(par[6]*(varPhi_I(s,n-1,par)-1)))+(1-par[10])*(par[5]+(1-par[5])*exp(par[6]*(varPhi_N(s, n-1, par)-1)))-1)))
  }
}

varPhi_N<-function(s, n, par){
  if(n==1){
    return(par[1]+(1-par[1])*exp(par[4]*(par[5]+(1-par[5])*exp(par[6]*(s-1))-1)))
    }else{
    return(par[1] + (1-par[1])*exp(par[4]*(par[10]*(par[5]+(1-par[5])*exp(par[6]*(varPhi_I(s, n-1, par)-1)))
                                           +(1-par[10])*(par[5]+(1-par[5])*exp(par[6]*(varPhi_N(s, n-1, par)-1)))-1)))
  }
}

# need to add n = 1 and n = 2 cases
psi_M<-function(s, n, par){###par[7] = d_M, par[8] = phi_B, par[9] = theta_B
  return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)))*psi_prod_helper(s, n, par))
}

dpsi_M<-function(s, n, par){
  return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)))*psi_sum(s, n, par))
}
# this one starts at i = 2 as weell
# psi_sum <- function(s, n, par){
#   ### i = 1
#   i <- 0
#   sum = 0
#   while(i<=(n-2)){
#     i = i+1
#     sum = sum + ((1-par[8])*par[9]*(1-par[7])*dvarPhi_N(s, i, par)*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s, i, par)-1)))*psi_prod_helper(s, n, par)
#     # i think this line is set up wrong. the t!=r applies to the product not the sum. If we want to take it away later 
#     #  like this we will need to divide to undo the unwanted product, not subtract
#     sum = sum - ((1-par[8])*par[9]*(1-par[7])*dvarPhi_N(s, i, par)*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s, i, par)-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s,i,par)-1)))
#   }
#   return(sum)
# }

psi_sum <- function(s, n, par){
  ### i = 1
  i <- 0
  sum = 0
  while(i<=(n-2)){
    i = i+1
    num <- ((1-par[8])*par[9]*(1-par[7])*dvarPhi_N(s, i, par)*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s, i, par)-1)))*psi_prod_helper(s, n, par)
    den <- (par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s,i,par)-1)))
    sum <- sum + (num / den)
  }
  return(sum)
}
# this function is not called anywhere else is it?
psi_prod<-function(s, n, par){
  if(n == 0){
    return(s)
  }else if(n==1){
    return(par[8]+(1-par[8])*exp(par[9]*(s-1)))
  }else if(n==2){
    return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1))))
  }else{
    return((par[8]+(1-par[8])*exp(par[9]*(s-1)))*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*s-1)))*psi_prod_helper(s,n,par))
  }
}
# I believe the way this is written the first iteration of the while loop is i = 2. 
#  moving the i + 1 to after the prod line or initializing with i = 0 could 
#  fix the issue
psi_prod_helper<-function(s, n, par){
  ### i = 1
  i <- 0
  prod = 1
  while(i<=(n-2)){
    i=i+1
    prod = prod*(par[8]+(1-par[8])*exp(par[9]*(par[7]+(1-par[7])*varPhi_N(s,i,par)-1)))
  }
  return(prod)
}
