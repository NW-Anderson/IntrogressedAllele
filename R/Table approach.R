###############################
#     params for testing      #
###############################
# ngen <- 40
# par <-   c(.69,   .5,    1,      6,       .5,     1,      .57,  0,    100,    .05)
# pop = 613
# tol = 10^-3 
# num = 100
###############################
#          Functions          #
###############################

# equation 13 in manuscript
varPhi_I<- function(s, n, par, varPhi_I.array){
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
    return(phi_A + (1-phi_A)*exp(theta_A*(phi_I + (1-phi_I)*
                                            exp(theta_I*(varPhi_I.array[n-1,which(colnames(varPhi_I.array) == paste('s =', s))]-1))-1)))
  }
}

# equation 14
varPhi_N<-function(s, n, par, varPhi_I.array, varPhi_N.array){
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
    return(phi_A + (1-phi_A)*exp(theta_A*(h*(phi_N+(1-phi_N)*
                                               exp(theta_N*(varPhi_I.array[n-1,which(colnames(varPhi_I.array)== paste('s =', s))]-1)))+
                                            (1-h)*(phi_N+(1-phi_N)*exp(theta_N*
                                                                         (varPhi_N.array[n-1,which(colnames(varPhi_N.array)== paste('s =', s))]
                                                                          -1)))-1)))
  }
}

# equation 19
psi_M<-function(s, n, par, varPhi_N.array){
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
    return((phi_B+(1-phi_B)*exp(theta_B*(s-1)))*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*s-1)))*
             psi_prod_helper(s, n, par, varPhi_N.array))
  }
}

# product inside of equation 19 and 24
psi_prod_helper<-function(s, n, par, varPhi_N.array){
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
    prod = prod*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*varPhi_N.array[i,which(colnames(varPhi_N.array)== paste('s =', s))]-1)))
  }
  return(prod)
}

# equation 20
Psi <- function(s, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array){
  return(varPhi_I.array[n,which(colnames(varPhi_I.array)== paste('s =', s))] * 
           (varPhi_N.array[n,which(colnames(varPhi_N.array)== paste('s =', s))])^pop * 
           Psi_M.array[n,which(colnames(Psi_M.array)== paste('s =',s))])
}

# equation 22
dvarPhi_I<-function(s, n, par, varPhi_I.array, dvarPhi_I.array){
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
    return(theta_I*theta_A*(1-phi_A)*(1-phi_I)*exp(theta_A*(phi_I+(1-phi_I)*exp(theta_I*(s-1))-1)+theta_I*(s-1)))
  }else{
    return(theta_I*theta_A*(1-phi_A)*(1-phi_I)*dvarPhi_I.array[n-1,which(colnames(dvarPhi_I.array)== paste('s =', s))] *
             exp(theta_A*(phi_I+(1-phi_I)*exp(theta_I*(varPhi_I.array[n-1,which(colnames(varPhi_I.array)== paste('s =', s))]-1))-1)+
                   theta_I*(varPhi_I.array[n-1, which(colnames(varPhi_I.array)== paste('s =', s))]-1)))
  }
}

# equation 23
dvarPhi_N<-function(s, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array){
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
    # terms are in different order but are equal as far as I can tell
    return((1-phi_A)*theta_A*(1-phi_N)*theta_N*h*exp(theta_N*(s-1))*exp(theta_A*(phi_N+(1-phi_N)*exp(theta_N*(s-1))-1)))
  }else{
    return((1-phi_A)*theta_A*(h*(1-phi_N)*theta_N*dvarPhi_I.array[n-1,which(colnames(dvarPhi_I.array) == paste('s =', s))]*
                                exp(theta_N*(varPhi_I.array[n-1,which(colnames(varPhi_I.array)== paste('s =', s))]-1))+
                                (1-h)*(1-phi_N)*theta_N*dvarPhi_N.array[n-1, which(colnames(dvarPhi_N.array) == paste('s =', s))] *
                                exp(theta_N*(varPhi_N.array[n-1,which(colnames(varPhi_N.array)== paste('s =', s))]-1)))*
             exp(theta_A*(h*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_I.array[n-1,which(colnames(varPhi_I.array) == paste('s =', s))]-1)))+
                            (1-h)*(phi_N+(1-phi_N)*exp(theta_N*(varPhi_N.array[n-1, which(colnames(varPhi_N.array)== paste('s =', s))]-
                                                                  1)))-1)))
  }
}

# equation 24
dpsi_M<-function(s, n, par, varPhi_N.array, dvarPhi_N.array){
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
    return((phi_B+(1-phi_B)*exp(theta_B*(s-1)))*(phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*s-1)))*
             psi_sum(s, n, par, varPhi_N.array, dvarPhi_N.array))
  }
}

# sum inside of equation 24
psi_sum <- function(s, n, par, varPhi_N.array, dvarPhi_N.array){
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
    num <- ((1-phi_B)*theta_B*(1-d_M)*dvarPhi_N.array[i,which(colnames(dvarPhi_N.array) == paste('s =', s))]*
              exp(theta_B*(d_M+(1-d_M)*varPhi_N.array[i, which(colnames(varPhi_N.array) == paste('s =', s))]-1)))*
      psi_prod_helper(s, n, par, varPhi_N.array)
    den <- (phi_B+(1-phi_B)*exp(theta_B*(d_M+(1-d_M)*varPhi_N.array[i,which(colnames(varPhi_N.array) == paste('s =', s))]-1)))
    sum <- sum + (num / den)
  }
  return(sum)
}

# equation 21
dPsi <- function(s, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array){
  term1 <- dvarPhi_I.array[n, which(colnames( dvarPhi_I.array) == paste('s =', s))]*
    (varPhi_N.array[n, which(colnames(varPhi_N.array) == paste('s =', s))]^pop) * 
    Psi_M.array[n, which(colnames(Psi_M.array) == paste('s =', s))]
  
  term2 <- varPhi_I.array[n,which(colnames(varPhi_I.array) == paste('s =', s))]* pop*
    (varPhi_N.array[n,which(colnames(varPhi_N.array) == paste('s =', s))]^(pop-1))* 
    dvarPhi_N.array[n,which(colnames(dvarPhi_N.array) == paste('s =', s))]*
    Psi_M.array[n,which(colnames(Psi_M.array) == paste('s =', s))]
  
  term3 <- varPhi_I.array[n,which(colnames(varPhi_I.array)== paste('s =', s))]*
    (varPhi_N.array[n, which(colnames(varPhi_N.array) == paste('s =', s))]^pop)* 
    dPsi_M.array[n, which(colnames(dPsi_M.array)== paste('s =', s))]
  return(term1 + term2 + term3)
} 
################################
#          Table Maker         #
################################
SortTable <- function(tbl){
  svals <- sort(as.numeric(substr(colnames(tbl),5,nchar(colnames(tbl)))))
  new.tbl <- tbl[,paste('s =', svals)]
  return(list('tbl' = new.tbl,
              'svals' = svals))
}

MakeSeedTable <- function(ngen, par, pop){
  ##### Creating empty Tables #####
  varPhi_I.array <- array(dim = c(ngen, 5))
  colnames(varPhi_I.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(varPhi_I.array) <- paste('gen = ', 1:ngen, sep = '')
  
  varPhi_N.array <- array(dim = c(ngen, 5))
  colnames(varPhi_N.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(varPhi_N.array) <- paste('gen = ', 1:ngen, sep = '')
  
  Psi_M.array <- array(dim = c(ngen, 5))
  colnames(Psi_M.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(Psi_M.array) <- paste('gen = ', 1:ngen, sep = '')
  
  Psi.array <- array(dim = c(ngen, 5))
  colnames(Psi.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(Psi.array) <- paste('gen = ', 1:ngen, sep = '')
  
  dvarPhi_I.array <- array(dim = c(ngen, 5))
  colnames(dvarPhi_I.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(dvarPhi_I.array) <- paste('gen = ', 1:ngen, sep = '')
  
  dvarPhi_N.array <- array(dim = c(ngen, 5))
  colnames(dvarPhi_N.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(dvarPhi_N.array) <- paste('gen = ', 1:ngen, sep = '')
  
  dPsi_M.array <- array(dim = c(ngen, 5))
  colnames(dPsi_M.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(dPsi_M.array) <- paste('gen = ', 1:ngen, sep = '')
  
  dPsi.array <- array(dim = c(ngen, 5))
  colnames(dPsi.array) <- paste('s = ', c('0', '0.25', '0.5', '0.75', '1'), sep = '')
  rownames(dPsi.array) <- paste('gen = ', 1:ngen, sep = '')
  
  ##### Filling in Tables #####
  for (n in 1:ngen) {
    varPhi_I.array[n,1] <- varPhi_I(0, n, par, varPhi_I.array)
    varPhi_I.array[n,2] <- varPhi_I(0.25, n, par, varPhi_I.array)
    varPhi_I.array[n,3] <- varPhi_I(0.5, n, par, varPhi_I.array)
    varPhi_I.array[n,4] <- varPhi_I(0.75, n, par, varPhi_I.array)
    varPhi_I.array[n,5] <- varPhi_I(1, n, par, varPhi_I.array)
    
    varPhi_N.array[n,1] <- varPhi_N(0, n, par, varPhi_I.array, varPhi_N.array)
    varPhi_N.array[n,2] <- varPhi_N(0.25, n, par, varPhi_I.array, varPhi_N.array)
    varPhi_N.array[n,3] <- varPhi_N(0.5, n, par, varPhi_I.array, varPhi_N.array)
    varPhi_N.array[n,4] <- varPhi_N(0.75, n, par, varPhi_I.array, varPhi_N.array)
    varPhi_N.array[n,5] <- varPhi_N(1, n, par, varPhi_I.array, varPhi_N.array)
    
    Psi_M.array[n,1] <- psi_M(0, n, par, varPhi_N.array)
    Psi_M.array[n,2] <- psi_M(0.25, n, par, varPhi_N.array)
    Psi_M.array[n,3] <- psi_M(0.5, n, par, varPhi_N.array)
    Psi_M.array[n,4] <- psi_M(0.75, n, par, varPhi_N.array)
    Psi_M.array[n,5] <- psi_M(1, n, par, varPhi_N.array)
    
    Psi.array[n,1] <- Psi(0, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array)
    Psi.array[n,2] <- Psi(0.25, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array)
    Psi.array[n,3] <- Psi(0.5, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array)
    Psi.array[n,4] <- Psi(0.75, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array)
    Psi.array[n,5] <- Psi(1, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array)
    
    dvarPhi_I.array[n,1] <- dvarPhi_I(0, n, par, varPhi_I.array, dvarPhi_I.array)
    dvarPhi_I.array[n,2] <- dvarPhi_I(0.25, n, par, varPhi_I.array, dvarPhi_I.array)
    dvarPhi_I.array[n,3] <- dvarPhi_I(0.5, n, par, varPhi_I.array, dvarPhi_I.array)
    dvarPhi_I.array[n,4] <- dvarPhi_I(0.75, n, par, varPhi_I.array, dvarPhi_I.array)
    dvarPhi_I.array[n,5] <- dvarPhi_I(1, n, par, varPhi_I.array, dvarPhi_I.array)
    
    dvarPhi_N.array[n,1] <- dvarPhi_N(0, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array)
    dvarPhi_N.array[n,2] <- dvarPhi_N(0.25, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array)
    dvarPhi_N.array[n,3] <- dvarPhi_N(0.5, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array)
    dvarPhi_N.array[n,4] <- dvarPhi_N(0.75, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array)
    dvarPhi_N.array[n,5] <- dvarPhi_N(1, n, par, varPhi_I.array, varPhi_N.array, dvarPhi_I.array, dvarPhi_N.array)
    
    dPsi_M.array[n,1] <- dpsi_M(0, n, par, varPhi_N.array, dvarPhi_N.array)
    dPsi_M.array[n,2] <- dpsi_M(0.25, n, par, varPhi_N.array, dvarPhi_N.array)
    dPsi_M.array[n,3] <- dpsi_M(0.5, n, par, varPhi_N.array, dvarPhi_N.array)
    dPsi_M.array[n,4] <- dpsi_M(0.75, n, par, varPhi_N.array, dvarPhi_N.array)
    dPsi_M.array[n,5] <- dpsi_M(1, n, par, varPhi_N.array, dvarPhi_N.array)
    
    dPsi.array[n,1] <- dPsi(0, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array)
    dPsi.array[n,2] <- dPsi(0.25, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array)
    dPsi.array[n,3] <- dPsi(0.5, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array)
    dPsi.array[n,4] <- dPsi(0.75, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array)
    dPsi.array[n,5] <- dPsi(1, n, pop, par, varPhi_I.array, varPhi_N.array, Psi_M.array, dvarPhi_I.array, dvarPhi_N.array, dPsi_M.array)
  }
  return(list('varPhi_I.array' = varPhi_I.array, 
              'varPhi_N.array' = varPhi_N.array, 
              'Psi_M.array' = Psi_M.array, 
              'Psi.array' = Psi.array,
              'dvarPhi_I.array' = dvarPhi_I.array,
              'dvarPhi_N.array' = dvarPhi_N.array, 
              'dPsi_M.array' = dPsi_M.array, 
              'dPsi.array' = dPsi.array))
}

MakeNewcolumn <- function(ngen, par, pop, s, tables){
  new.col <- array(dim = c(ngen,1))
  rownames(new.col) <-  paste('gen = ', 1:ngen, sep = '')
  colnames(new.col) <- paste('s =', s)
  
  tables$varPhi_I.array <- cbind(tables$varPhi_I.array, new.col)
  tables$varPhi_N.array <- cbind(tables$varPhi_N.array, new.col)
  tables$Psi_M.array <- cbind(tables$Psi_M.array, new.col)
  tables$Psi.array <- cbind(tables$Psi.array, new.col)
  tables$dvarPhi_I.array <- cbind(tables$dvarPhi_I.array, new.col)
  tables$dvarPhi_N.array <- cbind(tables$dvarPhi_N.array, new.col)
  tables$dPsi_M.array <- cbind(tables$dPsi_M.array, new.col)
  tables$dPsi.array <- cbind(tables$dPsi.array, new.col)
  for (n in 1:ngen) {
    tables$varPhi_I.array[n,ncol(tables$varPhi_I.array)] <- varPhi_I(s, n, par, 
                                                                     tables$varPhi_I.array)
    
    tables$varPhi_N.array[n,ncol(tables$varPhi_I.array)] <- varPhi_N(s, n, par, 
                                                                     tables$varPhi_I.array, 
                                                                     tables$varPhi_N.array)
    
    tables$Psi_M.array[n,ncol(tables$Psi_M.array)] <- psi_M(s, n, par, 
                                                            tables$varPhi_N.array)
    
    tables$Psi.array[n,ncol(tables$Psi.array)] <- Psi(s, n, pop, par, 
                                                      tables$varPhi_I.array, 
                                                      tables$varPhi_N.array, 
                                                      tables$Psi_M.array)
    
    tables$dvarPhi_I.array[n,ncol(tables$dvarPhi_I.array)] <- dvarPhi_I(s, n, par, 
                                                                        tables$varPhi_I.array, 
                                                                        tables$dvarPhi_I.array)
    
    tables$dvarPhi_N.array[n,ncol(tables$dvarPhi_N.array)] <- dvarPhi_N(s, n, par, 
                                                                        tables$varPhi_I.array, 
                                                                        tables$varPhi_N.array, 
                                                                        tables$dvarPhi_I.array, 
                                                                        tables$dvarPhi_N.array)
    
    tables$dPsi_M.array[n,ncol(tables$dPsi_M.array)] <- dpsi_M(s, n, par, 
                                                               tables$varPhi_N.array, 
                                                               tables$dvarPhi_N.array)
    
    tables$dPsi.array[n,ncol(tables$dPsi.array)] <- dPsi(s, n, pop, par, 
                                                         tables$varPhi_I.array, 
                                                         tables$varPhi_N.array, 
                                                         tables$Psi_M.array, 
                                                         tables$dvarPhi_I.array, 
                                                         tables$dvarPhi_N.array, 
                                                         tables$dPsi_M.array)
  }
  return(tables)
}
################################
#    Derivative Calculator     #
################################
FindDerivs <- function(tbl){
  secder <- array(dim = c(nrow(tbl),ncol(tbl)))
  rownames(secder) <- rownames(tbl)
  colnames(secder) <- colnames(tbl)
  
  temp <- SortTable(tbl)
  tbl <- temp[[1]]
  svals <- temp[[2]]
  rm(temp)
  for(r in 1:nrow(tbl)){
    for(c in 1:ncol(tbl)){
      if(c == 1){
        secder[r,c] <- threept(tbl[r,c],tbl[r,c+1], tbl[r,c+2],svals[c],svals[c+1],svals[c+2])
      }else if(c == ncol(tbl)){
        secder[r,c] <- threept(tbl[r,c],tbl[r,c-1], tbl[r,c-2],svals[c],svals[c-1],svals[c-2])
      }else{
        secder[r,c] <- threept(tbl[r,c],tbl[r,c-1], tbl[r,c+1],svals[c],svals[c-1],svals[c+1])
      }
    }
  }
  return(secder)
}

threept <- function(fx0,fx1,fx2,x0,x1,x2){
  return(fx0 * ((2*x0 -x1 - x2)/((x0-x1)*(x0-x2))) + 
    fx1 * ((2*x0 -x0 - x2)/((x1-x0)*(x1-x2))) +
    fx2 * ((2*x0 -x0 - x1)/((x2-x0)*(x2-x1))))
}

################################
#      Composite Simpsons      #
################################
CompositeSimpsonsWrapper <- function(Psi.array, deriv.array, secderiv.array, svals){
  difs <- svals[2:length(svals)] - svals[1:(length(svals)-1)]
  difs.index <-array(F,dim = c(length(svals), length(names(table(difs)))))
  rownames(difs.index) <- paste('s =', svals)
  colnames(difs.index) <- names(table(difs))
  for(i in 1:length(difs)){
    difs.index[i,which(as.character(difs[i]) == colnames(difs.index))] <- T
    difs.index[i + 1,which(as.character(difs[i]) == colnames(difs.index))] <- T
  }
  vars <- vector(length = nrow(deriv.array))
  for (n in 1:length(vars)) {
    total <- 0
    for (c in 1:ncol(difs.index)) {
      ssubset <- svals[difs.index[,c]]
      total <- total + CompositeSimpsons(svals = ssubset, 
                                         derivvals = deriv.array[n,], 
                                         secderivvals = secderiv.array[n,])
    }
    vars[n] <- 1/(1-Psi.array[n,1])* total
  }
  return(vars)
}

CompositeSimpsons <- function(svals, 
                  derivvals, 
                  secderivvals){
  a <- min(svals)
  b <- max(svals)
  n <- length(svals) - 1
  h <- (b-a)/n
  xi0 <- InsideIntegral(s = a,
                        der = derivvals[which(paste('s =', a) == names(derivvals))],
                        secder = secderivvals[which(paste('s =', a) == names(secderivvals))])+
    InsideIntegral(s = b,
                   der = derivvals[which(paste('s =', b) == names(derivvals))],
                   secder = secderivvals[which(paste('s =', b) == names(secderivvals))])
  xi1 <- 0
  xi2 <- 0
  for(i in 2:n){
    if(i %% 2 == 0) xi2 <- xi2 + InsideIntegral(s = svals[i],
                                                der = derivvals[which(paste('s =', svals[i]) == names(derivvals))],
                                                secder = secderivvals[which(paste('s =', svals[i]) == names(secderivvals))])
    if(i %% 2 == 1) xi1 <- xi1 + InsideIntegral(s = svals[i],
                                                der = derivvals[which(paste('s =', svals[i]) == names(derivvals))],
                                                secder = secderivvals[which(paste('s =', svals[i]) == names(secderivvals))])
  }
  return(h * (xi0 + 2*xi2 + 4* xi1)/3)
}

InsideIntegral <-  function(s,der,secder){
  if(s == 0) s <- 10^(-25)
  return(-log(s) * (s * secder + der))
}











