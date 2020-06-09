GetExp<-function(n, par, pop, tol, num){
  # i think we will need to add a special case for 0, not just enter 0 for n
  # coeff<- 1/(1-GetPGF(0,n,pop,par))
  coeff <- 1/(1-Psi(0,n,pop,par))
  APP = 0
  i = 1
  a = c()
  TOL = c()
  h = c()
  FA = c()
  FC = c()
  FB = c()
  S = c()
  L = c()
  a[i] = 0
  TOL[i] = 10*tol
  h[i] = .5
  FA[i] = GetPGF(0,n,pop,par)
  FC[i] = GetPGF(.5, n, pop, par)
  FB[i] = GetPGF(1, n, pop, par)
  S[i] = h[i]*(FA[i] + 4*FC[i]+FB[i])/3
  L[i] = 1
  v = c()
  while(i>0){
    FD = GetPGF(a[i] + (h[i])/2, n, pop, par)
    FE = GetPGF(a[i] + (3*h[i])/2, n, pop, par)
    S1 = h[i]*(FA[i] + 4*FD + FC[i])/6
    S2 = h[i]*(FC[i] + 4*FE + FB[i])/6
    v[1] = a[i]
    v[2] = FA[i]
    v[3] = FC[i]
    v[4] = FB[i]
    v[5] = h[i]
    v[6] = TOL[i]
    v[7] = S[i]
    v[8] = L[i]
    i = i-1
    if(abs(S1+S2-v[7])<v[6]){
      APP = APP + (S1+S2)
    }else{
      if(v[8] >= num){
        print("LEVEL EXCEEDED")
        break
      }else{
        i = i+1
        a[i] = v[1]+v[5]
        FA[i] = v[3]
        FC[i] = FE
        FB[i] = v[4]
        h[i] = v[5]/2
        TOL[i] = v[6]/2
        S[i] = S2
        L[i] = v[8] + 1
        
        i = i+1
        a[i] = v[1]
        FA[i] = v[2]
        FC[i] = FD
        FB[i] = v[3]
        h[i] = h[i-1]
        TOL[i] = TOL[i-1]
        S[i] = S1
        L[i] = L[i-1]
      }
    }
  }
  return(coeff*APP)
}
