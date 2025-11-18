nnd=function(acf,sampleT,lag.max,arord=Inf,b=NULL){
  # check if NND
  # Inf is used as the default so that bandtap by default moves towards IID
  nser=dim(acf)[2]
  if(length(arord)==1){rep(arord,nser)}
  
  acfstar=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
    gamma=c(1,acf[,i])
    Gamma=toeplitz(gamma)
    ei=eigen(Gamma)
    if(any(ei$values<=0)){ # NND
      if(arord[i]<sampleT^{1/3}){ # modeled by a "low" order AR, move towards that
        # using n^{1/3} (had log(n) and log(n)*2 previously), trying to get a balance for larger n as log(n) is too small
        rr=b[1:lag.max,i] # invert pacf
        R=toeplitz(c(1,rr))
        alpha=min(eigen(R)$values)
        beta=abs(min(ei$values))*(1+1/sampleT)
        cv=beta/(alpha+beta)
        acf[,i]=cv*rr+(1-cv)*acf[,i]
      }
      else{ # move towards IID 
        R=diag(rep(1,(lag.max+1)))
        alpha=min(eigen(R)$values)
        beta=abs(min(ei$values))*(1+1/sampleT)
        cv=beta/(alpha+beta)
        acf[,i]=(1-cv)*acf[,i]
      }
      Rfinal=toeplitz(c(1,acf[,i]))
    }
    return(acf[,i])
  }) # lag.max x nser
  return(acfstar)  
}