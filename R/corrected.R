######
#' @title corrected
#'
#' @description Calculated the corrected (partial) autocorrelation function
#' @param x a univariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample PACF. Defaults to \eqn{10log_{10}(N)} where N is the number of non-missing observations.
#' @param type 'character'. "correlation", "covariance" or "partial".
#' @param na.action function to be called to handle missing values. Default is na.fail, na.pass can be used.
#' @param demean  'logical'. If 'TRUE' (the default), mean(x) is removed prior to estimation.
#' @param lh vector of length 1 (value used for all lags), or length lag.max. Default uses formula in the description.
#' @param lambda controls the degree of shrinkage towards the target.  Could be a single value (repeated for all 1:lag.max), a single vector of length lag.max (repeated for all nser), or a lag.max x nser matrix. Default is data driven choice.
#' @param target the unbiased (partial) autocorrelation function from a (model) assumption.  Could be a single value (repeated for all 1:lag.max), a single vector of length lag.max (repeated for all nser), or a lag.max x nser matrix. Default is data driven choice.
#' @param ... additional arguments passed to plotting.
#'
#' @return An object of type acf with the following elements:
#' \describe{
#' \item{\code{acf}}{A max.lag x nseries x nseries array containing the estimated penalized acf/pacf.}
#' \item{\code{type}}{Character vector returning the type argument requested.}
#' \item{\code{n.used}}{Numeric of the number of points used for estimation after na.action has been applied.}
#' \item{\code{lag}}{A max.lag x nseries x nseries array containing the lags at which the acf/pacf is estimated.}
#' \item{\code{series}}{The name of the time series.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' \item{\code{lh}}{Matrix \code{lag.max} x \code{nseries} of the lh used in the estimate when \code{penalized==TRUE}.}
#' \item{\code{lambda}}{Matrix \code{lag.max} x \code{nseries} of the lambda used in the estimate when \code{penalized==TRUE}.}
#' \item{\code{target}}{Matrix \code{lag.max} x \code{nseries} of the target used in the estimate when \code{penalized==TRUE}.}
#' }
#'
#'
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' 
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Example for penalized acf estimation
#' corrected(data)
#' }
#' @keywords internal
#' @importFrom stats na.fail
#' @importFrom stats as.ts
#' @importFrom stats toeplitz
#####

corrected = function(x, lag.max = NULL, type = c("correlation", "covariance", 
          "partial"), na.action = na.fail, demean = TRUE, 
          lh = NULL, lambda = NULL, target = NULL,...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  
  cov=FALSE
  if(type=="covariance"){
    cov=TRUE
    type=="correlation"
    # set the flag and proceed as if you wanted to calculate correlation
  }

  if(is.null(lag.max)){lag.max <- floor(10 * (log10(sampleT) - log10(nser)))}

  h=1:lag.max
  
  # if pacf then the next line does pacf anyway
  acf=stats::acf(x,lag.max=max(lag.max,ceiling(log(sampleT))+1),type=type,plot=FALSE,na.action=na.action,demean=demean,...)
  # taking a (potentially) larger lag.max here as we need the larger acf values for the bias correction terms.
  if(type=="partial"){
    tmpacf=acf$acf[1:lag.max,,,drop=FALSE]
  }# need to take a copy so dimensions are the same as we want the output, original is needed for k lags for bias
  else{
    tmpacf=acf$acf[2:(lag.max+1),,,drop=FALSE]
    bacf=acf$acf[2:(ceiling(log(sampleT))+1),,,drop=FALSE] # needed for the bias correction calculation, upper bound here matches j below
  } # need to take a copy so we can remove the lag0 when we do acf and it matches with pacf
  
  j=1:ceiling(log(sampleT))
  if(type=="partial"){j=j-1}
  
  if(length(lh)==1){lh=matrix(rep(lh,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(lh)==lag.max){ # same lh vector for all series
    lh = matrix(lh, nrow = lag.max,ncol=nser)
  }
  else if(length(lh)==nser){ # one value per series, repeat for lags
    lh=matrix(rep(lh,each=lag.max),nrow=lag.max)
  }
  else if(is.null(lh)){ # none supplied so use default
    lh=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
      el = sqrt(log(sampleT)/(sampleT))
      el=el*sqrt(1-tmpacf[,i,i]^2)
      return(el)
    }) # lh is a matrix lag.max x nser
  }
  else if(is.null(dim(lh))){stop("lh must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")}
  else if(any(dim(lh)!=c(lag.max,nser))){ # not something we expect so stop
    stop("lh must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
  }
  if(any(!is.numeric(lh))){stop("lh must be numeric")}
  if(lag.max==1){lh=matrix(lh,nrow=1)}
  else if(nser==1){lh=matrix(lh,ncol=1)}
  
  nserIndexM=matrix(1:nser,ncol=1)
  # bias correction calculation
  if(type=="partial"){
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      b=tmpacf[,i,i]+((h+1)*tmpacf[,i,i]+(tmpacf[,i,i]+1+h%%2==0))/sampleT
              b=sign(b)*pmin(abs(b),0.99)
      return(b)
    }) # returns lag.max x nser matrix
  }
  else{
    # estimate AR order, if "low" then move toward an inverted PACF rather than bias corrected
    arord=apply(nserIndexM, MARGIN = 1, FUN = function(i) {
      return(ar.penyw(x[,i])$order)
    })
    
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      if(arord[i]<(sampleT)^{1/3}){
        # if approximated well by a "low" order AR process
        # using n^{1/3} (had log(n) and log(n)*2 previously), trying to get a balance for larger n as log(n) is too small
        return(acf(x[,i],plot=FALSE,lag.max=lag.max,estimate="invertpacf")$acf[-1,,])
      }
      b=tmpacf[,i,i]+(h*tmpacf[,i,i]+(1-tmpacf[,i,i])*(1+2*sum((1-j/sampleT)*bacf[j,i,i])))/sampleT
      # use bacf which is longer as we need more terms for the truncated infinite sum, if we just use lag.max it may not be enough terms
      b = sign(b) * pmin(abs(b), 0.99)
      return(b)
    }) # returns lag.max x nser matrix
  }
  # above doesn't necessarily return a matrix so make one if needed
  if(lag.max==1){b=matrix(b,nrow=1)}
  else if(nser==1){b=matrix(b,ncol=1)}
  
  # bias correct if larger than lh, otherwise shrink
  
  if(!(is.numeric(lambda) || is.null(lambda))){stop("lambda must be numeric")}
  if(any(!lambda >= 0)){stop("lambda must be positive")}
  if(length(lambda)==1){lambda=matrix(rep(lambda,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(lambda)==lag.max){ # same lambda vector for all series
    lambda = matrix(lambda, nrow = lag.max,ncol=nser)
  }
  else if(length(lambda)==nser){ # one value per series, repeat for lags
    lambda=matrix(rep(lambda,each=lag.max),nrow=lag.max)
  }
  else if(is.null(lambda)){ # none supplied so use default
    lambda = apply(nserIndexM,MARGIN=1,FUN=function(i){
      ind = (abs(tmpacf[,i,i])>lh[,i])
      lambda = (!ind)*h*10*log10(sampleT) *lh[,i]*(lh[,i]-abs(tmpacf[,i,i]))/abs(tmpacf[,i,i])^{3}+ # shrink more aggressively for larger lags
        (ind)*(abs(tmpacf[,i,i])-lh[,i])*(1-lh[,i])/(1-abs(tmpacf[,i,i]))^2*10*log10(sampleT) # movement towards target doesn't depend on lag
      return(lambda)
    }) # lambda is a matrix lag.max x nser
  }
  else if(ncol(lambda)!=nser | nrow(lambda)!=lag.max){
    stop("lambda must be a lag.max x nser matrix")
  }
  else{ # not something we expect so stop
    stop("lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
  }
  if(lag.max==1){lambda=matrix(lambda,nrow=1)}
  else if(nser==1){lambda=matrix(lambda,ncol=1)}
  
  
  if(!(is.numeric(target) || is.null(target))){stop("target must be numeric")}
  if(length(target)==1){target=matrix(rep(target,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(target)==lag.max){ # same target vector for all series
    target = matrix(target, nrow = lag.max,ncol=nser)
  }
  else if(length(target)==nser){ # one value per series, repeat for lags
    target = matrix(rep(target,each=lag.max),nrow=lag.max)
  }
  else if(is.null(target)){ # none supplied so use default
    target = apply(nserIndexM,MARGIN=1,FUN=function(i){
      ind = (abs(tmpacf[,i,i])>lh[,i])
      target = b[,i]*ind 
      return(target)
    }) # target is a matrix lag.max x nser
  }
  else if(ncol(target)!=nser | nrow(target)!=lag.max){
    stop("target must be a lag.max x nser matrix")
  }
  else{ # not something we expect so stop
    stop("target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
  }
  if(any(abs(target) > 1)){stop("target must be between 1 and -1")}

  if(lag.max==1){target=matrix(target,nrow=1)}
  else if(nser==1){target=matrix(target,ncol=1)}
  
  weights=lambda/(1+lambda) # doesn't drop dimensions

  acfstar=apply(nserIndexM,MARGIN=1,FUN=function(i){
    acfstar=(1-weights[,i])*tmpacf[,i,i] + weights[,i]*target[,i]
            acfstar=acfstar[1:lag.max,drop=FALSE] # needed if lag.max=1
    return(acfstar)
  }) # lag.max x nser

  if(lag.max==1){acfstar=matrix(acfstar,nrow=1)}
  else if(nser==1){acfstar=matrix(acfstar,ncol=1)}
  
  # check if NND
  if(type!="partial"){
    acfstar=nnd(acfstar,nser,lag.max,arord,b) # check if non-negative definite and if not make it so
  }

  if(type=="partial"){
    acf$acf=acf$acf[1:lag.max,,,drop=FALSE]
    acf$lag=acf$lag[1:lag.max,,,drop=FALSE]
  }
  else{
    acf$acf=acf$acf[1:(lag.max+1),,,drop=FALSE]
    acf$lag=acf$lag[1:(lag.max+1),,,drop=FALSE]
  }
  
  for(i in 1:nser){
    if(type=="partial"){
      acf$acf[,i,i]=acfstar[,i]
    }
    else{
      acf$acf[-1,i,i]=acfstar[,i]
    }
  }
  acf$lh=lh
  acf$lambda=lambda
  acf$target=target
  return(acf)
}

