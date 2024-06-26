######
#' @title acf
#'
#' @description The Penalized (Partial) Autocorrelation Function (ACF/PACF) Estimation
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample ACF/PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param type type of acf to be computed. Allow values are "correlation" (the default) or "partial".
#' @param plot  'logical'. If 'TRUE' (the default) the penalized/sample ACF/PACF is plotted.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized ACF/PACF is computed; if 'FALSE' the sample ACF/PACF is computed.
#' @param lh sequence of threshold values across h. Could be a single value (repeated for all h), a single vector of length h (repeated for all nser), or a h x nser matrix. Default is data driven.
#' @param estimate character vector of the estimation method for the ACF, options are "direct" (default) or "invertpacf".  "invertpacf" is preferred when the data can be approximated by a low order (<3) ARMA model.
#' @param ... additional arguments for specific methods.
#'
#' @return An object of penalized/sample ACF/PACF estimation with the following elements:
#' \describe{
#' \item{\code{acf}}{A max.lag x nseries x nseries array containing the estimated penalized acf/pacf.}
#' \item{\code{type}}{Character vector returning the type argument requested.}
#' \item{\code{n.used}}{Numeric of the number of points used for estimation after na.action has been applied.}
#' \item{\code{lag}}{A max.lag x nseries x nseries array containing the lags at which the acf/pacf is estimated.}
#' \item{\code{series}}{The name of the time series.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' \item{\code{penalized}}{Logical indicating if the acf/pacf returned is penalized.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Examples for penalized ACF/PACF and sample ACF/PACF
#' acf(data)
#' acf(data, penalized = FALSE)
#' acf(data, type ="partial")
#' acf(data, type ="partial", penalized = FALSE)
#'
#' x1 <- arima.sim(n=100, model=list(ar=0.5))
#' x2 <- arima.sim(n=100, model=list(ar=0.1))
#' x3 <- arima.sim(n=100, model=list(ar=0.9))
#' x <- cbind(x1, x2, x3)
#'
#' acf(x)
#' acf(x, penalized = FALSE)
#' acf(x, type ="partial")
#' acf(x, type ="partial", penalized = FALSE)
#' }
#' @export
#' @importFrom stats as.ts
#' @importFrom stats na.fail
#' @importFrom stats var
#####

acf <-
function (x, lag.max = NULL, type = c("correlation", "covariance", 
    "partial"), plot = TRUE, na.action = na.fail, demean = TRUE, 
    penalized=TRUE,lh=NULL,estimate="direct",...){
    type <- match.arg(type)
    if (type == "partial") {
        m <- match.call()
        m[[1L]] <- quote(pacf)
        m$type <- NULL
        return(eval(m, parent.frame()))
    }
    if(estimate=="direct"){
      if(!penalized){ # not penalised so run usual acf
        acf=stats::acf(x,lag.max,type,plot=F,na.action,demean,...)
        acf$penalized=FALSE
        acf$lh=NULL
      }
      else{ #run penalised estimation
        acf=corrected(x,lag.max,type,na.action,demean,lh,...)
        acf$penalized=TRUE
      }
      acf$estimate="direct"
    }
    else if(estimate=="invertpacf"){
      x <- na.action(as.ts(x))
      x <- as.matrix(x)
      if (!is.numeric(x)) 
        stop("'x' must be numeric")
      sampleT=nrow(x)
      nser=ncol(x)
      if (is.null(lag.max)) 
        lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
      if(is.matrix(lh)){
        if(dim(lh)[1]!=lag.max){
          stop("lh is a matrix with the incorrect dimension, must have nrow(lh)=lag.max.")
        }
        else if(dim(lh)[2]!=nser){
          stop("lh is a matrix with the incorrect dimension, must have ncol(lh)=nser.")
        }
      }
      
      # First get the approximate AR order by AIC
      xpacf=pacf(x,lag.max=lag.max,plot=F,na.action=na.action,penalized=penalized,lh=lh)
      AICpen <- apply(matrix(1:lag.max,ncol=1), MARGIN=1,
                      FUN=function(i){
                      sampleT*log(apply(x,MARGIN=2,FUN=var)* # nser length of variances
                      apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
                        prod(1-xpacf$acf[1:i,j,j]^2) # nser length of products
                      })) + 2*i
      })
      if(nser==1){
        AICpen=matrix(AICpen,nrow=1)
      }
      AICpen <- cbind(sampleT * log(apply(x,MARGIN=2,FUN=var)),AICpen)
      
      xarorder <- apply(AICpen,MARGIN=1,FUN=which.min)-1
      if(any(xarorder>3)){
        warning("Approximated AR order is larger than 3 in atleast one of the series. Returning the penalized direct estimator instead.")
      }
      for(i in 1:nser){ # zero the pacf above the fitted ar lag
        if(xarorder[i]!=lag.max){
          xpacf$acf[(xarorder[i]+1):lag.max,i,i]=0
        }
      }
      
      
      # Then calculate the AR coefficients for that order
      arcoef=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
        if(xarorder[i]==0){
          return(NULL)
        }
        if(is.matrix(lh)){
          arcoef=DLpencoef(x[,i], lag.max = xarorder[i], na.action=na.action,penalized=penalized,lh=lh[,i],return.mat=TRUE)
        }
        else{
          arcoef=DLpencoef(x[,i], lag.max = xarorder[i], na.action=na.action,penalized=penalized,lh=lh,return.mat=TRUE)
        }
        if(xarorder[i]==1){
          arcoef=matrix(arcoef,nrow=1)
        }
        return(arcoef)
      },simplify=FALSE) # list of length nser with each element being xarorder[i] x xarorder[i]
      
      # Then convert those into the acf
      acf=acf(x,lag.max=lag.max,type=type,plot=FALSE,na.action=na.action,demean=demean,penalized=penalized,lh=lh)
      xacf=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
        if(xarorder[i]==0){ # Independence so no AR terms
          xacf=rep(0,lag.max)
          return(xacf)
        }
        xacf=xpacf$acf[,i,i] #initialize
        if(lag.max<=1){
          return(xacf)
        }
        prodacf=cumprod((1-xpacf$acf[1:(lag.max-1),i,i]^2)) * xpacf$acf[2:lag.max,i,i]
        for(j in 2:min(lag.max,xarorder[i]+1)){
          xacf[j]=prodacf[j-1] + arcoef[[i]][1:(j-1),j-1]%*%xacf[(j-1):1]
        }
        if(xarorder[i]==lag.max){return(xacf)}
        for(j in (xarorder[i]+2):lag.max){
          xacf[j]=arcoef[[i]][,xarorder[i]]%*%xacf[(j-1):(j-xarorder[i])]
        }
        return(xacf)
      })
      if(nser==1){
        xacf=array(xacf,dim=c(lag.max,1))
      }
      else if(lag.max==1){
        xacf=array(xacf,dim=c(1,nser))
      }
      
      for(i in 1:nser){
        acf$acf[-1,i,i]=xacf[,i]
      }
      acf$penalized=penalized
      acf$estimate="invertpacf"
     }
    else{stop("The estimate argument can only take values 'direct' or 'invertpacf'.")}

    if(plot){
      plot(acf, ...)
      invisible(acf)
    }
    else return(acf)
}
