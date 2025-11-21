######
#' @title Banded and Tapered estimation of the autocorrelation function
#'
#' @description It is well known that the default acf sample estimates are biased.
#' This function provides the banded and tapered estimation of the acf.  The function
#' is designed to be called from the main acf function via the \code{estimate="bandtap"}
#' argument.  This function does not calculate confidence intervals for the estimate, users
#' who want confidence intervals should use the \code{taperedacf} function from the \code{forecast}
#' package.
#' @param acf original biased acf (likely from the stats package).
#' @param sampleT length of the underlying dataset.
#' 
#' @details
#' The banded and tapered method identifies an appropriate cut-off, $r$, beyond which, the acf is 
#' likely noise.  Up to this point, the acf is returned as entered, between $r$ and $2r$ the 
#' estimate is tapered towards zero, after lag $2r$ the acf is set to exactly zero.  This method
#' works well for processes that have a natural cut-off such as the acf of an MA process.
#' 
#' @return A vector of the same length as \code{acf}.
#' 
#' @references 
#' McMurry, T.L. and Politis, D.N. (2010), Banded and tapered estimates for autocovariance 
#' matrices and the linear process bootstrap. \emph{Journal of Time Series Analysis}, 
#' 31: 471-482.
#' 
#' @examples
#' \dontrun{
#' ### AR(1)
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ar=0.8))
#' dataacf=acf(data,penalized=FALSE)$acf[,,1] # usual stats::acf() estimate
#' databandtap=PenalizedCorr:::bandtap(dataacf,length(data))
#' ts.plot(dataacf)
#' abline(h=0,col='grey')
#' lines(1:length(databandtap),databandtap,col='red')
#' dataacf==databandtap # first 10 are the same, 20 onwards is exactly zero
#' 
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ma=0.8))
#' dataacf=acf(data,penalized=FALSE)$acf[,,1] # usual stats::acf() estimate
#' databandtap=PenalizedCorr:::bandtap(dataacf,length(data))
#' ts.plot(dataacf)
#' lines(1:length(databandtap),databandtap,col='red')
#' dataacf==databandtap 
#' # first 2 are the same (lag 0 and lag 1), lag 2 is tapered, then zero from lag 3.
#' }
#####

bandtap <- function (acf, sampleT){
  
  lh <- 2*sqrt(log10(sampleT)/sampleT)
  #compares the absolute value of each acf with the threshold.
  thresh = abs(acf) < lh
  rle = rle(thresh)
  #finds position in the run-length encoding where a value repeats 5 times
  rle5 <- which(rle$lengths >= 5)
  if(length(rle5) == 0){
    return(acf)
  }
  r <- sum(rle$lengths[1:(rle5[which(rle$values[rle5] == TRUE)[1]]-1)])+1 # the first time atleast 5 in a row are TRUE
  w <- c(rep(1,r), (2-((r+1):(2*r)/r)), rep(0, length(acf)-(2*r)))
  return (acf*w)
}
