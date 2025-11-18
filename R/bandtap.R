######
#' @title bandtap
#'
#' @description Calculated the corrected (partial) autocorrelation function
#' @param acf incoming acf value
#' @param sampleT incoming length of x from acf funtion
#' @param ... additional arguments passed to plotting.
#' 
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' McMurray & Politis (2010) forecast package for CI
#####

bandtap <- function (acf, sampleT, ...){
  
  lh <- 2*sqrt(log10(sampleT)/sampleT)
  #compares the absolute value of each acf with the threshold.
  thresh = abs(acf) < lh
  rle = rle(thresh)
  #finds position in the run-length encoding where a value repeats 5 times
  rle5 <- which(rle$lengths == 5)
  if(length(rle5) == 0){
    return(acf)
  }
  r <- which(rle$values[rle5] == TRUE)[1]
  w <- c(rep(1,r), (2-((r+1):(2*r)/r)), rep(0, length(acf)-(2*r)))
  return (acf*w)
}
