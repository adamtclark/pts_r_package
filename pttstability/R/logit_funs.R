#' Logit
#'
#' Returns the logit transformation of x
#' @param x a number, vector, matrix, etc. to be transformed from (0, 1) to (-inf inf) by the logit transform
#' @keywords logit
#' @return transformed result - impossible values are replaced with NA, without warnings
#' @export

logit<-function(x) {
  suppressWarnings(res<-(-log(1/x-1)))
  res[!is.finite(res)]<-NA
  res
}


#' Inverse logit
#'
#' Returns the inverse logit transformation of x
#' @param x a number, vector, matrix, etc. to be transformed from (-inf, inf) to (0 1) by the inverse logit transform
#' @keywords logit
#' @return transformed result
#' @export

ilogit<-function(x) {
  1/(1+exp(-x))
}

