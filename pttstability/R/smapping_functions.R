#' Make an embedding block from timeseries data
#'
#' Returns a matrix X, where columns are time-delayed embeddings of Y, with number of embeddings
#' specified by embedding dimension E.
#' See help file for the S_map_Sugihara1994 function for examples.
#' @param Y a timeseries vector from which to build the embedding.
#' @param E a positive integer, specifying the embedding dimension
#' @param lib an optional matrix of library positions, for specifying cases where
#' Y is a composite timeseries made up of multiple separate observations (e.g. spatial replicates).
#' Matrix should have two columns, with the first row in each column specifying the start of the
#' timeseries section, and the second column specifying the end.
#' @return a matrix of time-delayed embeddings
#' @export

makeblock = function(Y, E, lib = NULL) {
  if(is.null(lib)) {
    X = rep(1, length(Y))
    for(i in 0:(E-1)) {
      X = cbind(X, c(rep(NA, i+1), Y[-c((length(Y)-i):length(Y))]))
    }
  } else {
    X = NULL
    for(i in 1:nrow(lib)) {
      Ytmp = Y[lib[i,1]:lib[i,2]]
      Xtmp = rep(1, length(Ytmp))
      for(i in 0:(E-1)) {
        Xtmp = cbind(Xtmp, c(rep(NA, i+1), Ytmp[-c((length(Ytmp)-i):length(Ytmp))]))
      }
      X = rbind(X, Xtmp)
    }
  }
    
  colnames(X) = paste("col", 0:(ncol(X)-1), sep = "")
  return(X)
}

#' Process S-mapping coefficients
#'
#' Processes s-mapping coefficients from S_map_Sugihara1994 into a matrix of form C1, C2, C3, ... C0, where C0 is the intercept,
#' C1 is the current time step t, C2 is timestep t-1, C3 is timestep t-2, and so on.
#' Rows correspond to the time step used to produce the prediction, e.g. row 4 is used to calculate
#' predicted value for time step 5. This is the format expected by the EDMfun0 function.
#' See help file for the S_map_Sugihara1994 function for examples.
#' @param smap_coefs a matrix of s-map coefficients, taken from the S_map_Sugihara1994 function.
#' @return a matrix of s-mapping coefficients
#' @export

process_scof <- function(smap_coefs) {
  if(ncol(smap_coefs)>1) {
    smap_coefs_out = smap_coefs[,c(2:ncol(smap_coefs),1)]
    colnames(smap_coefs_out) <- paste("C", c(1:(ncol(smap_coefs)-1), 0), sep="")
  } else {
    colnames(smap_coefs_out) = "C0"
  }
  
  smap_coefs_out
}

#' Apply S-mapping algorithm from Sugihara 1994
#'
#' Carries out an S-mapping analysis, following the algorithm outlined in Sugihara (1994).
#' 
#' @param Y a timeseries vector from which to build the embedding.
#' @param E a positive integer, specifying the embedding dimension
#' @param theta a positive numeric scalar, specifying the nonlinearity parameter for the analysis. A value of 0 indicates
#' a fully linear analysis; higher numbers indicate greater nonlinearity.
#' @param X an optional matrix of time-delayed embeddings to use for the analysis
#' @param lib an optional matrix of library positions, for specifying cases where
#' Y is a composite timeseries made up of multiple separate observations (e.g. spatial replicates).
#' Matrix should have two columns, with the first row in each column specifying the start of the
#' timeseries section, and the second column specifying the end.
#' @param trimNA a logical specifying whether NA values should be removed from Y and X - defaults to FALSE
#' @return a list, including the timeseries used for S-mapping (Y), the delay embedding matrix used for S-mapping (X),
#' a vector of predictions (Y_hat), a matrix of S-mapping coefficients (C), the standard errors for the S-mapping
#' coefficients (C_SE), and goodness of fit metrics R-squared (R2) and root mean square error (RMSE).
#' @source Sugihara, G. (1994). Nonlinear forecasting for the classification of natural time-series. Philos. Trans. R. Soc. -Math. Phys. Eng. Sci., 348, 477–495.
#' @export
#' @examples
#' # create an example timeseries
#' n = 100
#' set.seed(1234)
#' datout<-makedynamics_general(n = n+2,
#'                              pdet=log(c(0.8,1)),
#'                              proc = -2.5,
#'                              detfun = detfun0_sin)
#' plot(datout$true, type = "l") # plot timeseries
#' Y = datout$true # extract true values
#' 
#' # run s-mapping
#' sout = S_map_Sugihara1994(Y = Y, E = 2, theta = 0.5)
#' s_coef = process_scof(sout$C) # process coefficients from the S-mapping output
#' 
#' # find best E/theta
#' fitout = data.frame(E = 1:5, theta = NA, RMSE = NA)
#'
#' for(i in 1:nrow(fitout)) {
#'    E = fitout$E[i]
#'    Ytmp = Y[-c(1:E)]
#'    optout = optimize(f = function(x) {S_map_Sugihara1994(Ytmp, E, x)$RMSE}, interval = c(0,10))
#'   
#'   fitout$theta[i] = optout$minimum # get best theta for given E
#'   fitout$RMSE[i] = optout$objective # get error
#' }
#' ps = which.min(fitout$RMSE)
#' 
#' E = fitout$E[ps] # get best E
#' theta = fitout$theta[ps] # get best theta
#' X = makeblock(Y, E) # get X for analysis
#' Y = Y[-c(1:E)] # trim NA values (corresponding to positions in X)
#' X = X[(E+1):nrow(X),] # trim NA values
#' sout = S_map_Sugihara1994(Y = Y, E = E,
#'   theta = theta, X = X) # run S-mapping for best paramter combination
#' sout$R2 # look at R-squared
#' 
#' # check fit
#' plot(sout$Y_hat, Y)
#' abline(a=0, b=1, lty=2)

S_map_Sugihara1994 = function(Y, E, theta, X = NULL, lib = NULL, trimNA = FALSE) {
  if(is.null(X)) {
    X = makeblock(Y = Y, E = E, lib = lib)
  }
  
  if(trimNA) {
    ps = which(!(is.na(rowSums(X)) | is.na(Y)))
    X = X[ps,]
    Y = Y[ps]
  }
  
  n = length(Y)
  C = matrix(nrow = n, ncol = E+1)
  A = matrix(nrow = n-1, ncol = E+1)
  B = matrix(nrow = n-1, ncol = 1)
  
  C_SE = matrix(nrow = n, ncol = E+1)
  
  for(i in 1:n) {
    d = sqrt(colSums((c(X[i,])-t(X))^2, na.rm=TRUE))
    w = exp(-theta*d/mean(d, na.rm=TRUE))
    
    A[] = (w*X)[-i,]
    B[] = (w*Y)[-i]
    
    if(sum(rowSums(is.finite(A))+is.finite(B))>0) {
      tmp = summary(lm(B~A-1))$coefficients
      C[i,] = unname(tmp[,1])
      C_SE[i,] = unname(tmp[,2])
    } else {
      C[i,] = NA
      C_SE[i,] = NA
    }
    
    # check for correctness
    #plot(B, A%*%C[i,])
    #abline(a=0, b = 1, lty = 2)
  }
  
  Y_hat = rowSums(C*X)
  
  #goodness of fit
  RSS = sum((Y-Y_hat)^2, na.rm=TRUE)
  R2 = 1-RSS/sum((Y-mean(Y, na.rm=TRUE))^2, na.rm=TRUE)
  RMSE = sqrt(RSS/n)
  return(list(Y = Y, X = X, Y_hat = Y_hat, C = C, C_SE = C_SE, R2 = R2, RMSE = RMSE))
}
