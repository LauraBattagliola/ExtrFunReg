### Functions involved in the optimal local bandwith choice for the conditional CDF estimator

diff.F <- function(D, Y, h, h1, n, kernCDFI , n.y = 200){
  Y.min <- Y[which.min(Y)]
  Y.max <- Y[which.max(Y)]
  y.seq <- seq(Y.min, Y.max, length.out=n.y+1)
  Fhat1 <- sapply(y.seq, cdf_ksm_single, D=D, Y=Y, h=rep(h1,n), n=n,kernCDFI=kernCDFI)
  Fhat <- sapply(y.seq, cdf_ksm_single, D=D, Y=Y, h=rep(h,n), n=n,kernCDFI=kernCDFI)

  res <- list()
  res[[1]] <- Fhat
  res[[2]] <- Fhat1
  res[[3]] <- Fhat1 - Fhat
  names(res) <- c("Fhat", "Fhat1", "diff")
  return(res)
}

small_ball_diff <- function(D_row, h, s){
  ind <- which(D_row<=h)
  return(length(ind)/length(D_row))
}

V_diff <- function(D_row, h, kappa, s){
  n <- length(D_row)
  phi.hat <- small_ball_diff(D_row,h, s)
  if(phi.hat!=0) return(kappa*(log(n)/(n*phi.hat)))
  else return(Inf)
}

#--------------------------------------------------------------------------------------

#' Optimal Bandwidth Calculation for Conditional CDF Estimator
#'
#'
#' @param X          Matrix of functional input data where each row represents an observation
#'                   and each column represents the evaluation of the function on a certain point of the grid.
#' @param s          A vector indicating the common grid on which the functions are evaluated
#' @param Y          Response variable corresponding to the rows of X.
#' @param n.h        Number of bandwidths to consider in the optimization process.
#'                   Default is 20.
#' @param kernCDFI   Kernel choice indicated by an integer used in the difference
#'                   function calculation. Default is 3, corresponding to Epanechnikov
#' @param kappa      A scaling factor which controls the trade-off between estimated bias and variance. Default is 1.
#'
#' @return A vector of optimal bandwidths, one for each observation in the dataset.
#'
#' @author Maria Laura Battagliola
#'
#' @references Battagliola, M. L., and Bladt, M. (2024).
#'	Extremile scalar-on-function regression with application to climate scenarios.
#'	\emph{Under review}.
#' @references Chagny, G., and Roche, A. (2014).
#'	Adaptive and minimax estimation of the cumulative distribution function given a functional covariate.
#'	\emph{Electronic Journal of Statistics} 8 2352 – 2404.
#'
#' @examples
#' # Simulating some data
#' n <- 200
#' S <- 100
#' tau <- 0.5
#' simData <- simData(n=n, tau=tau, S=S)
#' X <- simData[[2]]$X
#' Y <- simData[[2]]$Y
#'
#' # Calculating optimal bandwidths for the conditional CDF
#' s <- seq(0,1,length.out=S)
#' n.h <- 10
#' optimal_bandwidths <- optband.F(X, s, Y, n.h)
#' print(optimal_bandwidths)
#'
#'
optband.F<- function(X, s, Y, n.h=20, kernCDFI=3, kappa=1){

  n <- dim(X)[1]
  h.opt <- rep(NA,n)

  D = as.matrix(parDist(X, method = "euclidean")) / sqrt(length(s) - 1)

  h.min <- min(D[which(D!=0)])
  h.max <- max(D)
  h.vec <- seq(h.min, h.max, length.out=n.h)

  for(i in 1:n){

    vi <- rep(NA,n.h)
    count <- 1

    for(j in 1:n.h){
      vi[j] <- V_diff(D_row=D[i,], h=h.vec[j], s=s, kappa=kappa)
    }


    d <- D[i,-i]
    Yi <- Y[-i]
    A.vec <- matrix(0,n.h)
    Y.min <- Yi[which.min(Yi)]
    Y.max <- Yi[which.max(Yi)]


    ai <- rep(NA, n.h)
    count <- 1

    for(h in h.vec){
      for(j in 1:n.h){

        hi <- max(h, h.vec[j])

        dF<- diff.F( D=d, Y=Y, h=hi, h1=h.vec[j], n=n,kernCDFI=kernCDFI)[[3]]
        if(sum(which(dF!=0))==0) norm.diff <- 0
        else norm.diff <- trapezoidal_integration(dF^2, Y.min, Y.max, length(dF)-1)

        a <- norm.diff - vi[j]

        if(a>=0) A.vec[j] <- a

      }

      ai[count]<- A.vec[which.max(A.vec)]
      count <- count+1
    }

    h.opt[i] <- h.vec[which.min(ai+vi)]
  }

  return(h.opt)

}



