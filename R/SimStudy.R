#' Conduct a Simulation Study for Extremile Regression
#'
#' This function performs a simulation study to assess the effectiveness of extremile
#' functional regression. It generates datasets, computes extremiles, and optionally compares
#' these to results obtained from functional linear regression when certain conditions are met.
#'
#' @param n Integer, the number of observations to simulate in each dataset.
#' @param S Integer, the number of points in regular grid over [0,1].
#' @param gamma Numeric, the heteroskedasticity parameter.
#' @param tau.gen Numeric, the extremile level used for generating data.
#' @param sd.eps Numeric, standard deviation of the error term, defaults to 0.25.
#' @param tau.est Numeric, the extremile level used for estimation in the simulation.
#' @param percent Numeric, the percentage of variance explained when computing the eigenfunctions, defaults to 0.95.
#' @param B Integer, the number of simulation iterations.
#' @param my_seed Integer, seed for random number generation to ensure reproducibility.
#' @param hFMat Matrix, optional, pre-specified bandwidths for each iteration. If NULL,
#'        bandwidths are calculated within the function. Defaults to NULL.
#' @param comparison.fllr Logical, indicates whether to perform a comparison with
#'        functional linear regression. Only conducted when \code{tau.est} is 0.5. Defaults to FALSE.
#' @param kappa Numeric, a scaling factor used in the bandwidth selection, defaults to 1.
#' @param n.hF Integer, the number of points for the bandwidth function, defaults to 20.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{True extremiles} Matrix of true extremiles used in data generation.
#'   \item \code{Estimated extremiles} Matrix of extremiles estimated from the simulated data.
#'   \item \code{MSE} Vector of mean squared errors calculated for each simulation.
#'   \item \code{CDF bandwidths} Matrix of calculated bandwidths for each simulation.
#'   \item \code{Estimated means} (Optional) Matrix of estimated means from functional linear regression,
#'         included only if \code{comparison.fllr} is TRUE and \code{tau.est} is 0.5.
#' }
#'
#' @details
#' The function iteratively generates data using \code{simData}, then applies \code{ExtrFunReg}
#' to estimate extremiles. If specified, it also compares these estimates to those obtained via
#' functional linear regression (\code{fllr}) when \code{tau.est} is 0.5. This function is useful
#' for evaluating the performance of estimation methods under varying simulation conditions.
#'
#' @author Maria Laura Battagliola
#'
#' @references Battagliola, M. L., and Bladt, M. (2024).
#'	Extremile scalar-on-function regression with application to climate scenarios.
#'	\emph{Under review}.
#'
#' @examples
#' results <- SimStudy(n = 50, S = 100, gamma = 0.5, tau.gen = 0.9, tau.est = 0.9, B = 10, my_seed = 123)
#' print(results$MSE)
#'


SimStudy <- function(n, S, gamma, tau.gen, sd.eps = 0.25, tau.est, percent=0.95, B, my_seed, hFMat=NULL, comparison.fllr=FALSE, kappa=1, n.hF=20){

  extr.true <-  extr.est <- hF.mat <-  matrix(NA, B, n)
  set.seed(my_seed)
  my_seeds <- round(runif(B,min=1, max=10^6), digits=0)
  if(tau.est==0.5 & comparison.fllr) mean.est <- matrix(NA, B, n)
  for(b in 1:B){

    if (b %% 10 == 0) {
      cat("Running simulation iteration:", b, "of", B, "\n")
    }

    dat <- simData(n=n, tau=tau.gen, S=S, gamma=gamma, myseed=my_seeds[b], sd.eps=sd.eps )
    LEARN = dat[[2]]$X
    PRED = LEARN
    Responses = dat[[2]]$Y
    extr.true[b,]=dat[[2]]$xi
    s = dat[[3]]
    if(is.null(hFMat)){
      res.efr <- ExtrFunReg(Responses, LEARN, PRED, s, tau=tau.est, percent = percent, kappa=kappa, n.hF=n.hF)
    }
    if(!is.null(hFMat)){
      res.efr <- ExtrFunReg(Responses, LEARN, PRED, s, tau=tau.est, percent = percent, hF = hFMat[b,], kappa=kappa, n.hF=n.hF)
    }
    extr.est[b,] <- res.efr$Estimated.extremiles
    hF.mat[b,] <- res.efr$hF
    if(tau.est==0.5 & comparison.fllr){
      res.fllr <- fllr(Responses, LEARN, PRED, percent = percent)
      mean.est[b,] <- res.fllr$Estimated.responses
    }
  }

  res <- list()
  res[[1]] <- extr.true
  res[[2]] <- extr.est
  res[[3]] <- rowMeans((extr.true - extr.est)^2)
  res[[4]] <- hF.mat
  names(res) <- c("True extremiles", "Estimated extremiles", "MSE", "CDF bandwidths")

  if(tau.est==0.5 & comparison.fllr){
  res[[5]] <- mean.est
  names(res) <- c("True extremiles", "Estimated extremiles", "MSE", "CDF bandwidths", "Estimated means")
  }

  return(res)
}

#------------------------------------------------------------------------
SimStudy.fllr <- function(n, S, gamma, tau.gen, sd.eps = 0.25, tau.est, percent=0.95, B, my_seed){

  mean.est <-  matrix(NA, B, n)
  set.seed(my_seed)
  my_seeds <- round(runif(B,min=1, max=10^6), digits=0)

  for(b in 1:B){

    if (b %% 10 == 0) {
      cat("Running simulation iteration:", b, "of", B, "\n")
    }

    dat <- simData(n=n, tau=tau.gen, S=S, gamma=gamma, myseed=my_seeds[b], sd.eps=sd.eps )
    LEARN = dat[[2]]$X
    PRED = LEARN
    Responses = dat[[2]]$Y
    s = dat[[3]]

      res.fllr <- fllr(Responses, LEARN, PRED, percent = percent)
      mean.est[b,] <- res.fllr$Estimated.responses

  }

  res <- mean.est

  return(res)
}


