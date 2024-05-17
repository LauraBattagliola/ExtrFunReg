#' Simulate Functional Data for Extremile Regression Models
#'
#'
#' @param n Integer, the number of observations to simulate.
#' @param tau Numeric, extremile level, used in the calculation of certain model components.
#' @param gamma Numeric, the heteroskedasticity parameter, defaults to 0.
#' @param S Integer, the number of points in regular grid over [0,1] to evaluate the functions over.
#' @param myseed Integer, seed for random number generation to ensure reproducibility, defaults to 1234.
#' @param sd.eps Numeric, standard deviation of the error term, defaults to 0.25.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{tau} Extremile level used in the simulation.
#'   \item \code{data} Data frame containing the response variable \code{Y}, the true extremile \code{xi},
#'         and the matrix of regressors \code{X}.
#'   \item \code{grid} The sequence of points in the functional coordinate used to construct \code{X}.
#' }
#'
#' @details
#' The function constructs the functional data \code{X} by combining sinusoidal basis functions
#' with random coefficients. The response \code{Y} is modeled as the integral of the product
#' of \code{X} and a functional coefficient \code{beta(s)} plus an error term that is scaled
#' by both a random component and the \code{gamma} parameter.
#'
#' @author Maria Laura Battagliola
#'
#' @references Battagliola, M. L., and Bladt, M. (2024).
#'	Extremile scalar-on-function regression with application to climate scenarios.
#'	\emph{Under review}.
#'
#' @examples
#' simulated_data <- simData(n = 100, tau = 0.9, gamma = 0.5, S = 200)
#' par(mfrow=c(1,2))
#' plot(simulated_data$data$Y, simulated_data$data$xi, main = "Response vs. Extremile")
#' abline(a=0, b=1)
#' matplot(simulated_data$grid,t(simulated_data$data$X),type="l", main="Functional covariates")
#'


simData <- function(n, tau, gamma=0, S, myseed=1234, sd.eps=0.25){
  set.seed(myseed)
  my_seeds <- round(runif(5,min=1, max=10^6), digits=0)

  ## Funtional coordinate in [0,1] wlog
  s <- seq(0,1,length.out=S)

  ## Create X
  # Basis
  phi1 <- sin(pi*s)
  phi2 <- cos(10*pi*s)
  phi3 <- sin(30*pi*s)
  phi4 <- cos(40*pi*s)
  phi5 <- s


  # Coefficients
  set.seed(my_seeds[1])
  c0 <- rnorm(n, mean=0, sd=1)
  set.seed(my_seeds[2])
  c1 <- rnorm(n, mean=2, sd=0.25)
  set.seed(my_seeds[3])
  c2 <- rnorm(n, mean=0, sd=0.25)
  set.seed(my_seeds[4])
  c3 <- rnorm(n, mean=0, sd=0.05)
  set.seed(my_seeds[5])
  c4 <- rnorm(n, mean=0, sd=0.05)


  phi <- cbind(rep(1,S), phi1, phi2, phi3, phi4)
  c<- rbind(c0, c1, c2, c3, c4)
  X <- t(phi%*%c)

  ## Functional coefficient beta(s)
  beta <- function(s) return(2*cos(2*pi*s))
  ## Integral of X and beta
  int.Xbeta <- apply(X,1,L2inprod,  w=unlist(lapply(s,beta)), t=s)


  ## Epsilon
  eps <- rnorm(n, sd=sd.eps)
  ## Extremile of epsilon
  xi.eps <- integrate(true.xi.std.Gauss, lower=0, upper=1, tau=tau)$value
  ## Compute sigma
  sigma <- 1 + gamma*apply(X, 1, trapezoidal_integration, a=0, b=1, n=S-1)
  ## Compute response Y
  Y <-int.Xbeta + eps*sigma
  ## Compute true extremile
  xi.Y <-  int.Xbeta + xi.eps*sigma*sd.eps

  ## Put results in in data.frame
  data <- data.frame(Y=Y,xi = xi.Y)
  data$X <- X

  res <- list()
  res[[1]] <- tau
  res[[2]] <- data
  res[[3]] <- s
  names(res) <- c("tau", "data", "grid")
  return(res)
}

