% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtrFunReg.R
\name{ExtrFunReg}
\alias{ExtrFunReg}
\title{Extremile Scalar-on-Function Regression}
\usage{
ExtrFunReg(
  Responses,
  LEARN,
  PRED,
  s,
  tau,
  kNN.grid.length = 30,
  percent = 0.95,
  BASIS = NULL,
  J.est = TRUE,
  fullgrid = NULL,
  hF = NULL,
  kappa = 1,
  n.hF = 20
)
}
\arguments{
\item{Responses}{Vector of response values corresponding to the learning set.}

\item{LEARN}{Matrix of predictors in the learning set. every row corresponds to one functional covariate.}

\item{PRED}{Matrix of predictors in the prediction set; may be identical to LEARN.}

\item{s}{Vector of functional coordinates, assumed to be a regular dense grid.}

\item{tau}{Numeric, the extremile level of interest.}

\item{kNN.grid.length}{Integer, number of grid points for k-nearest neighbor calculation.}

\item{percent}{Numeric, percentage of variance explained in the eigenfunction expansion.}

\item{BASIS}{Optional matrix of basis functions; if NULL, they are estimated from the data.}

\item{J.est}{Logical, whether to estimate the number of basis functions; defaults to TRUE.}

\item{fullgrid}{Logical, whether to use a full grid; defaults to using a full grid if LEARN
has 250 or fewer rows.}

\item{hF}{Optional pre-specified bandwidths; if NULL, they are calculated within the function.}

\item{kappa}{Numeric, a scaling factor used in bandwidth optimization.}

\item{n.hF}{Integer, the number of bandwidth points to consider in optimization.}
}
\value{
A list containing:
\itemize{
  \item \code{Estimated.extremiles} Vector of estimated extremile values.
  \item \code{ESTIMATED.DERIV} Matrix of estimated functional derivatives.
  \item \code{min.crit.reg} Minimum criterion value from the regression.
  \item \code{knn.opt} Optimal number of nearest neighbors used.
  \item \code{dim.opt} Optimal number of dimensions used in the basis expansion.
  \item \code{hF} Bandwidths used in the regression.
  \item \code{Predicted.extremiles} (optional) Vector of predicted extremile values for out-of-sample data.
  \item \code{PREDICTED.DERIV} (optional) Matrix of predicted functional derivatives for out-of-sample data.
}
}
\description{
Conducts extremile functional regression to estimate conditional extremiles when the covariates are functional.
This function accommodates both in-sample and out-of-sample predictions if specified.
}
\details{
This function is a wrapper of \code{fllr} in the omonimous package.
\code{ExtrFunReg} estimates the conditional extremiles using local bandwidth computations for both the local linear regression and conditional CDF estimator, involved in the regression weigths.
In particular, the estimation of the extremiles is carried out in a similar fashion of \code{fllr}, where first the eigenfunctions are computed and then the leave-one-out cross validation procedure
is used to find the optimal bandwidths for the local linear mean regression. We rescale such bandwidths (with constant depending on \code{tau}) to obtain the bandwiths for the extremile regression.
Moreover, if one was to carry out out-of-sample prediction, the optimal bandwiths for the conditional CDF estimator computed for the estimation are used in the prediction.
Finally, notice that we call "derivative" the finite dimensional representation of coefficient beta in Battagliola, M. L., and Bladt, M. (2024).
This makes sense in the case when \code{tau}=0.5, namely when performing
local linear mean regression as Ferraty, F., and Nagy, S. (2021). However, we keep the same notation of \code{fllr} even for different extremile levels.
}
\examples{
n <- 100
S <- 200
tau <- 0.9
simulated_data <- simData(n=n, S=S, tau=tau)
Responses <- simulated_data$data$Y
LEARN <- simulated_data$data$X
PRED <- LEARN  # for in-sample predictions
s <- simulated_data$grid
tau <- 0.95
result <- ExtrFunReg(Responses, LEARN, PRED, s, tau)
print(result$Estimated.extremiles)

}
\references{
Battagliola, M. L., and Bladt, M. (2024).
Extremile scalar-on-function regression with application to climate scenarios.
\emph{Under review}.

Ferraty, F., and Nagy, S. (2021).
Scalar-on-function local linear regression and beyond. .
\emph{Biometrika} 109 439-455.

Daouia, A., Gijbels I. and Stupfler, G. (2022).
Extremile regression.
\emph{Journal of the American Statistical Association} 117 1579–1586.

Chagny, G., and Roche, A. (2014).
Adaptive and minimax estimation of the cumulative distribution function given a functional covariate.
\emph{Electronic Journal of Statistics} 8 2352 – 2404.
}
\author{
Maria Laura Battagliola
}
