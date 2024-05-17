#' Extremile Scalar-on-Function Regression
#'
#' Conducts extremile functional regression to estimate conditional extremiles when the covariates are functional.
#' This function accommodates both in-sample and out-of-sample predictions if specified.
#'
#' @param Responses Vector of response values corresponding to the learning set.
#' @param LEARN Matrix of predictors in the learning set. every row corresponds to one functional covariate.
#' @param PRED Matrix of predictors in the prediction set; may be identical to LEARN.
#' @param s Vector of functional coordinates, assumed to be a regular dense grid.
#' @param tau Numeric, the extremile level of interest.
#' @param kNN.grid.length Integer, number of grid points for k-nearest neighbor calculation.
#' @param percent Numeric, percentage of variance explained in the eigenfunction expansion.
#' @param BASIS Optional matrix of basis functions; if NULL, they are estimated from the data.
#' @param J.est Logical, whether to estimate the number of basis functions; defaults to TRUE.
#' @param fullgrid Logical, whether to use a full grid; defaults to using a full grid if LEARN
#'        has 250 or fewer rows.
#' @param hF Optional pre-specified bandwidths; if NULL, they are calculated within the function.
#' @param kappa Numeric, a scaling factor used in bandwidth optimization.
#' @param n.hF Integer, the number of bandwidth points to consider in optimization.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Estimated.extremiles} Vector of estimated extremile values.
#'   \item \code{ESTIMATED.DERIV} Matrix of estimated functional derivatives.
#'   \item \code{min.crit.reg} Minimum criterion value from the regression.
#'   \item \code{knn.opt} Optimal number of nearest neighbors used.
#'   \item \code{dim.opt} Optimal number of dimensions used in the basis expansion.
#'   \item \code{hF} Bandwidths used in the regression.
#'   \item \code{Predicted.extremiles} (optional) Vector of predicted extremile values for out-of-sample data.
#'   \item \code{PREDICTED.DERIV} (optional) Matrix of predicted functional derivatives for out-of-sample data.
#' }
#'
#' @details
#' This function is a wrapper of \code{fllr} in the omonimous package.
#' \code{ExtrFunReg} estimates the conditional extremiles using local bandwidth computations for both the local linear regression and conditional CDF estimator, involved in the regression weigths.
#' In particular, the estimation of the extremiles is carried out in a similar fashion of \code{fllr}, where first the eigenfunctions are computed and then the leave-one-out cross validation procedure
#' is used to find the optimal bandwidths for the local linear mean regression. We rescale such bandwidths (with constant depending on \code{tau}) to obtain the bandwiths for the extremile regression.
#' Moreover, if one was to carry out out-of-sample prediction, the optimal bandwiths for the conditional CDF estimator computed for the estimation are used in the prediction.
#' Finally, notice that we call "derivative" the finite dimensional representation of coefficient beta in Battagliola, M. L., and Bladt, M. (2024).
#' This makes sense in the case when \code{tau}=0.5, namely when performing
#' local linear mean regression as Ferraty, F., and Nagy, S. (2021). However, we keep the same notation of \code{fllr} even for different extremile levels.
#'
#' @author Maria Laura Battagliola
#'
#' @references Battagliola, M. L., and Bladt, M. (2024).
#'	Extremile scalar-on-function regression with application to climate scenarios.
#'	\emph{Under review}.
#'	@references Ferraty, F., and Nagy, S. (2021).
#'	Scalar-on-function local linear regression and beyond. .
#'	\emph{Biometrika} 109 439-455.
#'	@references Daouia, A., Gijbels I. and Stupfler, G. (2022).
#'	Extremile regression.
#'	\emph{Journal of the American Statistical Association} 117 1579–1586.
#'	@references Chagny, G., and Roche, A. (2014).
#'	Adaptive and minimax estimation of the cumulative distribution function given a functional covariate.
#'	\emph{Electronic Journal of Statistics} 8 2352 – 2404.
#'
#' @examples
#' n <- 100
#' S <- 200
#' tau <- 0.9
#' simulated_data <- simData(n=n, S=S, tau=tau)
#' Responses <- simulated_data$data$Y
#' LEARN <- simulated_data$data$X
#' PRED <- LEARN  # for in-sample predictions
#' s <- simulated_data$grid
#' tau <- 0.95
#' result <- ExtrFunReg(Responses, LEARN, PRED, s, tau)
#' print(result$Estimated.extremiles)
#'


ExtrFunReg = function(Responses, LEARN, PRED, s, tau, kNN.grid.length = 30,
                      percent = 0.95, BASIS = NULL,
                     J.est = TRUE, fullgrid=NULL, hF=NULL, kappa=1, n.hF=20)
{
  if(is.null(fullgrid)) fullgrid = (nrow(LEARN)<=250)
  # Check whether there are additional observations for predicting
  # corresponding responses
  outsample = T
  if(identical(LEARN,PRED)) outsample = F
  ## if LEARN is the same matrix as PRED, outsample = F and PRED = LEARN
  ##
  # Estimate the second order moment operator
  gridsize = ncol(LEARN)
  if(is.null(BASIS)){
    if(outsample){
      sample.size = nrow(LEARN) + nrow(PRED)

      COV = var(rbind(LEARN,PRED))/(gridsize-1) # added in R1
    }else{
      sample.size = nrow(LEARN)

      COV = var(LEARN)/(gridsize-1) # added in R1
    }
  }
  # Install the r package "parallelDist" to use "parDist" and compute
  # distances ||X_i - X_j||^2 between functional predictors
  # require(parallelDist)
  DIST = as.matrix(parDist(LEARN, method = "euclidean")) / sqrt(gridsize - 1)

  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
  if(is.null(BASIS)){
    # Compute estimated first "dim.max" eigenfunctions
    res.eig = eigen(COV, sym = T)
    # Drop first eigenvalue and derive the approximating subspace expressing at
    # least "percent" of the remaining total variance
    cum.part.of.variance = cumsum(res.eig$values[-1]) / sum(res.eig$values[-1])
    dim.min = 1
    dim.max = min(sum(cum.part.of.variance < percent) + 2,gridsize)
    BASIS = t(res.eig$vectors[, 1:dim.max]) * sqrt(gridsize - 1)
  } else {
    # if J is estimated, dimensions are searched from one
    if(J.est) dim.min = 1 else dim.min = nrow(BASIS)
    dim.max = nrow(BASIS)
  }
  nlearn = nrow(LEARN)
  INNER.LEARN = tcrossprod(LEARN, BASIS) / (gridsize - 1)
  INNER.LEARN0=INNER.LEARN
  if(nlearn<50) warning("With less than 50 learning curves the cross-validation
                        for the kNN method may be unstable.")
  if(fullgrid) Knn.learn.ratio = c(0.02,1) else
    Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 1)
  # Build an adaptative grid of kNN (number of knearest neighbours used for
  # determining local bandwidths). If the minimum of the criterion ("CRIT.REG")
  # is reached for the bandwidth corresponding to "Knn.grid[kNN.grid.length]",
  # then the grid "Knn.grid" is updated
  kNN.iterate = TRUE # indicator stopping the iteration in kNN CV procedure
  knn.max = 2        # initiate the kNN CV procedure
  #knn.max = 1
  ratio = 1
  while(kNN.iterate){
    knn.min = min(max(round(nlearn * Knn.learn.ratio[ratio]), knn.max),nlearn-1)
    knn.max = min(max(round(nlearn * Knn.learn.ratio[ratio + 1]),
                      knn.min + kNN.grid.length  - 1),nlearn-1)
    ratio = ratio + 1
    if(ratio == length(Knn.learn.ratio)-1) kNN.iterate = FALSE # stop iteration
    if(knn.max == nlearn-1) kNN.iterate = FALSE                # stop iteration
    knn.range = knn.max - knn.min
    if(knn.range <= kNN.grid.length){
      # this now happens only if knn.max == nlearn
      Knn.grid = knn.min:knn.max
    }else{
      step = ceiling(knn.range / kNN.grid.length)
      Knn.grid = seq(from = knn.min, to = knn.max, by = step)
    }
    actual.kNN.grid.length = length(Knn.grid)
    rank.min = actual.kNN.grid.length
    ############################################################################
    # estimate regression operator with bandwidth maximizing cross-validation
    # criterion
    ############################################################################
    CRIT.REG = matrix(Inf, actual.kNN.grid.length, dim.max)

    H = matrix(nrow=nrow(INNER.LEARN),ncol=actual.kNN.grid.length)
    for(knni in 1:actual.kNN.grid.length)
      H[,knni] = DIST.SORTED.ROW.BY.ROW[, Knn.grid[knni] + 1]
    for(dim in dim.min:dim.max){
      CRIT.REG[,dim] = suppressWarnings(
        llsm_cv_singleR(INNER.LEARN[,1:dim,drop=FALSE],
                        DIST,Responses,H=H,kernel="Epanechnikov")$CV[,dim+1])
      # warnings indicate that the Cholesky decomposition failed, in case
      # when kNN is smaller than the dimension of the model. Then, the CV
      # criteria are are replaced by Inf
      CRIT.REG[Knn.grid<dim+2,dim] = Inf
      # numerical fix
      # in the case knn<dim, the weight matrix had to be singular for kernels
      # supported in [0,1]
    }
    dim.opt = (1:dim.max)[which.min(apply(CRIT.REG,2,min))]
    # plot(apply(CRIT.REG,2,min),type="b")


    if(rank.min != order(CRIT.REG[, dim.opt])[1]) break
  }
  # Determine corresponding criterion value for the regression
  min.crit.reg = min(CRIT.REG[, dim.opt])
  # Determine optimal number of k-nearest neighbours and dimension J for the
  # regression
  knn.reg = (Knn.grid)[order(CRIT.REG[, dim.opt])[1]]
  # Compute inner products between functional predictors and basis functions
  B = matrix(BASIS[1:dim.opt,],nrow=dim.opt)

  # hF <- optband.F(X=INNER.LEARN,s=s ,Y=Responses)

  INNER.LEARN = matrix(INNER.LEARN[,1:dim.opt],ncol=dim.opt)
  # Estimated.responses = llsm_leave(
  #   INNER.LEARN,INNER.LEARN,DIST,Responses,
  #   h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1]*bandwith.scaling(tau), kernel="Epanechnikov")[,1]
  #
  knn.opt = list(reg = knn.reg, deriv = knn.reg)

  if(is.null(hF)){
    hF <- optband.F(X=LEARN,s=s ,Y=Responses, kappa=kappa, n.h = n.hF)
  }


  res <- efr_leave(INNER.LEARN,INNER.LEARN,DIST,Responses,
             h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1]*bandwith.scaling(tau),
             kernel="Epanechnikov",  kernelCDF="Epanechnikov", tau=tau, hF=hF)
  xi.est <- res[,1]
  coef.est <- res[,2:dim(res)[2]]
  est.der = coef.est %*% B

  # Compute predicted responses + predicted fuctional derivative if "outsample"
  # is true (i.e. one has additional observations for the functional predictor)

  if(outsample){
    DIST = as.matrix(parDist(rbind(LEARN, PRED), method = "euclidean")) / sqrt(
      gridsize - 1)
    DIST.SORTED.PRED.BY.LEARN = t(apply(DIST[-(1:nlearn), 1:nlearn], 1, sort))
    # Compute inner products between functional predictors and basis functions
    INNER.PRED = tcrossprod(PRED, B) / (gridsize - 1)
    nout = nrow(PRED)


    res.pred <- efr(C=INNER.LEARN,Cnew=INNER.PRED,D=DIST[-(1:nlearn), 1:nlearn],Y=Responses,
               h = DIST.SORTED.PRED.BY.LEARN [, knn.reg + 1]*bandwith.scaling(tau),
               kernel="Epanechnikov",  kernelCDF="Epanechnikov", tau=tau, hF=hF)

    xi.pred <- res.pred[,1]
    coef.pred <- res.pred[,2:dim(res)[2]]
    pred.der = coef.pred %*% B


    return(list(Estimated.extremiles = xi.est,
                ESTIMATED.DERIV = est.der,
                Predicted.extremiles = xi.pred,
                PREDICTED.DERIV = pred.der,
                min.crit.reg = min.crit.reg,
                knn.opt = knn.opt,
                dim.opt = dim.opt,
                hF=hF))
  }else{
  return(list(Estimated.extremiles = xi.est,
              ESTIMATED.DERIV = est.der,
              min.crit.reg = min.crit.reg,
              knn.opt = knn.opt,
              dim.opt = dim.opt,
              hF=hF))
  }

}

