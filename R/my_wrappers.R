llsm_cv_singleR <- function(C,D,Y,H=NULL,kNN=NULL,nCV=NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  ## input:
  # C:      basis coefficients of the X functions (matrix n*J)
  # D:      matrix of distances of X from X functions (matrix n*n)
  # Y:      response vector (vector n)
  # H:      sequence of bandwidths to cross-validate for
  # in CV for derivative, only such H are taken so that at least 10 functions are in
  # the h-neighborhood of any function x / otherwise we get very unstable, unreliable
  # results
  # kNN:    proportion of k-nearest neighbors to be used in bandwidth definition (in the interval (0,1))
  # quantile.type: type option in quantile function (see ?quantile)

  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(H) & is.null(kNN)) stop("Neither candidate bandwidths nor nearest neighbour proportions specified.")
  if((!is.null(H)) & (!is.null(kNN))) stop("Specify either bandwidths or nearest neighbour proportions, not both.")
  # construct full matrix of candidate bandwidths, different for each function
  if(is.null(H)) Heval = length(kNN)
  if(is.null(kNN)){ if(is.null(dim(H))) Heval = length(H) else Heval = ncol(H) }
  Hf = matrix(nrow=nrow(D),ncol=Heval) # each function has now its own sequence of candidate bw
  for(i in 1:nrow(D)){
    if(is.null(H)){
      Hf[i,] = quantile(D[i,-i],kNN,type=quantile.type)
      Hf[i,kNN==1] = Inf
    }
    if(is.null(kNN)){
      if(is.null(dim(H))) Hf[i,] = H else Hf[i,] = H[i,]
    }
  }

  Hf[!is.finite(Hf)] = -1  # -1 is the code for usual linear regression


  if(is.null(C)){
    iC = matrix(1,nrow=nrow(D))
    iCnew = 1
  } else {
    iC<-cbind(1,C)
  }
  if(is.null(nCV)) nCV = nrow(D)
  if(nCV>nrow(D)) nCV = nrow(D)

  n=nrow(D)
  J = ncol(iC)
  CV = rep(0,Heval*ncol(iC))
  CVB = rep(0,Heval*ncol(iC))

  suppressWarnings(llsm_cv_single_cpp(iC, D, Y, n, J, Hf, Heval, CV, CVB, nCV, kernI))


  CVM = matrix(CV,nrow=Heval)
  CVM[is.na(CVM)] = Inf
  CVB = matrix(CVB,nrow=Heval)
  CVB[is.na(CVB)] = Inf


  Hf[Hf==-1] = Inf # back to the usual coding for Inf
  if(!is.null(H)){ # H cross-validation
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")

    return(list(CV = CVM, CVB = CVB, kernel = kernel,
                Hf=Hf, iC=iC, Heval=Heval, D=D, J=J))
  }
  if(is.null(H)){ # kNN cross-validation
    rownames(CVM) = rownames(CVB) = paste("kNN=",round(kNN,2),sep="")
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")

    return(list(CV=CVM, CVB = CVB, kernel = kernel))
  }
}


