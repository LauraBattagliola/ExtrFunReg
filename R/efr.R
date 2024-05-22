efr <- function(C, Cnew, D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7,
                tau, kernelCDF=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), hF) {
  if(is.null(dim(D))){
    D = matrix(D,nrow=1)
  }
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel,krn)
  kernI = which(kernel==krn)

  krnCDF = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernelCDF = match.arg(kernelCDF,krnCDF)
  kernCDFI = which(kernelCDF==krnCDF)

  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D)) # if all functions should be used perform usual (linear) regression
  }
  # if(h==Inf) return(llsm(C,Cnew,rep(1,length(D)),Y,h=-1,allJ = allJ, kernel = kernel)) # if h==Inf return usual regression output
  h[!is.finite(h)] = -1
  if(is.null(C)){
    iC = matrix(1,nrow=ncol(D))
    iCnew = matrix(1,nrow=nrow(D),1)
  } else {
    iC<-cbind(1,C)
    iCnew = matrix(NA,nrow(D),ncol(iC))
    iCnew[,1] = 1
    iCnew[,2:ncol(iC)] = Cnew
  }


  res <- efr_single(C=iC, Cnew=iCnew, D=D, Y=Y, h=h, hF=hF,n=nrow(iC), m=nrow(D), J=ncol(iC), kernI=kernI, tau=tau, kernCDFI=kernCDFI)
  return(res)
}


efr_leave <- function(C, Cnew, D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7,
                tau, kernelCDF=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), hF) {
  if(is.null(dim(D))){
    D = matrix(D,nrow=1)
  }
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel,krn)
  kernI = which(kernel==krn)

  krnCDF = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernelCDF = match.arg(kernelCDF,krnCDF)
  kernCDFI = which(kernelCDF==krnCDF)

  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D)) # if all functions should be used perform usual (linear) regression
  }
  # if(h==Inf) return(llsm(C,Cnew,rep(1,length(D)),Y,h=-1,allJ = allJ, kernel = kernel)) # if h==Inf return usual regression output
  h[!is.finite(h)] = -1
  if(is.null(C)){
    iC = matrix(1,nrow=ncol(D))
    iCnew = matrix(1,nrow=nrow(D),1)
  } else {
    iC<-cbind(1,C)
    iCnew = matrix(NA,nrow(D),ncol(iC))
    iCnew[,1] = 1
    iCnew[,2:ncol(iC)] = Cnew
  }


  res <- efr_single_leave(iC, iCnew, D, Y, h, hF,nrow(iC), nrow(D), ncol(iC), kernI, tau, kernCDFI)
  return(res)
}
