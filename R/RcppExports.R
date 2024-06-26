# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

metrl2_cpp <- function(A, B) {
    .Call(`_ExtrFunReg_metrl2_cpp`, A, B)
}

metrl2b_cpp <- function(A, B) {
    .Call(`_ExtrFunReg_metrl2b_cpp`, A, B)
}

kern <- function(t, k) {
    .Call(`_ExtrFunReg_kern`, t, k)
}

J_tau <- function(tau, t) {
    .Call(`_ExtrFunReg_J_tau`, tau, t)
}

kern_CDF <- function(t, k) {
    .Call(`_ExtrFunReg_kern_CDF`, t, k)
}

cdf_ksm_single <- function(D, Y, y, h, n, kernCDFI) {
    .Call(`_ExtrFunReg_cdf_ksm_single`, D, Y, y, h, n, kernCDFI)
}

multiplyMatrices <- function(A, B) {
    .Call(`_ExtrFunReg_multiplyMatrices`, A, B)
}

multiplyMatrixVector <- function(M, v) {
    .Call(`_ExtrFunReg_multiplyMatrixVector`, M, v)
}

solveLinearSystem <- function(A, b) {
    .Call(`_ExtrFunReg_solveLinearSystem`, A, b)
}

llsm_single_cpp <- function(C, Cnew, D, Y, h, n, m, J, kernI) {
    .Call(`_ExtrFunReg_llsm_single_cpp`, C, Cnew, D, Y, h, n, m, J, kernI)
}

llsm_cv_single_cpp <- function(C, D, Y, n, J, H, nH, CV, CVB, nCV, kernI) {
    invisible(.Call(`_ExtrFunReg_llsm_cv_single_cpp`, C, D, Y, n, J, H, nH, CV, CVB, nCV, kernI))
}

efr_single <- function(C, Cnew, D, Y, h, hF, n, m, J, kernI, tau, kernCDFI) {
    .Call(`_ExtrFunReg_efr_single`, C, Cnew, D, Y, h, hF, n, m, J, kernI, tau, kernCDFI)
}

efr_single_leave <- function(C, Cnew, D, Y, h, hF, n, m, J, kernI, tau, kernCDFI) {
    .Call(`_ExtrFunReg_efr_single_leave`, C, Cnew, D, Y, h, hF, n, m, J, kernI, tau, kernCDFI)
}

