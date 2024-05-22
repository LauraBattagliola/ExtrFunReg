
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#' Scaling Factor for the Local Bandwidth
#'
#' Computes the scaling factor for the local bandwidth based on a given extremile level.
#'
#' @param tau Numeric, representing the extremile level, which is used to tailor the
#'            scaling factor. This should be a value between 0 and 1.
#'
#' @return Numeric, the computed scaling factor for the local bandwidth. This factor
#'         is used to adjust bandwidth found via cross-validation.
#'
#' @author Maria Laura Battagliola
#'
#' @references Battagliola, M. L., and Bladt, M. (2024).
#'	Extremile scalar-on-function regression with application to climate scenarios.
#'	\emph{Under review}.
#' @references Daouia, A., Gijbels I. and Stupfler, G. (2022).
#'	Extremile regression.
#'	\emph{Journal of the American Statistical Association} 117 1579–1586.
#'
#' @examples
#' tau_value <- 0.9
#' scaling_factor <- bandwith.scaling(tau = tau_value)
#' print(scaling_factor)
#'
bandwith.scaling <- function(tau){

  f1 <- function(t, tau) {
    return(sapply(t, function(ti) qnorm(ti) * J_tau(tau = tau, t = ti)))
  }

  f2 <- function(t, tau, mu) {
    return(sapply(t, function(ti) (qnorm(ti) - mu)^2 * J_tau(tau = tau, t = ti)))
  }


  mu <- integrate(f1, lower=0, upper=1, tau=tau)$value

  V <- integrate(f2,  lower=0, upper=1, tau=tau, mu=mu)$value

return((4*tau*(1-tau)*V*(J_tau(tau,tau))^2)^(1/5))

}

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

## Find theoretical extremile of Gaussian centered in 0 (used in simulation study)
Jtau0 <- function(t ,tau){
  if(tau>0 & tau<=0.5) {
    s_tau <- log(0.5)/log(1-tau)
    return(s_tau*(1-t)^(s_tau-1))
  }
  if(tau>0.5 & tau<=1){
    r_tau <- log(0.5)/log(tau)
    return(r_tau*t^(r_tau-1))
  }

}

true.xi.Gauss <- function(t,tau, sd){
  return(qnorm(t, sd=sd.eps)*Jtau0(t,tau))
}

true.xi.std.Gauss <- function(t,tau){
  return(qnorm(t, sd=1)*Jtau0(t,tau))
}

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

## Trapeziodal integration
trapezoidal_integration <- function(f_values, a, b, n) {

  h <- (b - a) / n

  if (length(f_values) != n + 1) {
    stop("The number of function values must be n + 1")
  }

  integral_approx <- (h / 2) * (f_values[1] + 2 * sum(f_values[2:n]) + f_values[n + 1])

  return(integral_approx)
}

## Norm in L2
L2norm <- function(x, t){
  S <- length(t)
  return(sqrt(trapezoidal_integration(f_values = x^2, a=t[1], b=t[S], n=S-1)))
}

## Inner prod in L2
L2inprod <- function(x,w,t){
  S <- length(t)
  return(trapezoidal_integration(f_values = x*w, a=t[1], b=t[S], n=S-1))
}
