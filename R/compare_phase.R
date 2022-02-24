#' Compare dose-response functions
#'
#' Plot monophasic and biphasic dose-response functions based on a set of input parameters.
#' 
#' @export
#' 
#' 
#' @author Phil J. Bouchet

compare_phase <- function(lower = 60,
                          upper = 215,
                          mu = 140,
                          phi = 25,
                          sigma = 20,
                          psi = 0.5,
                          omega = 0.5,
                          nu1 = 90,
                          nu2 = 178,
                          alpha = 140,
                          tau = 10,
                          npts = 20,
                          linecol = c("steelblue", "firebrick")){
  
  dose.range <- seq(lower, upper, 0.5)
  
  par(mfrow = c(1, 2))
  
  # BIPHASIC
  
  d1 <- truncnorm::dtruncnorm(x = dose.range, mean = nu1, sd = tau, a = lower, b = alpha)
  d2 <- truncnorm::dtruncnorm(x = dose.range, mean = nu2, sd = tau, a = alpha, b = upper)
  plot(dose.range, d1, col = linecol[1], type = 'l', ylab = 'density', xlab = "Acoustic dose (dB)")
  lines(dose.range, d2, col = linecol[1], lty = 2, ylab = 'density', xlab = "Acoustic dose (dB)")
  
  d3 <- truncnorm::dtruncnorm(x = dose.range, mean = mu, sd = sqrt(phi^2+sigma^2), a = lower, b = upper)
  lines(dose.range, d3, col = linecol[2], lty = 1, ylab = 'density', xlab = "Acoustic dose (dB)")
  
  p.response.individ <- matrix(data = 0, nrow = npts, ncol = length(dose.range))
  pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))], mean = psi, sd = omega))
  
  # Dose-response curve corresponding to lower mixture component
  p.response.individ.lower <- 
    truncnorm::ptruncnorm(q = dose.range, 
                          a = lower, 
                          b = alpha,
                          mean = nu1, 
                          sd = tau)
  
  # Dose-response curve corresponding to higher mixture component
  p.response.individ.upper <- 
    truncnorm::ptruncnorm(q = dose.range, 
                          a = alpha,
                          b = upper,
                          mean = nu2,
                          sd = tau)
  
  for (j in 1:npts) { p.response.individ[j, ] <- pi.individ[j] * p.response.individ.lower +
    (1 - pi.individ[j]) * p.response.individ.upper}
  
  p.response <- apply(X = p.response.individ, MARGIN = 2, FUN = mean)
  
  plot(dose.range, seq(0, 1, length.out = length(dose.range)), type ='n', xlab = "Acoustic dose (dB)", ylab = "p(response)")
  lines(dose.range, p.response, col = linecol[1])
  
  # MONOPHASIC
  lines(dose.range, truncnorm::ptruncnorm(q = dose.range, a = lower, b = upper, mean = mu, sd = sqrt(phi^2+sigma^2)), col = linecol[2])
  
  
  
}