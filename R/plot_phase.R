#' Compare dose-response functions
#'
#' Plot monophasic and biphasic dose-response functions based on a set of input parameters.
#' 
#' @export
#' 
#' 
#' @author Phil J. Bouchet

plot_phase <- function(rjdat = NULL,
                       plot.mono = TRUE,
                       plot.bi = TRUE,
                       lower = 60,
                       upper = 215,
                       mu = 140,
                       phi = 25,
                       sigma = 20,
                       psi = 0.5,
                       omega = 0.5,
                       nu1 = 90,
                       nu2 = 178,
                       alpha = 140,
                       tau1 = 10,
                       tau2 = 10,
                       npts = 20,
                       linecol = c("steelblue", "firebrick")){
  
  x.lab <- "Acoustic dose (dB)"
  y.lab <- "p(response)"
  
  if(!is.null(rjdat)){
    
    lower <- rjdat$dat$param$dose.range[1]
    upper <- rjdat$dat$param$dose.range[2]
    
    rjtrace <- do.call(rbind, rjdat$trace)
    
    if("mu.1" %in% colnames(rjtrace)) mu <- mean(rjtrace[, "mu.1"])
    if("phi" %in% colnames(rjtrace)) phi <- mean(rjtrace[, "phi"])
    if("sigma" %in% colnames(rjtrace)) sigma <- mean(rjtrace[, "sigma"])
    if("psi" %in% colnames(rjtrace)) psi <- mean(rjtrace[, "psi"])
    if("omega" %in% colnames(rjtrace)) omega <- mean(rjtrace[, "omega"])
    if("nu.lower.1" %in% colnames(rjtrace)) nu1 <- mean(rjtrace[, "nu.lower.1"])
    if("nu.upper.1" %in% colnames(rjtrace)) nu2 <- mean(rjtrace[, "nu.upper.1"])
    if("alpha.1" %in% colnames(rjtrace)) alpha <- mean(rjtrace[, "alpha.1"])
    if("tau.lower" %in% colnames(rjtrace)) tau1 <- mean(rjtrace[, "tau.lower"])
    if("tau.upper" %in% colnames(rjtrace)) tau2 <- mean(rjtrace[, "tau.upper"])
    
  }
  
  bp <- function(L, U, nu1, tau1, nu2, tau2, psi, omega, by = by){
    
    B <- 1E5
    
    # Alpha uniform between nu1 and nu2
    alpha <- runif(B, nu1, nu2)
    
    psi_i <- rnorm(B, psi, omega)
    pi_ij <- pnorm(psi_i)
    k_ij <- rbinom(B, 1, pi_ij)
    mu_1i <- truncnorm::rtruncnorm(B, L, alpha, nu1, tau1)
    mu_2i <- truncnorm::rtruncnorm(B, alpha, U, nu2, tau2)
    t_ij <- ifelse(k_ij, mu_1i, mu_2i)
    
    x <- seq(L, U, by = by)
    n.x <- length(x)
    y <- numeric(n.x)
    for(i in 1:n.x){
      y[i] <- sum(t_ij <= x[i]) / B
    }
    return(y)
  }
  
  dose.range <- seq(lower, upper, 0.5)
  
  par(mfrow = c(1, 2))
  
  # BIPHASIC
  if(plot.bi){
    d1 <- truncnorm::dtruncnorm(x = dose.range, mean = nu1, sd = tau1, a = lower, b = alpha)
    d2 <- truncnorm::dtruncnorm(x = dose.range, mean = nu2, sd = tau2, a = alpha, b = upper)
    plot(dose.range, d1, col = linecol[1], type = 'l', ylab = 'density', xlab = x.lab)
    lines(dose.range, d2, col = linecol[1], lty = 2, ylab = 'density', xlab = x.lab)
  }
  
  if(plot.mono){
    
    d3 <- truncnorm::dtruncnorm(x = dose.range, mean = mu, sd = sqrt(phi^2+sigma^2), a = lower, b = upper)
    if(!plot.bi)  plot(dose.range, d3, col = linecol[2], lty = 1, type = 'l', ylab = 'density', xlab = x.lab) else 
      lines(dose.range, d3, col = linecol[2], lty = 1, ylab = 'density', xlab = x.lab)
  }
  
  if(plot.bi){
    
    p.response.individ <- matrix(data = 0, nrow = npts, ncol = length(dose.range))
    pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))], mean = psi, sd = omega))
    
    # Dose-response curve corresponding to lower mixture component
    p.response.individ.lower <- 
      truncnorm::ptruncnorm(q = dose.range, 
                            a = lower, 
                            b = alpha,
                            mean = nu1, 
                            sd = tau1)
    
    # Dose-response curve corresponding to higher mixture component
    p.response.individ.upper <- 
      truncnorm::ptruncnorm(q = dose.range, 
                            a = alpha,
                            b = upper,
                            mean = nu2,
                            sd = tau2)
    
    for (j in 1:npts) { p.response.individ[j, ] <- pi.individ[j] * p.response.individ.lower +
      (1 - pi.individ[j]) * p.response.individ.upper}
    
    p.response <- apply(X = p.response.individ, MARGIN = 2, FUN = mean)
    
    plot(dose.range, seq(0, 1, length.out = length(dose.range)), type ='n', xlab = x.lab, ylab = y.lab)
    lines(dose.range, bp(lower, upper, nu1, tau1, nu2, tau2, psi, omega, by = 0.5), col = linecol[1], lwd = 1.5)
    lines(dose.range, p.response, col = "grey")
  }

  # MONOPHASIC
  if(plot.mono){
    if(!plot.bi) plot(dose.range, truncnorm::ptruncnorm(q = dose.range, a = lower, b = upper, mean = mu, sd = sqrt(phi^2+sigma^2)), col = linecol[2], type = 'l', xlab = x.lab, ylab = y.lab) else lines(dose.range, truncnorm::ptruncnorm(q = dose.range, a = lower, b = upper, mean = mu, sd = sqrt(phi^2+sigma^2)), col = linecol[2])
  }
  
  legend("topleft", bty = "n", inset = c(0.01, 0.03), c(ifelse(plot.bi, "Biphasic", ""), ifelse(plot.mono, "Monophasic", "")), col = linecol, lty = c(ifelse(plot.bi, 1, NA), ifelse(plot.mono, 1, NA)))
  
}