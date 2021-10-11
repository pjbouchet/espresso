lower = 60
upper = 215
dose.range <- seq(lower, upper, 0.5)
mu = 140
phi = 25
sigma = 20

psi = 0.2
omega = 1.5
npts = 50
nu1 = 130
nu2 = 155
alpha = mean(c(nu1, nu2))
# tau = sqrt(phi^2+sigma^2)
tau = 22

p.response.individ <- matrix(data = 0, nrow = npts, ncol = length(dose.range))

pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))],
                              mean = psi,  sd = omega))
                 
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
lines(dose.range, truncnorm::ptruncnorm(q = dose.range, a = lower, b = upper, mean = mu, sd = sqrt(phi^2+sigma^2)))
lines(dose.range, p.response, col = "blue")
