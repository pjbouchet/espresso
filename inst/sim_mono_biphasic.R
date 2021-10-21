# Strong biphasic pattern
lower = 60
upper = 215
dose.range <- seq(lower, upper, 0.5)
mu = 140
phi = 25
sigma = 20

psi = 0.5
omega = 0.5
npts = 20
nu1 = 90
nu2 = 178
alpha = 140
tau = 10

# Uncertainty between biphasic / monophasic (??)
lower = 60
upper = 215
dose.range <- seq(lower, upper, 0.5)
mu = 140
phi = 25
sigma = 20

psi = 0
omega = 2
npts = 20
nu1 = 115
nu2 = 165
alpha = 140
tau = 25

par(mfrow = c(1, 2))
d1 <- truncnorm::dtruncnorm(x = dose.range, mean = nu1, sd = tau, a = lower, b = alpha)
d2 <- truncnorm::dtruncnorm(x = dose.range, mean = nu2, sd = tau, a = alpha, b = upper)
plot(dose.range, d1, col = "lightblue", type = 'l', ylab = 'density')
lines(dose.range, d2, col = "blue", ylab = 'density')
# lines(dose.range, pnorm(psi) * d1 + (1-pnorm(psi)) * d2, ylab = 'density')

p.response.individ <- matrix(data = 0, nrow = npts, ncol = length(dose.range))

pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))], mean = psi, sd = omega))
                 
# for(j in pi.individ){
#   lines(dose.range, j * d1 + (1-j) * d2, ylab = 'density', col = "grey")
# }


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
