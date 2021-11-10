# Test 1 – when nu1 and nu2 are close, would expect 0.5/0.5

mydat <- simulate_data(biphasic = TRUE,
                       n.species = 1, 
                       n.whales = 12,
                       min.trials = 2,
                       max.trials = 4,
                       psi = 0,
                       omega = 1,
                       alpha = 140,
                       nu = list(c(139, 141)),
                       tau = c(32, 32),
                       Lc = c(60, 60.5),
                       Rc = c(211, 215))

# Test 2 - expect favouring monophasic when tau is big, even when nu1 and nu2 are far apart

mydat <- simulate_data(biphasic = TRUE,
                       n.species = 1, 
                       n.whales = 12,
                       min.trials = 2,
                       max.trials = 4,
                       psi = 0,
                       omega = 1,
                       alpha = 140,
                       nu = list(c(100, 180)),
                       tau = c(32, 32),
                       Lc = c(60, 60.5),
                       Rc = c(211, 215))

# Test 3 - expect favouring biphasic when tau is small

mydat <- simulate_data(biphasic = TRUE,
                       n.species = 1, 
                       n.whales = 12,
                       min.trials = 2,
                       max.trials = 4,
                       psi = 0,
                       omega = 1,
                       alpha = 140,
                       nu = list(c(100, 180)),
                       tau = c(10, 10),
                       Lc = c(60, 60.5),
                       Rc = c(211, 215))

# Test 4 - when few animals in one mixture, might favour monophasic

mydat <- simulate_data(biphasic = TRUE,
                       n.species = 1, 
                       n.whales = 12,
                       min.trials = 2,
                       max.trials = 4,
                       psi = -2,
                       omega = 1,
                       alpha = 140,
                       nu = list(c(100, 180)),
                       tau = c(10, 10),
                       Lc = c(60, 60.5),
                       Rc = c(211, 215))