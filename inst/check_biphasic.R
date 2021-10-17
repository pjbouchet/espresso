load("../biphasic/mcmat.rda")
load("../biphasic/biphasic_sim.Rdata")

# All tibble columns shown, no colouring negative numbers, etc.
options(tibble.width = Inf) 
options(pillar.neg = FALSE) 
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

rm(d_binom, rtnorm, dtnorm, alpha, glance, hexa2hex)

devtools::load_all(".")
# library(espresso)

# Comparison with biphasic_mcmc_LT2
mydat <- simulate_data(biphasic = TRUE,
                       n.species = 1, 
                       n.whales = 10,
                       max.trials = 4,
                       nu = list(c(105, 149)),
                       tau = c(20, 20),
                       psi = 0.5,
                       omega = 1,
                       alpha = 125,
                       covariates = list(range = 1, sonar = c(0, -10), exposed = c(0, 5)),
                       Lc = c(60, 60.5),
                       Rc = c(150, 165),
                       seed = 122)

mydat.config <- configure_rjMCMC(dat = mydat, 
                                 model.select = FALSE, 
                                 covariate.select = FALSE, 
                                 function.select = FALSE,
                                 biphasic = TRUE,
                                 n.rep = 100)


n.chains = 3
n.iter = 1000
n.burn = 5000
do.update = FALSE
dat = mydat.config

if(do.update){
  last.iter <- purrr::map(.x = seq_along(dat), .f = ~glance(dat = dat, which.chain = .x, f = "update"))
} else {
  last.iter <- vector(mode = "list", length = n.chains)
}
rj.list <- purrr::map(.x = seq_len(n.chains), 
                      .f = ~setup_rjMCMC(rj.input = dat,
                                         n.chains = n.chains,
                                         n.burn = n.burn,
                                         n.iter = n.iter,
                                         p.split = dat$config$move$prob[1],
                                         p.merge = dat$config$move$prob[2],
                                         moves = dat$config$move$moves,
                                         m = dat$config$move$freq,
                                         move.ratio = dat$config$move$ratio,
                                         do.update = do.update,
                                         start.values = last.iter[[.x]]))
rj <- rj.list[[sample(1:3, 1, TRUE)]]

sum(y_ij, na.rm = TRUE); sum(rj$dat$obs$y_ij, na.rm = TRUE)

rj$t.ij[1, ] <- mcmat$t.ij[1, ]
rj$mu.ij[1, , 1] <- mcmat$mu.ij$`1`[1, ]
rj$mu.ij[1, , 2] <- mcmat$mu.ij$`2`[1, ]
rj$psi.i[1, ] <- mcmat$psi.i[1, ]
rj$k.ij[1, ] <- mcmat$k.ij[1, ]
rj$pi.ij[1, ] <- mcmat$pi.ij[1, ]
rj$nu[1, ,] <- mcmat$nu[1, ]
rj$tau[1, ] <- mcmat$tau[1, ]
rj$alpha[1, ] <- mcmat$alpha[1]
rj$psi[1] <- mcmat$psi[1]
rj$omega[1] <- mcmat$omega[1]
rj$range[1, ] <- mcmat$range[1, ]
rj$exposed[1, ] <- c(0, mcmat$exposed[1, ])
rj$sonar[1, ] <- mcmat$sonar[1, ]


i = 2

if(rj$config$model.select & rj$dat$species$n > 1) { 
  
  # Propose jump and new parameters
  proposed.jump <- propose_jump(rj.obj = rj, 
                                move.type = rj$mcmc$move$m[i],
                                phase = rj$phase[i - 1])
  
  new.params <- proposal_rj(rj.obj = rj,
                            jump = proposed.jump,
                            phase = rj$phase[i - 1],
                            huelsenbeck = FALSE)
  
  # Priors
  logprior.new <- 
    uniform_prior(rj.obj = rj,
                  param.name = if(rj$phase[i - 1] == 1) "mu" else c("nu", "alpha"),
                  param = new.params)
  
  logprior.cur <- 
    uniform_prior(rj.obj = rj,
                  param.name = if(rj$phase[i - 1] == 1) "mu" else c("nu", "alpha"),
                  param = if(rj$phase[i - 1] == 1) 
                    list(out = list(mu = rj$mu[i - 1, ])) else 
                      list(out = list(nu = t(rj$nu[i - 1, ,]), 
                                      alpha = rj$alpha[i - 1, ])))
  
  # Likelihoods
  loglik.new <-
    likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
               param.name = if(rj$phase[i - 1] == 1) "mu" else
                 c("nu", "alpha"),
               rj.obj = rj,
               model = proposed.jump$model$id,
               values = new.params$out,
               included.cov = NULL,
               RJ = TRUE,
               lprod = TRUE)
  
  loglik.cur <- 
    likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
               param.name = if(rj$phase[i - 1] == 1) "mu" else
                 c("nu", "alpha"),
               rj.obj = rj,
               model = rj$mlist[[rj$current.model]], 
               values = NULL,
               included.cov = NULL,
               RJ = TRUE, 
               lprod = TRUE)
  
  # Proposal densities
  logpropdens <- propdens_rj(rj.obj = rj,
                             param = new.params,
                             jump = proposed.jump,
                             phase = rj$phase[i - 1])
  
  # Posterior ratio
  lognum <- loglik.new + 
    logprior.new + 
    logpropdens[1] + 
    proposed.jump$p[2] + 
    proposed.jump$jacobian
  
  logden <- loglik.cur + 
    logprior.cur + 
    logpropdens[2] + 
    proposed.jump$p[1]
  
  if (runif(1) < exp(lognum - logden)) {
    
    new.model <- vec_to_model(input.vector = proposed.jump$model$id,
                              sp.names = rj$dat$species$names)
    
    rj$model[i] <- rj$current.model <- new.model
    
    # Add model to list
    if(!new.model %in% names(rj$mlist)) rj$mlist[[new.model]] <- 
      proposed.jump$model$id
    
    if(!new.model %in% rj$config$clust[[1]]$model){ 
      rj$config$clust[[1]] <- 
        rbind(rj$config$clust[[1]], 
              tibble::tibble(model = new.model, p = 0, p_scale = 0)) %>% 
        dplyr::mutate(p_scale = rescale_p(p))
    }
    
    if(rj$phase[i - 1] == 1){
      rj$mu[i, ] <- new.params$out$mu
      if(rj$config$function.select){
        rj$nu[i, ,] <- rj$nu[i - 1, ,]
        rj$alpha[i, ] <- rj$alpha[i - 1, ]
      }}
    
    if(rj$phase[i - 1] == 2){
      rj$nu[i, ,] <- t(new.params$out$nu)
      rj$alpha[i, ] <- new.params$out$alpha
      if(rj$config$function.select){
        rj$mu[i, ] <- rj$mu[i - 1, ]
      }}
    
    if(i > rj$mcmc$n.burn){
      rj$accept[[paste0(ifelse(rj$phase[i - 1] == 1, "mono.", "bi."),
                        "move.", proposed.jump$type)]] <- 
        rj$accept[[paste0(ifelse(rj$phase[i - 1] == 1, "mono.", "bi."),
                          "move.", proposed.jump$type)]] + 1 }
    
  } else { # If proposal not accepted
    
    rj$model[i] <- rj$current.model <- rj$model[i - 1]
    if(rj$config$function.select | rj$phase[i - 1] == 1){
      rj$mu[i, ] <- rj$mu[i - 1, ]}
    if(rj$config$function.select | rj$phase[i - 1] == 2){
      rj$nu[i, ,] <- rj$nu[i - 1, ,]
      rj$alpha[i, ] <- rj$alpha[i - 1, ]
    }
    
  }
  
} else { # If no model selection
  
  if(rj$config$function.select | rj$phase[i - 1] == 1){
    rj$mu[i, ] <- rj$mu[i - 1, ]}
  if(rj$config$function.select | rj$phase[i - 1] == 2){
    rj$nu[i, ,] <- rj$nu[i - 1, ,]
    rj$alpha[i, ] <- rj$alpha[i - 1, ]}
  
}

if(rj$config$function.select | rj$phase[i - 1] == 1){
  rj$iter["mu"] <- rj$iter["mu"] + 1}
if(rj$config$function.select | rj$phase[i - 1] == 2){
  rj$iter["nu"] <- rj$iter["nu"] + 1
  rj$iter["alpha"] <- rj$iter["alpha"] + 1}

if(rj$dat$covariates$n > 0){
  
  if(rj$config$covariate.select){
    
    # We consider all models to be equally likely 
    # when it comes to covariates.
    # Therefore, the priors on models cancel out.
    
    # Select a covariate at random and switch it on/off
    chosen.covariate <- sample(x = rj$dat$covariates$names, size = 1)
    proposed.covariates <- rj$include.covariates[i - 1, ]
    add.remove <- 1 - rj$include.covariates[i - 1, ][chosen.covariate]
    proposed.covariates[chosen.covariate] <- add.remove
    
    if(add.remove){
      
      # Coefficients for covariate terms
      beta.params <- sapply(X = chosen.covariate,
                            FUN = function(w)
                              as.matrix(rj[[w]])[i - 1, ],
                            simplify = FALSE, USE.NAMES = TRUE)
      
      # Propose new values when a covariate is added
      beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index] <- rnorm(n = rj$dat$covariates$fL[[chosen.covariate]]$nparam, mean = 0, sd = rj$config$prop$cov)
      
      
      # Prior
      logprior <- normal_prior(
        rj.obj = rj,
        param.name = chosen.covariate,
        param = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index])
      
      # Proposal density
      logpropdens <- sum(dnorm(x = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index], mean = 0, 
                               sd = rj$config$prop$cov, log = TRUE))
      
      # Ratio of prior / proposal density
      R <- logprior - logpropdens
      
    } else {
      
      beta.params <- list()
      beta.params[[chosen.covariate]] <- NULL
      
      # Prior
      logprior <- normal_prior(rj.obj = rj,
                               param.name = chosen.covariate, 
                               param = rj[[chosen.covariate]][i - 1, rj$dat$covariates$fL[[chosen.covariate]]$index])
      
      # Proposal density
      logpropdens <- sum(dnorm(x = rj[[chosen.covariate]][i - 1,                                                       rj$dat$covariates$fL[[chosen.covariate]]$index], 
                               mean = 0, sd = rj$config$prop$cov, log = TRUE))
      
      
      # Ratio of prior / proposal density
      R <- logpropdens - logprior
      
    }
    
    # Probabilities of addition / removal
    pmove <- log(c(switch(as.character(add.remove), 
                          "1" = 1/sum(rj$include.covariates[i - 1, ] == 0), 
                          "0" = 1/sum(rj$include.covariates[i - 1, ] == 1)), 
                   1/sum(proposed.covariates == add.remove)))
    
    # The `values` argument will be the proposed values when
    # adding a covariate and the current values when removing
    # one (but these will be ignored, as the corresponding
    # indicator in included.cov will be 0).
    # included.cov corresponds to the include.covariates 
    # matrix in rj - it is filled with either 1 or 0 
    # depending on whether a covariate is included or 
    # omitted from the model.
    
    loglik.new <- 
      likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                 param = "covariates",
                 rj.obj = rj,
                 model = rj$mlist[[rj$current.model]], 
                 values = stats::setNames(list(beta.params[[chosen.covariate]]), 
                                          chosen.covariate),
                 included.cov = proposed.covariates,
                 RJ = TRUE,
                 lprod = TRUE)
    
    loglik.cur <- 
      likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                 param = "covariates",
                 rj.obj = rj,
                 model = rj$mlist[[rj$current.model]],
                 values = NULL,
                 included.cov = rj$include.covariates[i - 1, ], 
                 RJ = TRUE,
                 lprod = TRUE)
    
    # Posterior ratio
    lognum <- loglik.new + pmove[2] + R
    logden <- loglik.cur + pmove[1]
    
    
    for (k in rj$dat$covariates$names) rj[[k]][i, ] <- rj[[k]][i - 1, ]
    
    if (runif(1) < exp(lognum - logden)) {
      
      rj$include.covariates[i, ] <- proposed.covariates
      if(add.remove) rj[[chosen.covariate]][i, ] <- 
          beta.params[[chosen.covariate]]
      if(i > rj$mcmc$n.burn) 
        rj$accept[["move.covariates"]] <- 
          rj$accept[["move.covariates"]] + 1
      
    } else {
      
      rj$include.covariates[i, ] <- rj$include.covariates[i - 1, ]
      
    }
    
  } else {
    
    for (k in rj$dat$covariates$names){
      rj[[k]][i, ] <- rj[[k]][i - 1, ]}
    
  }  # End if (covariate.select)
} # End if (n.covariates > 0)

rj$iter[c("covariates", rj$dat$covariates$names)] <- 
  rj$iter[c("covariates", rj$dat$covariates$names)] + 1


if(rj$config$function.select){
  
  # We can generalize the standard 'monophasic' dose-response function
  # by allowing the threshold to come from one of two truncated normal 
  # distributions, one with lower exposure values than the other. 
  # This gives a 'biphasic' function with a context-dependent part 
  # and a dose-dependent part.
  
  ff.prop <- proposal_ff(rj.obj = rj, from.phase = rj$phase[i - 1])
  
  # Likelihoods
  loglik.new <-
    likelihood(biphasic = ifelse(ff.prop$to.phase == 1, FALSE, TRUE),
               param.name = if(ff.prop$to.phase == 1) c("mu", "mu.i", "phi", "sigma") else c("nu", "alpha", "tau1", "tau2", "omega", "psi", "mu.ij", "psi.i", "k.ij"),
               rj.obj = rj,
               model = rj$mlist[[rj$current.model]],
               values = ff.prop,
               included.cov = NULL,
               RJ = TRUE,
               lprod = TRUE)
  
  loglik.cur <-
    likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
               param.name = NULL,
               rj.obj = rj,
               model = rj$mlist[[rj$current.model]],
               values = NULL,
               included.cov = NULL,
               RJ = TRUE,
               lprod = TRUE)
  
  # Priors
  logprior.new <- ff_prior(rj.obj = rj, 
                           param = ff.prop, 
                           phase = ff.prop$to.phase)
  
  logprior.cur <- ff_prior(rj.obj = rj, 
                           param = NULL,
                           phase = rj$phase[i - 1])
  
  # Proposal densities
  logpropdens <- propdens_ff(rj.obj = rj, param = ff.prop)
  
  # Posterior ratio (Jacobian is 1 and p is 1)
  lognum <- loglik.new + logprior.new + logpropdens[[rj$phase[i - 1]]]
  logden <- loglik.cur + logprior.cur + logpropdens[[ff.prop$to.phase]]
  
  # print(paste0(ff.prop$to.phase, ":", runif(1) < exp(lognum - logden)))
  
  if (runif(1) < exp(lognum - logden)) {
    
    rj$phase[i] <- ff.prop$to.phase
    
    if(i > rj$mcmc$n.burn){
      
      if(ff.prop$to.phase == 1)
        rj$accept[["to.monophasic"]] <- rj$accept[["to.monophasic"]] + 1
      
      if(ff.prop$to.phase == 2)
        rj$accept[["to.biphasic"]] <- rj$accept[["to.biphasic"]] + 1
    }
    
    
    if(ff.prop$to.phase == 1){
      
      rj$sigma[i] <- ff.prop$sigma
      rj$mu.i[i, ] <- ff.prop$mu.i
      rj$mu[i, ] <- ff.prop$mu
      rj$phi[i] <- ff.prop$phi
      
      rj$tau[i, ] <- rj$tau[i - 1, ]
      rj$omega[i] <- rj$omega[i - 1]
      rj$psi[i] <- rj$psi[i - 1]
      rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
      rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
      rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
      rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]
      
    } else {
      
      rj$alpha[i, ] <- ff.prop$alpha
      rj$nu[i, , 1] <- ff.prop$nu[1, ]
      rj$nu[i, , 2] <- ff.prop$nu[2, ]
      rj$tau[i, ] <- ff.prop$tau
      rj$omega[i] <- ff.prop$omega
      rj$psi[i] <- ff.prop$psi
      rj$mu.ij[i, ,] <- ff.prop$mu.ij
      rj$psi.i[i, ] <- ff.prop$psi.i
      rj$k.ij[i, ] <- ff.prop$k.ij
      rj$pi.ij[i, ] <- ff.prop$pi.ij
      
      rj$sigma[i] <- rj$sigma[i - 1]
      rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
      rj$phi[i] <- rj$phi[i - 1]
    }
    
  } else {
    
    rj$phase[i] <- rj$phase[i - 1]
    
    rj$sigma[i] <- rj$sigma[i - 1]
    rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
    rj$phi[i] <- rj$phi[i - 1]
    
    rj$tau[i, ] <- rj$tau[i - 1, ]
    rj$omega[i] <- rj$omega[i - 1]
    rj$psi[i] <- rj$psi[i - 1]
    rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
    rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
    rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
    rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]
    
  }
  
} else {
  
  rj$phase[i] <- rj$phase[i - 1]
  
  if(rj$config$function.select | rj$phase[i - 1] == 1){
    rj$sigma[i] <- rj$sigma[i - 1]
    rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
    rj$phi[i] <- rj$phi[i - 1]
  }
  
  if(rj$config$function.select | rj$phase[i - 1] == 2){
    rj$tau[i, ] <- rj$tau[i - 1, ]
    rj$omega[i] <- rj$omega[i - 1]
    rj$psi[i] <- rj$psi[i - 1]
    rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
    rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
    rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
    rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]
    
  }
  
}

rj$t.ij[i, ] <- rj$t.ij[i - 1, ] # t.ij stay the same
rj$iter["t.ij"] <- rj$iter["t.ij"] + 1

if(rj$config$function.select | rj$phase[i - 1] == 1){
  rj$iter[c("sigma", "mu.i", "phi")] <- 
    rj$iter[c("sigma", "mu.i", "phi")] + 1
}

if(rj$config$function.select | rj$phase[i - 1] == 2){
  rj$iter[c("tau", "omega", "psi", "mu.ij", "k.ij", "psi.i", "pi.ij")] <- 
    rj$iter[c("tau", "omega", "psi", "mu.ij", "k.ij", "psi.i", "pi.ij")] + 1
}

if(rj$dat$covariates$n > 0){
  
  rj$pi.ij[i, ] <- 
    pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id] + 
            cov_effects(rj.obj = rj))
  
} else {
  
  rj$pi.ij[i, ] <- pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id])
  
}

mcmat$pi.ij[i,] <- rj$pi.ij[i,]

results <- tibble::tibble(param = paste0(rep(c("alpha", "nu1", "nu2", "tau1", "tau2", "omega", "mu.ij", "psi.i", "k.ij"), each = 2), rep(c(".proposed", ".current"), times = 9)), espresso = NA, hard_coded = NA)

#'------------------------------
# // alpha ----
#'------------------------------
set.seed(134)
proposed.alpha <- proposal_mh(rj.obj = rj, param.name = "alpha")
results[1, 2] <- likelihood(biphasic = TRUE,
                                param.name = "alpha",
                                rj.obj = rj,
                                model = rj$mlist[[rj$current.model]],
                                values = proposed.alpha,
                                included.cov = NULL,
                                RJ = FALSE,
                                lprod = TRUE)
results[2, 2] <- likelihood(biphasic = TRUE,
                               param.name = "alpha",
                               rj.obj = rj,
                               model = rj$mlist[[rj$current.model]],
                               values = NULL,
                               included.cov = NULL,
                               RJ = FALSE,
                               lprod = TRUE)
proposed.alpha <- proposed.alpha[[1]]
results[1, 3] <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                              location = mcmat$nu[i - 1, 1], 
                              scale = mcmat$tau[i - 1, 1], 
                              L = L, U = proposed.alpha, log = TRUE)) + 
  sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ], 
             location = mcmat$nu[i - 1, 2], 
             scale = mcmat$tau[i - 1, 2], 
             L = proposed.alpha, U = U, log = TRUE))
results[2, 3] <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                             location = mcmat$nu[i - 1, 1], 
                             scale = mcmat$tau[i - 1, 1], 
                             L = L, U = mcmat$alpha[i - 1], log = TRUE)) + 
  sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ], 
             location = mcmat$nu[i - 1, 2], 
             scale = mcmat$tau[i - 1, 2], 
             L = mcmat$alpha[i - 1], U = U, log = TRUE))

mcmat$alpha[i] <- mcmat$alpha[i - 1]

proposed.nu <- proposal_mh(rj.obj = rj, param.name = "nu")
results[3, 2] <- likelihood(biphasic = TRUE,
                                param.name = "nu1",
                                rj.obj = rj,
                                model = rj$mlist[[rj$current.model]],
                                values = proposed.nu,
                                included.cov = NULL,
                                RJ = FALSE,
                                lprod = TRUE)
  
results[4, 2] <-  likelihood(biphasic = TRUE,
                               param.name = "nu1",
                               rj.obj = rj,
                               model = rj$mlist[[rj$current.model]],
                               values = NULL,
                               included.cov = NULL,
                               RJ = FALSE,
                               lprod = TRUE)
results[3, 3] <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                              location = proposed.nu[[1]][1,],
                              scale = mcmat$tau[i - 1, 1],
                              L = L, U = mcmat$alpha[i], log = TRUE))
results[4, 3]  <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                             location = mcmat$nu[i - 1, 1],
                             scale = mcmat$tau[i - 1, 1],
                             L = L, U = mcmat$alpha[i], log = TRUE))
mcmat$nu[i, 1] <- mcmat$nu[i - 1, 1]

 results[5,2] <- likelihood(biphasic = TRUE,
                                param.name = "nu2",
                                rj.obj = rj,
                                model = rj$mlist[[rj$current.model]],
                                values = proposed.nu,
                                included.cov = NULL,
                                RJ = FALSE,
                                lprod = TRUE)
  
  results[6,2] <- likelihood(biphasic = TRUE,
                               param.name = "nu2",
                               rj.obj = rj,
                               model = rj$mlist[[rj$current.model]],
                               values = NULL,
                               included.cov = NULL,
                               RJ = FALSE,
                               lprod = TRUE)
  results[5,3] <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                location = proposed.nu[[1]][2,],
                                scale = mcmat$tau[i - 1, 2],
                                L = mcmat$alpha[i], U = U, log = TRUE))
  results[6,3] <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                               location = mcmat$nu[i - 1, 2],
                               scale = mcmat$tau[i - 1, 2],
                               L = mcmat$alpha[i], U = U, log = TRUE))
  mcmat$nu[i, 2] <- mcmat$nu[i - 1, 2]
  
  proposed.tau <- proposal_mh(rj.obj = rj, param.name = "tau")
  results[7,2] <- likelihood(biphasic = TRUE,
                                param.name = "tau1",
                                rj.obj = rj,
                                model = rj$mlist[[rj$current.model]],
                                values = proposed.tau,
                                included.cov = NULL,
                                RJ = FALSE,
                                lprod = TRUE)
  
  results[8,2]<- likelihood(biphasic = TRUE,
                               param.name = "tau1",
                               rj.obj = rj,
                               model = rj$mlist[[rj$current.model]],
                               values = NULL,
                               included.cov = NULL,
                               RJ = FALSE,
                               lprod = TRUE)
  
  results[7,3] <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                                location = mcmat$nu[i, 1],
                                scale = proposed.tau[[1]][1],
                                L = L, U = mcmat$alpha[i], log = TRUE))
  results[8,3] <- sum(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ],
                               location = mcmat$nu[i, 1],
                               scale = mcmat$tau[i - 1, 1],
                               L = L, U = mcmat$alpha[i], log = TRUE))
  mcmat$tau[i, 1] <- mcmat$tau[i - 1, 1]
  
  results[9,2] <- likelihood(biphasic = TRUE,
                                param.name = "tau2",
                                rj.obj = rj,
                                model = rj$mlist[[rj$current.model]],
                                values = proposed.tau,
                                included.cov = NULL,
                                RJ = FALSE,
                                lprod = TRUE)
  
  results[10,2] <- likelihood(biphasic = TRUE,
                               param.name = "tau2",
                               rj.obj = rj,
                               model = rj$mlist[[rj$current.model]],
                               values = NULL,
                               included.cov = NULL,
                               RJ = FALSE,
                               lprod = TRUE)
  
  results[9, 3] <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                location = mcmat$nu[i, 2],
                                scale = proposed.tau[[1]][2],
                                L = mcmat$alpha[i], U = U, log = TRUE))
  results[10,3] <- sum(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                               location = mcmat$nu[i, 2],
                               scale = mcmat$tau[i - 1, 2],
                               L = mcmat$alpha[i], U = U, log = TRUE))
  mcmat$tau[i, 2] <- mcmat$tau[i - 1, 2]
  
  proposed.omega <- proposal_mh(rj.obj = rj, param.name = "omega")
  results[11,2] <- likelihood(biphasic = TRUE,
                                  param.name = "omega",
                                  rj.obj = rj,
                                  model = rj$mlist[[rj$current.model]],
                                  values = proposed.omega,
                                  included.cov = NULL,
                                  RJ = FALSE,
                                  lprod = TRUE)
    
    results[12,2] <- likelihood(biphasic = TRUE,
                                 param.name = "omega",
                                 rj.obj = rj,
                                 model = rj$mlist[[rj$current.model]],
                                 values = NULL,
                                 included.cov = NULL,
                                 RJ = FALSE,
                                 lprod = TRUE)
    results[11,3] <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = mcmat$psi[i - 1],
                                 sd = proposed.omega[[1]], log = TRUE))
    results[12,3] <- sum(dnorm(x = mcmat$psi.i[i - 1, ], mean = mcmat$psi[i - 1],
                                sd = mcmat$omega[i - 1], log = TRUE))
    
    mcmat$omega[i] <- mcmat$omega[i - 1]
    
    # Gibbs sampler as we have a Normal prior and a Normal likelihood
    tmp.calc <- rj$dat$whales$n / rj$omega[rj$iter["omega"]] ^ 2
    sigma2.post <- 1 / (rj$config$psi.gibbs[1] + tmp.calc) 
    mu.post <- sigma2.post * (rj$config$psi.gibbs[2] + 
                                mean(rj$psi.i[rj$iter["psi.i"], ]) * tmp.calc)
    rj$psi[i] <- mcmat$psi[i] <- rnorm(1, mean = mu.post, sd = sqrt(sigma2.post))
    
    proposed.mu.ij <- proposal_mh(rj.obj = rj, param.name = "mu.ij")
    
    results[13,2] <- sum(likelihood(biphasic = TRUE,
                                  param.name = "mu.ij",
                                  rj.obj = rj,
                                  model = rj$mlist[[rj$current.model]],
                                  values = proposed.mu.ij,
                                  included.cov = NULL,
                                  RJ = FALSE,
                                  lprod = FALSE))
    
    results[14,2] <- sum(likelihood(biphasic = TRUE,
                                 param.name = "mu.ij",
                                 rj.obj = rj,
                                 model = rj$mlist[[rj$current.model]],
                                 values = NULL,
                                 included.cov = NULL,
                                 RJ = FALSE,
                                 lprod = FALSE))
    
    # Set lower bounds for each of the two mixture components
    L.1 <- rep(L, n.trials)
    L.2 <- rep(mcmat$alpha[i], n.trials)
    # When right-censoring occurs, the lower bounds of the proposal distribution must be adjusted
    ind.1 <- which(is.censored == 1)[which(mcmat$k.ij[i - 1, which(is.censored == 1)] == 1)]
    ind.2 <- which(is.censored == 1)[which(mcmat$k.ij[i - 1, which(is.censored == 1)] == 2)]
    L.1[ind.1] <- purrr::map_dbl(.x = max.dose[ind.1], .f = ~max(.x, L))
    L.2[ind.2] <- purrr::map_dbl(.x = max.dose[ind.2], .f = ~max(.x, mcmat$alpha[i]))
    ##LT: Not sure we need to save these?  Perhaps drop L1 and L2 from the list of things saved?  
    ##LT: For one thing L1 is the same for each iteration!  (Indeed the above calculations for L1 could
    ##LT: be taken outside the sampling loop to help speed things up)  
    mcmat$L1[i, ] <- L.1
    mcmat$L2[i, ] <- L.2
    
    mu_ij_1 <- proposed.mu.ij[[1]][,1]
    mu_ij_2 <- proposed.mu.ij[[1]][,2]
    
    proposed.t.ij <- (2 - mcmat$k.ij[i - 1, ]) * mu_ij_1 + ((1 - (2 - mcmat$k.ij[i - 1, ])) * mu_ij_2)
    
    ##LT: There used to be code here creating a list for loglik.proposed etc -- but
    ##LT: since this was the same for every iteration of the i-loop I moved it outside the
    ##LT: loop; I also turned it from a list into a matrix as that's more efficient
    # The likelihood is evaluated using the standard bounds
    mu.loglik.proposed <- mu.loglik.current <- matrix(nrow = length(proposed.t.ij), ncol = 2)
    mu.loglik.proposed[, 1] <- unname(dtnorm(x = mu_ij_1, 
                                             location = mcmat$nu[i, 1], 
                                             scale = mcmat$tau[i, 1],
                                             L = L, 
                                             U = mcmat$alpha[i], 
                                             log = TRUE))
    mu.loglik.proposed[, 2] <- unname(dtnorm(x = mu_ij_2,
                                             location = mcmat$nu[i, 2],
                                             scale = mcmat$tau[i, 2],
                                             L = mcmat$alpha[i],
                                             U = U, 
                                             log = TRUE))
    ##LT: 1. This was not dealing with right-censored observations correctly.
    ##LT: They have a likelihood of 1 (so a zero value on the log scale) if the t.ij is more than the right censor point and zero (so a -Inf on log scale) otherwise.  
    ##LT: 2. the t.ij likelihood was being applied to both the mu1 and mu2 lists, which means the
    ##LT: observation likelihood was being applied twice -- I've modified it so that it just
    ##LT: applies to the relevant mixture component, depending on what k is pointing to
    loglik.obs <- dnorm(y_ij, mean = proposed.t.ij, sd = obs.sd, log = TRUE)
    loglik.obs[is.na(loglik.obs)] <- 0
    for(a in 1:n.trials) {
      mu.loglik.proposed[a, mcmat$k.ij[i - 1, a]] <- mu.loglik.proposed[a, mcmat$k.ij[i - 1, a]] + loglik.obs[a]
    }

    results[13, 3] <- sum(mu.loglik.proposed)  
    
    mu.loglik.current[, 1] <- unname(dtnorm(x = mcmat$mu.ij$`1`[i - 1, ], 
                                            location = mcmat$nu[i, 1], 
                                            scale = mcmat$tau[i, 1],
                                            L = L, 
                                            U = mcmat$alpha[i], 
                                            log = TRUE))
    mu.loglik.current[, 2] <- unname(dtnorm(x = mcmat$mu.ij$`2`[i - 1, ],
                                            location = mcmat$nu[i, 2],
                                            scale = mcmat$tau[i, 2],
                                            L = mcmat$alpha[i],
                                            U = U, 
                                            log = TRUE))
    loglik.obs <- dnorm(y_ij, mean = mcmat$t.ij[i-1, ], sd = obs.sd, log = TRUE)
    loglik.obs[is.na(loglik.obs)] <- 0
    for(a in 1:n.trials) {
      mu.loglik.current[a, mcmat$k.ij[i - 1, a]] <- mu.loglik.current[a, mcmat$k.ij[i - 1, a]] + loglik.obs[a]
    }
    
    results[14, 3] <- sum(mu.loglik.current)  
    
    mcmat$mu.ij[[1]][i, ] <- mcmat$mu.ij[[1]][i - 1, ]
    mcmat$mu.ij[[2]][i, ] <- mcmat$mu.ij[[2]][i - 1, ]
    mcmat$t.ij[i, ] <- ((2 - mcmat$k.ij[i - 1, ]) * mcmat$mu.ij$`1`[i, ]) +
      ((1 - (2 - mcmat$k.ij[i - 1, ])) * mcmat$mu.ij$`2`[i, ])
    
    proposed.psi.i <- proposal_mh(rj.obj = rj, param.name = "psi.i")
    
    results[15,2] <- sum(likelihood(biphasic = TRUE,
                                  param.name = "psi.i",
                                  rj.obj = rj,
                                  model = rj$mlist[[rj$current.model]],
                                  values = proposed.psi.i,
                                  included.cov = NULL,
                                  RJ = FALSE,
                                  lprod = FALSE))
    
    results[16,2] <- sum(likelihood(biphasic = TRUE,
                                 param.name = "psi.i", 
                                 rj.obj = rj,
                                 model = rj$mlist[[rj$current.model]],
                                 values = NULL,
                                 included.cov = NULL,
                                 RJ = FALSE,
                                 lprod = FALSE))

    proposed.psi.ij <- Reduce("+", append(list(proposed.psi.i[[1]][whale.id]),
                                          lapply(X = covariate.names, 
                                                 FUN = function(k){
                                                   apply(t(t(dummy.cov[[k]]) * 
                                                             qnorm(pnorm(mcmat[[k]][i - 1, ], 
                                                                         mean = priors[[k]][1],
                                                                         sd = priors[[k]][2]))), 1, sum)})))
    current.psi.ij <- Reduce("+", append(list(mcmat$psi.i[i - 1, whale.id]),
                                         lapply(X = covariate.names, 
                                                FUN = function(k){
                                                  apply(t(t(dummy.cov[[k]]) * 
                                                            qnorm(pnorm(mcmat[[k]][i - 1, ], 
                                                                        mean = priors[[k]][1],
                                                                        sd = priors[[k]][2]))), 1, sum)})))
    results[15,3] <- sum(dnorm(x = proposed.psi.i[[1]],
                             mean = mcmat$psi[i],
                             sd = mcmat$omega[i],
                             log = TRUE) +
      purrr::map_dbl(.x = seq_len(n.whales),
                     .f = ~ sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                        prob = pnorm(q = proposed.psi.ij), log = TRUE)[whale.id == .x])))
    
    results[16,3] <- sum(dnorm(x = mcmat$psi.i[i - 1, ],
                            mean = mcmat$psi[i],
                            sd = mcmat$omega[i],
                            log = TRUE) +
      purrr::map_dbl(.x = seq_len(n.whales),
                     .f = ~ sum(d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1,
                                        prob = pnorm(q = current.psi.ij), log = TRUE)[whale.id == .x])))
    
    mcmat$psi.i[i, ] <- mcmat$psi.i[i - 1, ]
    
    for(a in rj$dat$covariates$names){
    mcmat[[a]][i, ] <- mcmat[[a]][i - 1, ]
    }
    
    proposed.k.ij <- proposal_mh(rj.obj = rj, param.name = "k.ij")
    
    results[17, 2] <- sum(likelihood(biphasic = TRUE,
                                  param.name = "k.ij",
                                  rj.obj = rj,
                                  model = rj$mlist[[rj$current.model]],
                                  values = proposed.k.ij,
                                  included.cov = NULL,
                                  RJ = FALSE,
                                  lprod = FALSE))
    
    results[18,2] <- sum(likelihood(biphasic = TRUE,
                                 param.name = "k.ij",
                                 rj.obj = rj,
                                 model = rj$mlist[[rj$current.model]],
                                 values = NULL,
                                 included.cov = NULL,
                                 RJ = FALSE,
                                 lprod = FALSE))
    
    proposed.t.ij <- purrr::map_dbl(.x = seq_len(n.trials),
                                    .f = ~ mcmat$mu.ij[[proposed.k.ij[[1]][.x]]][i, .x])
    
    ##LT: Changed to use y_ij rather than mamat$y_ij  
    loglik.proposed <- d_binom(x = 2 - proposed.k.ij[[1]], size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
    loglik.obs <- dnorm(x = y_ij, mean = proposed.t.ij, sd = obs.sd, log = TRUE)
    loglik.obs[is.na(loglik.obs)] <- 0
    results[17,3] <- sum(loglik.proposed + loglik.obs)
    
    loglik.current <- d_binom(x = 2 - mcmat$k.ij[i - 1, ], size = 1, prob = mcmat$pi.ij[i, ], log = TRUE)
    loglik.obs <- dnorm(x = y_ij, mean = mcmat$t.ij[i, ], sd = obs.sd, log = TRUE)
    loglik.obs[is.na(loglik.obs)] <- 0
    results[18,3] <- sum(loglik.current + loglik.obs)
    
    
    
    
    
    results$diff <- round(results$hard_coded - results$espresso, 10)
    print(results)
    
    