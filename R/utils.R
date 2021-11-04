# Distributions -----------------------------------------------------------

#' The Truncated Normal Distribution
#'
#' Probability density for a truncated Normal distribution with mean equal to \code{location} and standard deviation equal to \code{scale}.
#'
#' @param x Vector of quantiles.
#' @param location Vector of means.
#' @param scale Vector of standard deviations.
#' @param log Logical. If \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param L Lower limit of the distribution.
#' @param U Upper limit of the distribution.

dtnorm <- function(x, location = 0, scale = 1, log = FALSE, L = -Inf, U = Inf) {

  d <- dnorm(x, location, scale, log = TRUE)
  denom <- log(pnorm(U, location, scale) - pnorm(L, location, scale))
  d <- d - denom
  d[(x < L) | (x > U)] <- -100000 # When input quantile is outside bounds
  d[is.infinite(d)] <- -100000 # When input location is outside bounds
  if(!log) d <- exp(d)
  return(d)
}

#' The Truncated Normal Distribution
#'
#' Generate random deviates from a truncated normal distribution with mean equal to \code{location} and standard deviation equal to \code{scale}.
#'
#' @param n Number of observations.
#' @param location Vector of means.
#' @param scale Vector of standard deviations.
#' @param log Logical. If \code{TRUE}, probabilities p are given as log(p).
#' @param L Lower limit of the distribution.
#' @param U Upper limit of the distribution.

rtnorm <- function(n, location, scale, L, U){
  
  location + scale * qnorm(pnorm(L, location, scale) + runif(n)*(pnorm(U, location, scale) - pnorm(L, location, scale)))
  
}

d_binom <- function(x, size, prob, log){
  d <- dbinom(x = x, size = size, prob = prob, log = log)
  d[is.infinite(d)] <- -100000 # To avoid Inf that cause numerical issues
  return(d)
}

# Likelihood --------------------------------------------------------------

#' Likelihood
#'
#' Calculate the log-likelihood of a monophasic dose-response model, as required by the rjMCMC sampler.
#' 
#' @param biphasic Logical. Indicates the type of model for which likelihoods should be calculated.
#' @param param.name Paramter name. Ignored if \code{values} is specified.
#' @param rj.obj rjMCMC object.
#' @param model Integer vector determining species groupings for the model of interest.
#' @param values List of proposed values, if available.
#' @param included.cov Boolean vector indicating which contextual covariates are included in the current iteration of the rjMCMC sampler.
#' @param RJ Logical. If \code{TRUE}, returns the likelihoods associated with between-model jumps. If \code{FALSE}, returns the likelihoods associated with the update step of the Metropolis sampler.
#' @param lprod If \code{TRUE}, returns the product of individual likelihoods. 

likelihood <- function(biphasic = FALSE,
                       param.name = NULL, 
                       rj.obj,
                       model, 
                       values = NULL, 
                       included.cov = NULL, 
                       RJ = FALSE,
                       lprod = TRUE){
  
  # Log-likelihoods
  loglikelihood <- 0
  
  if(!biphasic){
    
    # MONOPHASIC
    
    if(any(param.name == "t.ij") | RJ){
      LL.1 <- dnorm(x = rj.obj$dat$obs$y_ij, 
                    mean = if(any(param.name == "t.ij") & !is.null(values)) values$t.ij else
                      rj.obj$t.ij[rj.obj$iter["t.ij"],], 
                    sd = rj.obj$dat$obs$sd, 
                    log = TRUE)
      # Censored data have a likelihood of 1 (0 on log scale) if t.ij > Rc or < Lc
      # which we've made sure is the case by only proposing values of t.ij that meet these criteria
      LL.1[is.na(LL.1)] <- 0 
      if(lprod) LL.1 <- sum(LL.1)
      loglikelihood <- loglikelihood + LL.1
    }
    
    if(any(param.name %in% c("t.ij", "sigma", "mu.i", rj.obj$dat$covariates$names)) | RJ){
      
      if(rj.obj$dat$covariates$n > 0){
        
        if(is.null(included.cov)) included.cov <- rj.obj$include.covariates[rj.obj$iter["covariates"], ]
        
        if (sum(included.cov) == 0) {
          covnames <- "nil"
          cov.effects <- 0
        } else {
          included.cov <- included.cov[included.cov == 1]
          covvalues <- do.call(c, purrr::map(
            .x = names(included.cov),
            .f = ~ {
              if (!is.null(values) & .x %in% names(values)) values[[.x]] else rj.obj[[.x]][rj.obj$iter[.x], ]
            }
          ))
          cov.effects <- colSums(t(do.call(cbind, rj.obj$dat$covariates$dummy[names(included.cov)])) * covvalues)
        }
      } else {
        cov.effects <- 0
      }
      
      LL.2 <- dtnorm(
        x = if(any(param.name == "t.ij") & !is.null(values)) values$t.ij else rj.obj$t.ij[rj.obj$iter["t.ij"],],
        location = if(any(param.name == "mu.i") & !is.null(values)) 
          values$mu.i[rj.obj$dat$whales$id] + cov.effects else 
          rj.obj$mu.i[rj.obj$iter["mu.i"], rj.obj$dat$whales$id] + cov.effects,
        scale = if(any(param.name == "sigma") & !is.null(values)) values$sigma else 
          rj.obj$sigma[rj.obj$iter["sigma"]],
        L = rj.obj$dat$param$dose.range[1],
        U = rj.obj$dat$param$dose.range[2],
        log = TRUE)
      
      if (any(param.name == "mu.i") & !RJ) {
        LL.2 <- purrr::map_dbl(.x = seq_len(rj.obj$dat$whales$n), 
                               .f = ~ sum(LL.2[rj.obj$dat$whales$id == .x]))
      } 
      if (lprod) LL.2 <- sum(LL.2)
      loglikelihood <- loglikelihood + LL.2
    }
    
    if(any(param.name %in% c("mu.i", "mu", "phi")) | RJ){
      
      LL.3 <- dtnorm(
        x = if(any(param.name == "mu.i") & !is.null(values)) values$mu.i else rj.obj$mu.i[rj.obj$iter["mu.i"],],
        location = if(any(param.name == "mu") & !is.null(values)) values$mu[rj.obj$dat$species$id] else rj.obj$mu[rj.obj$iter["mu"], rj.obj$dat$species$id],
        scale = if(any(param.name == "phi") & !is.null(values)) values$phi else rj.obj$phi[rj.obj$iter["phi"]],
        L = rj.obj$dat$param$dose.range[1],
        U = rj.obj$dat$param$dose.range[2],
        log = TRUE)
      
      if (lprod) LL.3 <- sum(LL.3)
      loglikelihood <- loglikelihood + LL.3
    }
    
    if(lprod) return(sum(loglikelihood, na.rm = TRUE)) else return(unname(loglikelihood)) 
    
  } else {
    
    # BIPHASIC
    
    if(any(c("mu.ij", "k.ij") %in% param.name) | RJ){
      LL.1 <- dnorm(x = rj.obj$dat$obs$y_ij, 
                    mean = if("mu.ij" %in% param.name & !is.null(values) & !RJ) 
                      values$mu.ij[cbind(seq_len(nrow(values$mu.ij)), rj.obj$k.ij[rj.obj$iter["k.ij"], ])] else if("k.ij" %in% param.name & !is.null(values) & !RJ) rj.obj$mu.ij[rj.obj$iter["mu.ij"], ,][cbind(seq_along(values$k.ij), values$k.ij)] else rj.obj$t.ij[rj.obj$iter["t.ij"],], 
                    sd = rj.obj$dat$obs$sd, 
                    log = TRUE)
      
      LL.1[is.na(LL.1)] <- 0 
      if(lprod & !RJ) LL.1 <- sum(LL.1)
      loglikelihood <- loglikelihood + LL.1
    }
    
    if(any(c("mu.ij", "nu1", "tau1", "alpha") %in% param.name) | RJ){
      
      LL.2a <- 
        dtnorm(x = if("mu.ij" %in% param.name & !is.null(values)) 
          values$mu.ij[, 1] else rj.obj$mu.ij[rj.obj$iter["mu.ij"], , 1],
               location = if(any(c("nu1", "nu") %in% param.name) & !is.null(values)) values$nu[1, rj.obj$dat$species$trials] else rj.obj$nu[rj.obj$iter["nu"], rj.obj$dat$species$trials, 1],
               scale = if("tau1" %in% param.name & !is.null(values)) values$tau[1] else rj.obj$tau[rj.obj$iter["tau"], 1],
               L = rj.obj$dat$param$dose.range[1], 
               U = if("alpha" %in% param.name & !is.null(values)) values$alpha[rj.obj$dat$species$trials] else rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials],
               log = TRUE)
      
      if(lprod & !RJ) LL.2a <- sum(LL.2a)
      if(!"mu.ij" %in% param.name & !RJ) loglikelihood <- loglikelihood + LL.2a
      
    }
    
    if(any(c("mu.ij", "nu2", "tau2", "alpha") %in% param.name) | RJ){
      
      LL.2b <- 
        dtnorm(x = if("mu.ij" %in% param.name & !is.null(values)) 
          values$mu.ij[, 2] else rj.obj$mu.ij[rj.obj$iter["mu.ij"], , 2],
               location = if(any(c("nu2", "nu") %in% param.name) & !is.null(values)) values$nu[2, rj.obj$dat$species$trials] else rj.obj$nu[rj.obj$iter["nu"], rj.obj$dat$species$trials, 2],
               scale = if("tau2" %in% param.name & !is.null(values)) values$tau[2] else rj.obj$tau[rj.obj$iter["tau"], 2],
               L = if("alpha" %in% param.name & !is.null(values)) values$alpha[rj.obj$dat$species$trials] else rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials],
               U = rj.obj$dat$param$dose.range[2],
               log = TRUE)
      
      if(lprod & !RJ) LL.2b <- sum(LL.2b)
      
      if(!"mu.ij" %in% param.name & !RJ){
        
        loglikelihood <- loglikelihood + LL.2b
        
      } else {
        
        mu.ij.loglik <- cbind(LL.2a, LL.2b)
        
       # if(RJ | lprod){
       #   
       #   mu.ij.loglik <-
       #     mu.ij.loglik[cbind(seq_len(nrow(mu.ij.loglik)), rj.obj$k.ij[rj.obj$iter["k.ij"], ])] + loglikelihood
       # 
       #   loglikelihood <- sum(mu.ij.loglik)
       #   
       #   } else {
         
         mu.ij.loglik[cbind(seq_len(nrow(mu.ij.loglik)), rj.obj$k.ij[rj.obj$iter["k.ij"], ])] <-
           mu.ij.loglik[cbind(seq_len(nrow(mu.ij.loglik)), rj.obj$k.ij[rj.obj$iter["k.ij"], ])] + loglikelihood
         
        if(RJ | lprod) loglikelihood <- sum(mu.ij.loglik) else loglikelihood <- mu.ij.loglik
         
         # }
      }
    }
    
    # Note: psi updated via Gibbs sampling so not included here
    if(any(c("psi.i", "k.ij", rj.obj$dat$covariates$names) %in% param.name) | RJ){
      
      if(rj.obj$dat$covariates$n > 0){
        
      if(is.null(included.cov)) included.cov <- rj.obj$include.covariates[rj.obj$iter["covariates"], ]
      
      if (sum(included.cov) == 0) {
        covnames <- "nil"
        cov.effects <- 0
      } else {
        included.cov <- included.cov[included.cov == 1]
        covvalues <- do.call(c, purrr::map(
          .x = names(included.cov),
          .f = ~ {
            if (!is.null(values) & (is.list(values) & .x %in% names(values))) {
              qnorm(pnorm(q = values[[.x]],
                          mean = rj.obj$config$priors[.x, 1],
                          sd = rj.obj$config$priors[.x, 2]))
            } else { 
              qnorm(pnorm(q = rj.obj[[.x]][rj.obj$iter[.x], ], 
                          mean = rj.obj$config$priors[.x, 1], 
                          sd = rj.obj$config$priors[.x, 2]))
            }
          }
        ))
        cov.effects <- colSums(t(do.call(cbind, rj.obj$dat$covariates$dummy[names(included.cov)])) * covvalues)
      }
      
      } else {
        cov.effects <- 0
      }
      
      if("psi.i" %in% param.name & !is.null(values)) psi.ij <- values$psi.i[rj.obj$dat$whales$id] + cov.effects else psi.ij <- rj.obj$psi.i[rj.obj$iter["psi.i"], rj.obj$dat$whales$id] + cov.effects
      
      pi.ij <- pnorm(psi.ij)
      
    }
    
    if(any(c("omega", "psi.i") %in% param.name) | RJ){
      LL.3 <- 
      dnorm(x = if("psi.i" %in% param.name & !is.null(values)) values$psi.i else 
        rj.obj$psi.i[rj.obj$iter["psi.i"], ],
              mean = if("psi" %in% param.name & !is.null(values)) values$psi else 
                rj.obj$psi[rj.obj$iter["psi"]],
              sd = if("omega" %in% param.name & !is.null(values)) values$omega else
                rj.obj$omega[rj.obj$iter["omega"]],
              log = TRUE)
      if(lprod) LL.3 <- sum(LL.3)
      loglikelihood <- loglikelihood + LL.3
    }
    
    if(any(c("k.ij", "covariates", "psi.i", rj.obj$dat$covariates$names) %in% param.name) | RJ){
      
      LL.4 <- d_binom(x = if("k.ij" %in% param.name & !is.null(values)) 2 - values$k.ij else 
        2 - rj.obj$k.ij[rj.obj$iter["k.ij"], ], 
        size = 1, 
        prob = pi.ij, 
        log = TRUE)
      
      LL.5 <- purrr::map_dbl(.x = seq_len(rj.obj$dat$whales$n),
                             .f = ~ sum(LL.4[rj.obj$dat$whales$id == .x]))
      
      if(lprod) LL.4 <- sum(LL.4)
      if(any(c("k.ij", rj.obj$dat$covariates$names) %in% param.name) | RJ) 
        loglikelihood <- loglikelihood + LL.4
      
    }
    
    if("psi.i" %in% param.name & !RJ){
      if(lprod) LL.5 <- sum(LL.5)
      loglikelihood <- loglikelihood + LL.5 
    }
      
    #   LL.4a <- d_binom(x = if("k.ij" %in% param.name & !is.null(values)) 2 - values$k.ij else 
    #     2 - rj.obj$k.ij[rj.obj$iter["k.ij"], ], 
    #                    size = 1, 
    #                    prob = pi.ij, 
    #                    log = TRUE)
    #   
    #   if(lprod) LL.4a <- sum(LL.4a)
    #   if(any(c("k.ij", rj.obj$dat$covariates$names) %in% param.name)) 
    #     loglikelihood <- loglikelihood + LL.4a
    #   
    # }
    # 
    # if("psi.i" %in% param.name | RJ){
    #   LL.4b <- purrr::map_dbl(.x = seq_len(rj.obj$dat$whales$n),
    #                           .f = ~ sum(LL.4a[rj.obj$dat$whales$id == .x]))
    #   if(lprod) LL.4b <- sum(LL.4b)
    #   loglikelihood <- loglikelihood + LL.4b  
    # }
    
    if(lprod) return(sum(loglikelihood, na.rm = TRUE)) else return(unname(loglikelihood))
    
  }  
}

# Priors ----------------------------------------------------------------

#' Uniform prior
#'
#' Probability density for a Uniform distribution.
#' 
#' @param rj.obj rjMCMC object.
#' @param param.name Parameter name.
#' @param param Parameter values.

uniform_prior <- function(rj.obj,
                          param.name,
                          param = NULL) {
  
  if("mu" %in% param.name) {
  loglik.unif <- dunif(x = unique(param$out$mu), 
                       min = rj.obj$config$priors[param.name, 1],
                       max = rj.obj$config$priors[param.name, 2],
                       log = TRUE)
  
  } else {
    loglik.unif <- c(dunif(x = unique(param$out$nu[1, ]), 
                         min = rj.obj$config$priors["nu", 1],
                         max = unique(param$out$alpha),
                         log = TRUE),
                     dunif(x = unique(param$out$nu[2, ]), 
                           min = unique(param$out$alpha),
                           max = rj.obj$config$priors["nu", 2],
                           log = TRUE),
                     dunif(x = unique(param$out$alpha), 
                           min = rj.obj$config$priors["alpha", 1],
                           max = rj.obj$config$priors["alpha", 2],
                           log = TRUE))
  }
  
  if (any(abs(loglik.unif) == Inf)) loglik.unif <- -100000
  return(sum(loglik.unif))
}

ff_prior <- function(rj.obj, param = NULL, phase) {
  
  if(phase == 1) {

    if(is.null(param)) param <- list(mu = rj.obj$mu[rj.obj$iter["mu"], ],
                                     sigma = rj.obj$sigma[rj.obj$iter["sigma"]],
                                     phi = rj.obj$phi[rj.obj$iter["phi"]])
    
    loglik <- dunif(x = c(unique(param$mu), param$sigma, param$phi),
                         min = rj.obj$config$priors[c(rep("mu", length(unique(param$mu))), "sigma", "phi"), 1],
                         max = rj.obj$config$priors[c(rep("mu", length(unique(param$mu))), "sigma", "phi"), 2],
                         log = TRUE)
    
  } else if(phase == 2){
    
    if(is.null(param)){
      param <- list(alpha = unique(rj.obj$alpha[rj.obj$iter["alpha"], ]),
                    nu1 = unique(rj.obj$nu[rj.obj$iter["nu"], , 1]),
                    nu2 = unique(rj.obj$nu[rj.obj$iter["nu"], , 2]),
                    tau = rj.obj$tau[rj.obj$iter["tau"], ],
                    omega = rj.obj$omega[rj.obj$iter["omega"]],
                    psi = rj.obj$psi[rj.obj$iter["psi"]])
    } else {
      param$nu1 <- unique(param$nu[1, ])
      param$nu2 <- unique(param$nu[2, ])
      param$alpha <- unique(param$alpha)
    }

    loglik <- dunif(x = c(param$alpha, param$nu1, param$nu2, param$tau, param$omega),
          min = c(rj.obj$config$priors[rep("alpha", length(param$alpha)), 1],
                  rj.obj$config$priors[rep("nu", length(param$nu1)), 1],
                  param$alpha,
                  rj.obj$config$priors[c(rep("tau", 2), "omega"), 1]),
          max = c(rj.obj$config$priors[rep("alpha", length(param$alpha)), 2],
                  param$alpha,
                  rj.obj$config$priors[rep("nu", length(param$nu2)), 2],
                  rj.obj$config$priors[c(rep("tau", 2), "omega"), 2]),
          log = TRUE) + 
      dnorm(param$psi, mean = rj.obj$config$priors["psi", 1], sd = rj.obj$config$priors["psi", 2], log = TRUE)
    
  }
  
  if (any(abs(loglik) == Inf)) loglik <- -100000
  return(sum(loglik))
}

#' Normal prior on covariates
#'
#' Probability density for a Normal distribution.
#'
#' @param rj.obj rjMCMC object.
#' @param param.name Parameter name
#' @param param Parameter values

normal_prior <- function(rj.obj,
                         param.name, 
                         param = NULL) {
  
  loglik.norm <- dnorm(x = unique(param), 
                       mean = rj.obj$config$priors[param.name, 1], 
                       sd = rj.obj$config$priors[param.name, 2], 
                       log = TRUE)
  
  if (any(abs(loglik.norm) == Inf)) loglik.norm <- -100000
  return(sum(loglik.norm))
}
# Proposals ----------------------------------------------------------------

#' Proposal for between-model jumps
#'
#' Propose a new model.
#' 
#' @param move.type Type of between-model jump. (0): split/merge, (1): Bootstrap, (2): Cluster, (3): Random
propose_jump <- function(rj.obj, move.type, phase) {
  
  # Perform function checks
  if (!move.type %in% 0:3) stop("Unrecognised move type.")
  
  if (move.type > 0) {
    
    # Bootstrap
    if (move.type == 1) {
      
      # Rescale the probability values so that none are zero
      # if(any(rj.obj$config$clust[[1]]$p == 0)){
      #   
      #   zero.ind <- which(rj.obj$config$clust[[1]]$p == 0)
      #   positive.ind <- which(rj.obj$config$clust[[1]]$p > 0)
      #   
      #   # Use half of the minimum value
      #   scale.value <- min(rj.obj$config$clust[[1]]$p[positive.ind]) / 2
      #   
      #   rj.obj$config$clust[[1]]$p_scale[zero.ind] <- scale.value / length(zero.ind)
      #   
      #   rj.obj$config$clust[[1]]$p_scale[positive.ind] <- 
      #     rj.obj$config$clust[[1]]$p[positive.ind] - (scale.value * rj.obj$config$clust[[1]]$p[positive.ind])
      # }
      
      # Rescale probabilities so they sum to 1
      p1 <- rj.obj$config$clust[[1]]
      if(nrow(p1) >  1) p1 <- p1 %>% dplyr::filter(!model == rj.obj$current.model)
      p1 <- p1 %>% dplyr::mutate(p_scale = rescale_p(p))
      
      # Select new model
      new.model <- ifelse(nrow(p1) == 1, 1, sample(x = seq_len(nrow(p1)), size = 1, prob = p1$p)) 
      new.model <- rj.obj$mlist[[p1[new.model, ]$model]]

      # Rescale probabilities for reverse move so they sum to 1
      p2 <- rj.obj$config$clust[[1]]
      if(nrow(p2) > 1) p2 <- p2 %>% dplyr::filter(!model == vec_to_model(input.vector = new.model, 
                                             sp.names = rj.obj$dat$species$names))
      p2 <- p2 %>% dplyr::mutate(p_scale = rescale_p(p))
      
      p.jump <- c(p1$p_scale[which(p1$model == vec_to_model(input.vector = new.model, sp.names = rj.obj$dat$species$names))], p2$p_scale[which(p2$model == rj.obj$current.model)])
      
      J <- 1
      
      # Cluster
    } else if (move.type == 2) {
      
      keep.going <- TRUE
      
      while(keep.going){
        
        # Randomly select a number of groups
        new.cluster <- sample(x = rj.obj$config$clust[[2]]$cluster, size = 1, 
                              replace = TRUE, prob = rj.obj$config$clust[[2]]$p_scale)  
        
        # Generate a random grouping with this number of clusters - must be different to the current model
        new.model <- numeric(rj.obj$dat$species$n)
        new.model[sample(x = seq_len(rj.obj$dat$species$n),
                         size = new.cluster, replace = FALSE)] <- seq_len(new.cluster)
        new.model[new.model == 0] <- sample(x = seq_len(new.cluster), 
                                            size = sum(new.model == 0), replace = TRUE)
        new.model <- relabel(new.model)
        
        if(!identical(vec_to_model(input.vector = new.model,
                                   sp.names = rj.obj$dat$species$names), rj.obj$current.model)) 
          keep.going <- FALSE}
      
      # Probability of jump = prob(choose a cluster) x prob(choose model | cluster)
      # A Stirling number of the second kind (or Stirling partition number) is the
      # number of ways to partition a set of n objects into k <non-empty> subsets.
      p.jump <- c(
        rj.obj$config$clust[[2]][rj.obj$config$clust[[2]]$cluster == new.cluster, ]$p_scale *
          (1 / copula::Stirling2(n = rj.obj$dat$species$n, k = new.cluster)),
        
        rj.obj$config$clust[[2]][rj.obj$config$clust[[2]]$cluster == nb_groups(vec = rj.obj$mlist[[rj.obj$current.model]]), ]$p_scale *
          (1 / (copula::Stirling2(n = rj.obj$dat$species$n, 
                                     k = nb_groups(vec = rj.obj$mlist[[rj.obj$current.model]]))))
      )
      
      J <- 1 
      
      # Random
    } else if (move.type == 3) {
      
      model.space <- numbers::bell(rj.obj$dat$species$n) - 1
      
      keep.going <- TRUE
      while (keep.going) {
        new.model <- sample(x = seq_len(rj.obj$dat$species$n), 
                            size = rj.obj$dat$species$n, replace = TRUE) %>% relabel()
        if(!rj.obj$current.model == vec_to_model(input.vector = new.model, sp.names = rj.obj$dat$species$names)) keep.going <- FALSE}
      
      # Probability of jump = 1/(number of models - 1)
      p.jump <- c(1 / model.space, 1 / model.space)
      J <- 1 
    }
    
    n <- N_groups(rj.obj = rj.obj, vec = new.model)
    if (nb_groups(new.model) > nb_groups(rj.obj$mlist[[rj.obj$current.model]])) split.merge <- 1 
    if (nb_groups(new.model) < nb_groups(rj.obj$mlist[[rj.obj$current.model]])) split.merge <- 2
    if (nb_groups(new.model) == nb_groups(rj.obj$mlist[[rj.obj$current.model]])) split.merge <- 0
    
    species.from <- lapply(X = unique(rj.obj$mlist[[rj.obj$current.model]]),
                           FUN = function(x) which(rj.obj$mlist[[rj.obj$current.model]] == x))
    species.to <- lapply(X = unique(new.model), FUN = function(x) which(new.model == x))
    rjgroup <- NA
    
    # Split / merge
  } else {
    
    # Decide whether to split (1) or merge (2)
    if (rj.obj$current.model == vec_to_model(input.vector = rep(1, rj.obj$dat$species$n),
                                           sp.names = rj.obj$dat$species$names)) {
      
      split.merge <- 1
      p.splitmerge <- 1
      
    } else if (rj.obj$current.model == vec_to_model(input.vector = seq_len(rj.obj$dat$species$n),
                                                  sp.names = rj.obj$dat$species$names)) {
      
      split.merge <- 2
      p.splitmerge <- 1
      
    } else {
      
      # Choose between a split and a merge
      split.merge <- sample.int(n = 2, size = 1, prob = rj.obj$mcmc$move$prob)
      p.splitmerge <- rj.obj$mcmc$move$prob[split.merge]
      
    }
    
    # Propose a split
    if (split.merge == 1) {
      
      # Retrieve current model
      M <- rj.obj$mlist[[rj.obj$current.model]]
      
      # Randomly choose a group containing more than one species
      available.groups <- which(l_groups(M) > 1)
      p.choosegroup <- 1 / length(available.groups)
      rjgroup <- sample(x = available.groups, size = 1)
      
      # Identify species in the selected group
      species.from <- which(M == rjgroup)
      
      # Perform a random (single) split - i.e. only allowed to produce two groups out of one
      split.species <- split_count(species.from)
      
      # Probability of that particular split
      p.group <- ifelse(length(split.species) == 2, 1, 1 / length(split.species))
      species.to <- sample(x = split.species, size = 1) %>% purrr::flatten()
      
      # Update partitions
      M2 <- M
      M2[unlist(species.to)] <- unlist(sapply(X = seq_along(species.to), 
                                              FUN = function(x) rep(x, length(species.to[[x]])) + max(M2)))
      new.model <- relabel(M2)
      
      # Sample sizes per group
      n <- N_groups(rj.obj = rj.obj, vec = new.model)[unique(new.model[species.from])]
      if(!length(n) == 2) stop("Length mismatch in n <from: N_groups>")
      
      # Jacobian
      if(phase == 1){
        J <- (n[1] + n[2]) / (n[1] * n[2])
      } else if(phase == 2){
        # J <- ((n[1] + n[2])^2)/(n[1]^2 * n[2]^2)
        J <- ((n[1] + n[2])^3)/(n[1]^3 * n[2]^3)
        # Value of J can be verified by constructing the matrix of partial derivatives
        # m1 <- matrix(c(1, (1/n[1]), 1, -(1/n[2])), 2, 2, byrow = TRUE)
        # m2 <- matrix(0,2,2)
        # abs(det(rbind(cbind(m1, m2, m2), cbind(m2, m1, m2), cbind(m2, m2, m1))))
      }
      
      # Compute model jump probabilities
      p.jump <- c(p.splitmerge * p.choosegroup * p.group)
      
      # Reverse move
      M.reverse <- partitions::vec_to_eq(new.model)
      p.merge.reverse <- 1 / length(merge_count(M.reverse[]))
      if(p.splitmerge == 1) p.splitmerge.reverse <- 1 else p.splitmerge.reverse <- (1 - p.splitmerge)
      p.jump <- c(p.jump, p.splitmerge.reverse * p.merge.reverse)
      
    }
    
    # Propose a merge
    if (split.merge == 2) {
      
      # Retrieve current model
      M <- rj.obj$mlist[[rj.obj$current.model]]
      
      p.merge <- 1 / length(merge_count(partitions::vec_to_eq(M)[]))
      
      rjgroup <- sort(sample(x = seq_len(length(partitions::vec_to_eq(M)[])), size = 2, replace = FALSE))
      species.from <- lapply(X = rjgroup, FUN = function(x) which(M == x))
      
      # Update partitions
      species.to <- sort(do.call(c, partitions::vec_to_eq(M)[rjgroup]))
      names(species.to) <- NULL
      
      M2 <- M
      M2 <- partitions::vec_to_eq(M2)[]
      M2 <- M2[!seq_len(length(partitions::vec_to_eq(M)[])) %in% rjgroup]
      M2 <- append(M2, list(species.to))
      new.model <- format_group(M2) %>% .[["group"]] %>% relabel()
      
      # Compute model jump probabilities
      p.jump <- c(p.splitmerge * p.merge)
      
      # Sample sizes per group
      n <- N_groups(rj.obj = rj.obj, vec = M)[rjgroup]
      
      if(!length(n) == 2) stop("Length mismatch in n <from: N_groups>")
      
      # Jacobian
      if(phase == 1){
        J <- (n[1] * n[2]) / (n[1] + n[2])
      } else if(phase == 2){
        J <- (n[1]^3 * n[2]^3) / ((n[1] + n[2])^3)
        # Value of J can be verified by constructing the matrix of partial derivatives
        # m1 <- matrix(c(n[1]/sum(n), n[2]/sum(n), 0, -n[2]), 2, 2, byrow = TRUE)
        # m2 <- matrix(0,2,2)
        # abs(det(rbind(cbind(m1, m2, m2), cbind(m2, m1, m2), cbind(m2, m2, m1))))
      }
      
      # Reverse move
      available.groups.reverse <- which(l_groups(new.model) > 1)
      p.choosegroup.reverse  <- 1/length(available.groups.reverse)
      
      rjgroup.reverse <- which(purrr::map_lgl(.x = partitions::vec_to_eq(new.model)[], 
                                              .f = ~identical(.x, species.to)))
      
      species.from.reverse <- partitions::vec_to_eq(new.model)[[rjgroup.reverse]]
      split.species.reverse <- split_count(species.from.reverse)
      p.group.reverse <- ifelse(length(split.species.reverse) == 2, 1, 1 / length(split.species.reverse))
      
      if(p.splitmerge == 1) p.splitmerge.reverse <- 1 else p.splitmerge.reverse <- (1 - p.splitmerge)
      p.jump <- c(p.jump, p.splitmerge.reverse * p.choosegroup.reverse * p.group.reverse)
      p.jump
      
    }
  }
  
  return(list(
    type = move.type,
    move = split.merge,
    species = list(from = species.from, to = species.to),
    group = rjgroup,
    model = list(id = new.model, parts = partitions::vec_to_eq(new.model)),
    p = log(p.jump),
    n = n,
    jacobian = log(J)
  ))
}

#' Proposal for jumps between functional forms
#'
#' Generate proposed values for relevant model parameters
#' @param rj.obj Input list.
#' @param from.phase Jump from monophasic to biphasic (1) or from biphasic to monophasic (2).
#' @param scale.adj SD used in calls to rtnorm for the following parameters: alpha, nu, sigma, phi, mu

proposal_ff <- function(rj.obj, from.phase, prop.scale = list(alpha = 10, omega = 1, psi.i = 0.25)){
  
  fprop <- list()
  fprop$to.phase <- (2 - from.phase) + 1
  fprop$t.ij <- rj.obj$t.ij[rj.obj$iter["t.ij"], ]
  
  inits.bi <- purrr::map(
    .x = seq_len(nb_groups(rj.obj$mlist[[rj.obj$current.model]])),
    .f = ~ {
      censored <- rj.obj$dat$obs$censored[rj.obj$dat$species$trials %in%
                                            which(rj.obj$mlist[[rj.obj$current.model]] == .x)]
      y.obs <- rj.obj$dat$obs$y_ij[rj.obj$dat$species$trials %in%
                                     which(rj.obj$mlist[[rj.obj$current.model]] == .x)]
      Lc.obs <- rj.obj$dat$obs$Lc[rj.obj$dat$species$trials %in%
                                    which(rj.obj$mlist[[rj.obj$current.model]] == .x)]
      Rc.obs <- rj.obj$dat$obs$Rc[rj.obj$dat$species$trials %in%
                                    which(rj.obj$mlist[[rj.obj$current.model]] == .x)]
      
      y.obs[is.na(y.obs) & censored == 1] <-
        runif(
          n = sum(is.na(y.obs) & censored == 1),
          min = Rc.obs[is.na(y.obs) & censored == 1],
          max = rj.obj$dat$param$dose.range[2]
        )
      y.obs[is.na(y.obs) & censored == -1] <-
        runif(
          n = sum(is.na(y.obs) & censored == -1),
          min = rj.obj$dat$param$dose.range[1],
          max = Lc.obs[is.na(y.obs) & censored == -1]
        )
      
      sort(y.obs)
    })
  
  # MONOPHASIC TO BIPHASIC 
  
  if(fprop$to.phase == 2){
    
    fprop$k.ij <- rj.obj$k.ij[rj.obj$iter["k.ij"], ]
    
    # L.bound <- purrr::map_dbl(
    #   .x = seq_len(nb_groups(rj.obj$mlist[[rj.obj$current.model]])),
    #   .f = ~max(rj.obj$mu.ij[rj$iter["mu.ij"], , 1][which(rj.obj$mlist[[rj.obj$current.model]][rj.obj$dat$species$trials] == .x)]))
    # 
    # U.bound <- purrr::map_dbl(
    #   .x = seq_len(nb_groups(rj.obj$mlist[[rj.obj$current.model]])),
    #   .f = ~ min(rj.obj$mu.ij[rj$iter["mu.ij"], , 2][which(rj.obj$mlist[[rj.obj$current.model]][rj.obj$dat$species$trials] == .x)]))
    
    fprop$alpha <- rtnorm(n = length(inits.bi), 
                          location = sapply(X = inits.bi, FUN = function(a){min(a) + (max(a)-min(a))/2}), 
                          scale = prop.scale[["alpha"]], 
                          L = rj.obj$dat$param$dose.range[1], 
                          U = rj.obj$dat$param$dose.range[2])[rj.obj$mlist[[rj.obj$current.model]]]
    
    # k_ij
    fprop$k.ij <- 2 - (fprop$alpha[rj.obj$dat$species$trials] > fprop$t.ij)

    # pi_ij
    fprop$pi.ij <- rep(0.5, rj.obj$dat$trials$n)

    # psi_ij
    fprop$psi.ij <- qnorm(fprop$pi.ij) - cov_effects(rj.obj)
    
    # psi.i
    fprop$psi.i_mean <- sapply(X = 1:rj.obj$dat$whales$n, FUN = function(x) mean = mean(fprop$psi.ij[rj.obj$dat$whales$id == x]))
    
    fprop$psi.i <- rnorm(n = length(fprop$psi.i_mean),  mean = fprop$psi.i_mean, sd = prop.scale[["psi.i"]])
    
    # psi and omega
    fprop$psi <- mean(fprop$psi.i)
    
    fprop$omega <- rtnorm(n = 1, location = 1, scale = prop.scale[["omega"]],
                          L = rj.obj$config$priors["omega", 1],
                          U = rj.obj$config$priors["omega", 2])

    
    # mu_ij
    fprop$mu.ij <- matrix(data = NA, nrow = rj.obj$dat$trials$n, ncol = 2)
    fprop$mu.ij[cbind(1:nrow(fprop$mu.ij), fprop$k.ij)] <- fprop$t.ij
    
    L.alpha <- U.alpha <- fprop$alpha[rj.obj$dat$species$trials]
    
    L.alpha[rj.obj$dat$obs$censored == -1 & L.alpha > rj.obj$dat$obs$Lc] <- 
      rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1 & L.alpha > rj.obj$dat$obs$Lc]
    
    U.alpha[rj.obj$dat$obs$censored == 1 & U.alpha < rj.obj$dat$obs$Rc] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1 & U.alpha < rj.obj$dat$obs$Rc]

    fprop$na_1 <- is.na(fprop$mu.ij[, 1])
    fprop$na_2 <- is.na(fprop$mu.ij[, 2])
    
    fprop$mu.ij[is.na(fprop$mu.ij[, 1]), 1] <- 
      runif(n = sum(is.na(fprop$mu.ij[, 1])),
            min = rj.obj$dat$param$dose.range[1],
            max = L.alpha[is.na(fprop$mu.ij[, 1])])
    
    fprop$mu.ij[is.na(fprop$mu.ij[,2]), 2] <- 
      runif(n = sum(is.na(fprop$mu.ij[,2])),
            min = U.alpha[is.na(fprop$mu.ij[,2])],
            max = rj.obj$dat$param$dose.range[2])
    
    fprop$nu <- sapply(X = seq_len(rj.obj$dat$species$n),
           FUN = function(x){
             matrix(colMeans(fprop$mu.ij[rj.obj$dat$species$trials %in% 
                                    which(rj.obj$mlist[[rj$current.model]] == rj.obj$mlist[[rj$current.model]][x]), , drop = FALSE]))})
    
    fprop$tau <- apply(X = fprop$mu.ij, MARGIN = 2, sd) 
    
    fprop$mu <- rj.obj$mu[rj.obj$iter["mu"], ]
    fprop$sigma <- rj.obj$sigma[rj.obj$iter["sigma"]]
    fprop$phi <- rj.obj$phi[rj.obj$iter["phi"]]
    fprop$mu.i <- rj.obj$mu.i[rj.obj$iter["mu.i"], ]

  }
  
  # BIPHASIC TO MONOPHASIC 
  
  if(fprop$to.phase == 1){
    
    fprop$alpha <- rj.obj$alpha[rj.obj$iter["alpha"], ]
    L.alpha <- U.alpha <- fprop$alpha[rj.obj$dat$species$trials]
    L.alpha[rj.obj$dat$obs$censored == -1] <- rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]
    U.alpha[rj.obj$dat$obs$censored == 1] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]
    
    fprop$nu <- if(rj.obj$dat$species$n == 1) matrix(rj.obj$nu[rj.obj$iter["nu"], , ]) else
      t(rj.obj$nu[rj.obj$iter["nu"], , ])
    fprop$mu.ij <- rj.obj$mu.ij[rj.obj$iter["mu.ij"], , ]
    fprop$tau <- rj.obj$tau[rj.obj$iter["tau"], ]
    fprop$omega <- rj.obj$omega[rj.obj$iter["omega"]]
    fprop$psi <- rj.obj$psi[rj.obj$iter["psi"]]
    fprop$psi.i <- rj.obj$psi.i[rj.obj$iter["psi.i"], ]
    fprop$psi.ij <- fprop$psi.i[rj.obj$dat$whales$id] + cov_effects(rj.obj)
    fprop$pi.ij <- rj.obj$pi.ij[rj.obj$iter["pi.ij"], ]
    fprop$k.ij <- rj.obj$k.ij[rj.obj$iter["k.ij"], ]
    
    fprop$psi.i_mean <- sapply(X = 1:rj.obj$dat$whales$n, FUN = function(x) mean = mean(fprop$psi.ij[rj.obj$dat$whales$id == x]))
    
    fprop$t_ij <- fprop$t.ij - cov_effects(rj.obj)
    fprop$mu.i <- unique(sapply(X = dat$whales$id, FUN = function(n) mean(fprop$t_ij[dat$whales$id == n], na.rm = TRUE)))
    
    fprop$mu <- sapply(X = rj.obj$mlist[[rj$current.model]],
          FUN = function(sp) mean(fprop$mu.i[rj.obj$dat$species$id %in% which(rj.obj$mlist[[rj$current.model]] == sp)]))
    
    # This makes the moves to monophasic completely deterministic!
    # Need na.rm = TRUE to avoid numerical issues when individuals have only been subject to one exposure
    fprop$sigma <- mean(unique(sapply(X = dat$whales$id, FUN = function(n) 
      sd(fprop$t.ij[dat$whales$id == n], na.rm = TRUE))), na.rm = TRUE)
    
    fprop$phi <- mean(sapply(X = unique(rj.obj$dat$species$id), FUN = function(n) 
      sd(fprop$mu.i[rj.obj$dat$species$id == n], na.rm = TRUE)), na.rm = TRUE)
    
  }

  fprop$prop.scale <- prop.scale
  fprop$inits.bi <- inits.bi
  fprop$L.alpha <- L.alpha
  fprop$U.alpha <- U.alpha
  return(fprop) 
} # End proposal_ff

#' Proposal for between-model jumps
#'
#' Generate proposed values for relevant model parameters
#' @param rj.obj Input list.
#' @param jump Proposed model jump, as defined by \code{\link{propose_jump}}.
#' @param phase Monphasic (1) or biphasic (2).
#' @param huelsenbeck Defines the type of proposal.

proposal_rj <- function(rj.obj, jump, phase, huelsenbeck = TRUE) {

  # Monophasic
  
  if(phase == 1){
  
  # Current mu
  rj.means <- rj.obj$mu[rj.obj$iter["mu"], ]
  
  # Bootstrap moves
  if(jump$type == 1){
    
    rj.means <- median_data(rj.obj = rj.obj, model.id = jump$model$id)
    
    u <- rtnorm(n = length(unique(jump$model$id)), 
                location = 0, 
                scale = rj.obj$config$prop$dd,
                L = rj.obj$dat$param$dose.range[1] - unique(rj.means),
                U = rj.obj$dat$param$dose.range[2] - unique(rj.means))
    
    rj.means <- rj.means + u[jump$model$id]
    g_u <- NULL
    
    # Cluster and Random moves
  } else if (jump$type %in% c(2, 3)) {
    
    rj.means <- runif(n = nb_groups(jump$model$id), 
                      min = rj.obj$dat$param$dose.range[1],
                      max = rj.obj$dat$param$dose.range[2])[jump$model$id]
    g_u <- u <- NULL
    
    # Split / merge move
  } else if (jump$type == 0){
    
    # Split move
    if (jump$move == 1) {
      
      g_u <- numeric(2)
      b1 <- jump$n[1] * (rj.obj$dat$param$dose.range[2] - rj.means[jump$species$from[1]])
      b2 <- jump$n[1] * (rj.obj$dat$param$dose.range[1] - rj.means[jump$species$from[1]])
      b3 <- -jump$n[2] * (rj.obj$dat$param$dose.range[2] - rj.means[jump$species$from[1]])
      b4 <- -jump$n[2] * (rj.obj$dat$param$dose.range[1] - rj.means[jump$species$from[1]])
      b <- sort(c(b1, b2, b3, b4))
      g_u <- b[2:3]
      # g_u <- c(max(b[b<0]), min(b[b>0]))
      
      if(huelsenbeck){

      u <- runif(n = 1, min = g_u[1], max = g_u[2])
      
      } else {

        # Auxiliary variable for the change of dimensions
        u <- rtnorm(n = 1, location = 0, scale = rj.obj$config$prop$rj, L = g_u[1], U = g_u[2])
        
      }
      
      rj.means[jump$species$to[[1]]] <- rj.means[jump$species$to[[1]]] + (u / jump$n[1])
      rj.means[jump$species$to[[2]]] <- rj.means[jump$species$to[[2]]] - (u / jump$n[2])
      
    } # End move == 1
    
    # Merge move
    if (jump$move == 2) {
      
      # Means
      m <- sapply(X = jump$group, 
                  FUN = function(a) 
                    unique(rj.means[which(rj.obj$mlist[[rj.obj$current.model]] == a)]))
      
      # Calculate weighted average, using sample sizes as weights
      wm <- weighted.mean(x = m, w = jump$n)
      
      # current.means[now.together] <- wm
      rj.means[jump$species$to] <- wm
      
      # The auxiliary variable in a merge move is calculated deterministically
      g_u <- unname(c(-jump$n[1] * wm, jump$n[2] * wm))
      u <- (m[1] - wm) * jump$n[1]
      
      # Or: u <- (wm - m[2]) * jump$n[2]
      # Or: u = (m[2] - m[1]) * ((n[1]*n[2])/(n[1]+n[2]))
      
    } # End move == 2
  } # End random
  
  return(list(out = list(mu = rj.means), aux = list(u = u, g_u = g_u), n = jump$n, huelsenbeck = huelsenbeck))
  
  # Biphasic
  
  } else if(phase == 2){
    
    # Current mu
    rj.means <- rj.obj$nu[rj.obj$iter["nu"], ,]
    alpha.means <- rj.obj$alpha[rj.obj$iter["alpha"], ]
    
    # Bootstrap moves
    if(jump$type == 1){
      
      alpha.means <- median_data(rj.obj = rj.obj, model.id = jump$model$id) 
      z <- rtnorm(n = length(unique(jump$model$id)), location = 0, scale = rj.obj$config$prop$dd,
                  L = rj.obj$config$priors["alpha", 1] - alpha.means,
                  U = rj.obj$config$priors["alpha", 2] - alpha.means) 
      alpha.means <- alpha.means + z[jump$model$id]
      
      u <- rtnorm(n = length(unique(jump$model$id)), 
                 location = 0, 
                 scale = rj.obj$config$prop$dd, 
                 L = rj.obj$config$priors["nu", 1] - unique(alpha.means), 
                 U = 0)
      
      w <- rtnorm(n = length(unique(jump$model$id)), 
                  location = 0, 
                  scale = rj.obj$config$prop$dd, 
                  L = 0,
                  U = rj.obj$config$priors["nu", 2] - unique(alpha.means))

      rj.means <- matrix(c(alpha.means + u[jump$model$id], alpha.means + w[jump$model$id]), nrow = 2)
      
      # rj.med <- median_data(rj.obj = rj.obj, model.id = jump$model$id) 
      # rj.diff <- rj.med - rj.obj$alpha[rj.obj$iter["alpha"], ]
      # rj.complement <- rj.obj$alpha[rj.obj$iter["alpha"], ] - rj.diff
      # rj.means <- matrix(c(rj.med, rj.complement), nrow = 2, byrow = TRUE) %>% 
      #   apply(., MARGIN = 2, sort) + matrix(c(u, w), nrow = 2, byrow = TRUE)
      
      g_u <- g_w <- g_z <- NULL
      
      # Cluster and Random moves
    } else if (jump$type %in% c(2, 3)) {
      
      alpha.means <- runif(n = nb_groups(jump$model$id), 
                           min = rj.obj$config$priors["alpha", 1],
                           max = rj.obj$config$priors["alpha", 2])
      
      rj.means <- matrix(c(runif(n = nb_groups(jump$model$id),
                                     min = rj.obj$dat$param$dose.range[1],
                                     max = alpha.means)[jump$model$id],
                    runif(n = nb_groups(jump$model$id), 
                          min = alpha.means,
                          max = rj.obj$dat$param$dose.range[2])[jump$model$id]), nrow = 2)
      
      g_u <- g_w <- u <- w <-  NULL
      
      # Split / merge move
    } else if (jump$type == 0){
      
      # Split move
      if (jump$move == 1) {
        
        g_u <- g_w <- g_z <- numeric(2)
        
        # Alpha
        
        c1 <- min(jump$n[1] * (c(rj.obj$mu.ij[rj.obj$iter["mu.ij"], rj.obj$dat$species$trials %in% jump$species$to[[1]], 2],
                                 rj.means[, 2],
                                 rj.obj$config$priors["alpha", 2]) - 
                             alpha.means[jump$species$from[1]]))
        
        c2 <- max(jump$n[1] * (c(rj.obj$mu.ij[rj.obj$iter["mu.ij"], rj.obj$dat$species$trials %in% jump$species$to[[1]], 1],
                                 rj.means[, 1],
                                 rj.obj$config$priors["alpha", 1]) - 
                                 alpha.means[jump$species$from[1]]))
        
        c3 <- min(-jump$n[2] * (c(rj.obj$mu.ij[rj.obj$iter["mu.ij"], rj.obj$dat$species$trials %in% jump$species$to[[2]], 1],
                                 rj.means[, 1],
                                 rj.obj$config$priors["alpha", 1]) - 
                                 alpha.means[jump$species$from[1]]))
        
        c4 <- max(-jump$n[2] * (c(rj.obj$mu.ij[rj.obj$iter["mu.ij"], rj.obj$dat$species$trials %in% jump$species$to[[2]], 2],
                                  rj.means[, 2],
                                  rj.obj$config$priors["alpha", 2]) - 
                              alpha.means[jump$species$from[1]]))

        # g_z <- sort(c(max(c2, c4), min(c1, c3)))
        g_z <- sort(c(c1, c2, c3, c4))[2:3]

        if(huelsenbeck){
          z <- runif(n = 1, min = g_z[1], max = g_z[2])
        } else {
          z <- rtnorm(n = 1, location = 0, scale = rj.obj$config$prop$rj, L = g_z[1], U = g_z[2])
        }
        
        alpha.means[jump$species$to[[1]]] <- alpha.means[jump$species$to[[1]]] + (z / jump$n[1])
        alpha.means[jump$species$to[[2]]] <- alpha.means[jump$species$to[[2]]] - (z / jump$n[2])
        
        # nu_1
        
        a1 <- jump$n[1] * (alpha.means[jump$species$from[1]] - rj.means[jump$species$from[1], 1])
        a2 <- jump$n[1] * (rj.obj$dat$param$dose.range[1] - rj.means[jump$species$from[1], 1])
        a3 <- -jump$n[2] * (alpha.means[jump$species$from[1]] - rj.means[jump$species$from[1], 1])
        a4 <- -jump$n[2] * (rj.obj$dat$param$dose.range[1] - rj.means[jump$species$from[1], 1])

        g_u <- sort(c(a1, a2, a3, a4))[2:3]
        
        b1 <- jump$n[1] * (rj.obj$dat$param$dose.range[2] - rj.means[jump$species$from[1], 2])
        b2 <- jump$n[1] * (alpha.means[jump$species$from[1]] - rj.means[jump$species$from[1], 2])
        b3 <- -jump$n[2] * (alpha.means[jump$species$from[1]] - rj.means[jump$species$from[1], 2])
        b4 <- -jump$n[2] * (rj.obj$dat$param$dose.range[2] - rj.means[jump$species$from[1], 2])

        g_w <- sort(c(b1, b2, b3, b4))[2:3]
        
        if(huelsenbeck){

          u <- runif(n = 1, min = g_u[1], max = g_u[2])
          w <- runif(n = 1, min = g_w[1], max = g_w[2])
          
        } else {
          
          # Auxiliary variable for the change of dimensions
          u <- rtnorm(n = 1, location = 0, scale = rj.obj$config$prop$rj, L = g_u[1], U = g_u[2])
          w <- rtnorm(n = 1, location = 0, scale = rj.obj$config$prop$rj, L = g_w[1], U = g_w[2])
          
        }
        
        rj.means[jump$species$to[[1]], 1] <- rj.means[jump$species$to[[1]], 1] + (u / jump$n[1])
        rj.means[jump$species$to[[2]], 1] <- rj.means[jump$species$to[[2]], 1] - (u / jump$n[2])
        rj.means[jump$species$to[[1]], 2] <- rj.means[jump$species$to[[1]], 2] + (w / jump$n[1])
        rj.means[jump$species$to[[2]], 2] <- rj.means[jump$species$to[[2]], 2] - (w / jump$n[2])
        
      } # End move == 1
      
      # Merge move
      if (jump$move == 2) {
        
        # Means
        m <- purrr::map(.x = 1:2,
                        .f = ~sapply(X = jump$group, FUN = function(a) 
                      unique(rj.means[which(rj.obj$mlist[[rj.obj$current.model]] == a), .x])))
        
        m <- append(m, list(sapply(X = jump$group, FUN = function(a) 
          unique(alpha.means[which(rj.obj$mlist[[rj.obj$current.model]] == a)]))))
        
        # Calculate weighted average, using sample sizes as weights
        wm <- purrr::map_dbl(.x = m, .f = ~weighted.mean(x = .x, w = jump$n))
        
        # current.means[now.together] <- wm
        rj.means[jump$species$to, 1] <- wm[1]
        rj.means[jump$species$to, 2] <- wm[2]
        alpha.means[jump$species$to] <- wm[3]
        
        # The auxiliary variable in a merge move is calculated deterministically
        g_u <- unname(c(-jump$n[1] * wm[1], jump$n[2] * wm[1])) 
        g_w <- unname(c(-jump$n[1] * wm[2], jump$n[2] * wm[2]))
        g_z <- unname(c(-jump$n[1] * wm[3], jump$n[2] * wm[3]))
        
        u <- (m[[1]][1] - wm[1]) * jump$n[1]
        w <- (m[[2]][1] - wm[2]) * jump$n[1]
        z <- (m[[3]][1] - wm[3]) * jump$n[1]
        
        # Or: u <- (wm - m[2]) * jump$n[2]
        # Or: u = (m[2] - m[1]) * ((n[1]*n[2])/(n[1]+n[2]))
        
      } # End move == 2
    } # End random

    return(list(out = list(nu = t(rj.means), alpha = alpha.means),
                aux = list(u = u, 
                           w = w, 
                           z = z,
                           g_u = g_u, 
                           g_w = g_w, 
                           g_z = g_z),
                n = jump$n, 
                huelsenbeck = huelsenbeck))
    
  }
}

#' Proposal for within-model jumps
#'
#' Generate proposed values for relevant model parameters
#'
#' @param rj.obj RJ object.
#' @param param.name Parameter name.
#' @param seed Random seed. Only used for testing purposes.

proposal_mh <- function(rj.obj, param.name, seed = NULL) {
    
  # Set the random seed
  if(!is.null(seed)) set.seed(seed) 
  
  # Initialise p.bounds
  p.bounds <- NULL
  
  # Species groups
  ng <- unique(rj.obj$mlist[[rj.obj$current.model]])
  
  # Number of proposal values
  if(param.name %in% c("t.ij", "k.ij", "mu.ij")) N <- rj.obj$dat$trials$n else 
    if(param.name %in% c("mu.i", "psi.i")) N <- rj.obj$dat$whales$n else 
      if(param.name %in% c("mu", "alpha", "nu")) N <- length(ng) else 
          if(param.name %in% c("tau")) N <- 2 else 
            if(param.name %in% c("phi", "sigma", "omega", "psi")) N <- 1 else
              if(param.name %in% rj.obj$dat$covariates$names) N <- rj.obj$dat$covariates$fL[[param.name]]$nparam
  
  if(param.name == "k.ij") {
    
    pval <- (1 - rbinom(n = N, size = 1, prob = rj.obj$config$prop$mh[["k.ij"]])) + 1
    
  } else {
    
    # Mean(s) for proposal(s)
    if (param.name %in% c("t.ij", "mu.i", "alpha", "mu", "psi.i", "tau"))
      m <- rj.obj[[param.name]][rj.obj$iter[param.name], ]
    
    if (param.name %in% c("sigma", "phi", "omega", "psi"))
      m <- rj.obj[[param.name]][rj.obj$iter[param.name]]
    
    if (param.name %in% c("nu", "mu.ij"))
      m <- rj.obj[[param.name]][rj.obj$iter[param.name], , ]
    
    if (param.name %in% rj.obj$dat$covariates$names) 
      m <- rj.obj[[param.name]][rj.obj$iter[param.name], rj.obj$dat$covariates$fL[[param.name]]$index]
    
    # Generate proposals
    # BIPHASIC
    
    if(param.name == "mu.ij"){
      
      p.bounds <- cbind(rep(rj.obj$dat$param$dose.range[1], rj.obj$dat$trials$n),
                        rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials],
                        rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials],
                        rep(rj.obj$dat$param$dose.range[2], rj.obj$dat$trials$n))
      
      p.bounds[, 1][rj$k.ij[rj$iter["k.ij"],] == 1 & 
                      rj.obj$dat$obs$censored == 1 & 
                      rj.obj$dat$obs$Rc > p.bounds[, 1]] <- 
        rj.obj$dat$obs$Rc[rj$k.ij[rj$iter["k.ij"],] == 1 & 
                            rj.obj$dat$obs$censored == 1 & 
                            rj.obj$dat$obs$Rc > p.bounds[, 1]]
      
      p.bounds[, 3][rj$k.ij[rj$iter["k.ij"],] == 2 & 
                      rj.obj$dat$obs$censored == 1 & 
                      rj.obj$dat$obs$Rc > p.bounds[, 3]] <- 
        rj.obj$dat$obs$Rc[rj$k.ij[rj$iter["k.ij"],] == 2 & 
                            rj.obj$dat$obs$censored == 1 & 
                            rj.obj$dat$obs$Rc > p.bounds[, 3]]

      p.bounds[, 2][rj$k.ij[rj$iter["k.ij"],] == 1 & 
                      rj.obj$dat$obs$censored == -1 & 
                      rj.obj$dat$obs$Lc < p.bounds[, 2]] <- 
        rj.obj$dat$obs$Lc[rj$k.ij[rj$iter["k.ij"],] == 1 & 
                            rj.obj$dat$obs$censored == -1 & 
                            rj.obj$dat$obs$Lc < p.bounds[, 2]]
      
      p.bounds[, 4][rj$k.ij[rj$iter["k.ij"],] == 2 & 
                      rj.obj$dat$obs$censored == -1 & 
                      rj.obj$dat$obs$Lc < p.bounds[, 4]] <- 
        rj.obj$dat$obs$Lc[rj$k.ij[rj$iter["k.ij"],] == 2 & 
                            rj.obj$dat$obs$censored == -1 & 
                            rj.obj$dat$obs$Lc < p.bounds[, 4]]
      
      # Sense checks
      if(sum(p.bounds[, 2] < p.bounds[, 1]) > 0) stop("Inconsistent bounds for mu_ij")
      if(sum(p.bounds[, 4] < p.bounds[, 3]) > 0) stop("Inconsistent bounds for mu_ij")
      
      # data.frame(p.bounds,
      #            alpha = unname(rj.obj$alpha[rj.obj$iter["alpha"],rj.obj$dat$species$trials]),
      #            k = rj.obj$k.ij[rj.obj$iter["k.ij"], ],
      #            censored = rj$dat$obs$censored,
      #            Rc = rj$dat$obs$Rc,
      #            Lc = rj$dat$obs$Lc, row.names = NULL)
      
      
    # L.1 <- rep(rj.obj$dat$param$dose.range[1], rj.obj$dat$trials$n)
    # U.1 <- L.2 <- rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials]
    # U.2 <- rep(rj.obj$dat$param$dose.range[2], rj.obj$dat$trials$n)
    # 
    # L.2[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2 & rj.obj$dat$obs$censored == 1] <-
    #   apply(X = cbind(L.2[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2 & rj.obj$dat$obs$censored == 1],
    #                   rj.obj$dat$obs$Rc[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2 & rj.obj$dat$obs$censored == 1]),
    #         MARGIN = 1, max)
    # 
    # U.1[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1 & rj.obj$dat$obs$censored == -1] <-
    #   apply(X = cbind(U.1[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1 & rj.obj$dat$obs$censored == -1],
    #                   rj.obj$dat$obs$Lc[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1 & rj.obj$dat$obs$censored == -1]),
    #         MARGIN = 1, min)
    # 
    # data.frame(alpha = rj.obj$alpha[rj.obj$iter["alpha"],rj.obj$dat$species$trials],
    #                                 L.1, U.1, L.2, U.2, k = rj.obj$k.ij[rj.obj$iter["k.ij"], ],censored = rj$dat$obs$censored, Rc = rj$dat$obs$Rc, Lc = rj$dat$obs$Lc)
    
    } else if(param.name == "t.ij" & rj.obj$config$function.select) {
      
      p.bounds <- cbind(rep(rj.obj$dat$param$dose.range[1], rj.obj$dat$trials$n),
                      rep(rj.obj$dat$param$dose.range[2], rj.obj$dat$trials$n))
      
      p.bounds[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2, 1] <- rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials][rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2]
      
      p.bounds[rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1, 2] <- rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials][rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1]
      
      p.bounds[, 1][rj.obj$dat$obs$censored == 1 & rj.obj$dat$obs$Rc > p.bounds[, 1]] <- 
        rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1 & rj.obj$dat$obs$Rc > p.bounds[, 1]]
      
      p.bounds[, 2][rj.obj$dat$obs$censored == -1 & rj.obj$dat$obs$Lc < p.bounds[, 2]] <- 
        rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1 & rj.obj$dat$obs$Lc < p.bounds[, 2]]
      
      # Sense checks
      if(sum(p.bounds[, 2] < p.bounds[, 1]) > 0) stop("Inconsistent bounds for t_ij")
      
      # L.lower <- L.upper <- rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials]
      # 
      # L.lower[rj.obj$dat$obs$Rc > L.lower & rj.obj$dat$obs$censored == 1] <-
      #   rj.obj$dat$obs$Rc[rj.obj$dat$obs$Rc > L.lower & rj.obj$dat$obs$censored == 1]
      # 
      # L.lower[rj.obj$dat$obs$Lc < L.lower & rj.obj$dat$obs$censored == -1] <-
      #   rj.obj$dat$obs$Lc[rj.obj$dat$obs$Lc < L.lower & rj.obj$dat$obs$censored == -1]
      # 
      # L.upper[rj.obj$dat$obs$Rc > L.upper & rj.obj$dat$obs$censored == 1] <-
      #   rj.obj$dat$obs$Rc[rj.obj$dat$obs$Rc > L.upper & rj.obj$dat$obs$censored == 1]
      # 
      # L.upper[rj.obj$dat$obs$Lc < L.upper & rj.obj$dat$obs$censored == -1] <-
      #   rj.obj$dat$obs$Lc[rj.obj$dat$obs$Lc < L.upper & rj.obj$dat$obs$censored == -1]
      
     # Check that all is correct
     # data.frame(alpha = rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials],
     #            L = cbind(rj.obj$dat$param$dose.range[1], L.lower)[cbind(1:rj.obj$dat$trials$n, rj.obj$k.ij[rj.obj$iter["k.ij"],])], U = cbind(L.upper, rj.obj$dat$param$dose.range[2])[cbind(1:rj.obj$dat$trials$n, rj.obj$k.ij[rj.obj$iter["k.ij"],])], k = rj$k.ij[rj.obj$k.ij[rj.obj$iter["k.ij"], ]], censored = rj$dat$obs$censored)
     
    } else {
      
      pn <- ifelse(param.name %in% c("t.ij", "mu.i"), "mu", param.name)
      lower.limit <- rep(rj.obj$config$priors[pn, 1], N)
      upper.limit <- rep(rj.obj$config$priors[pn, 2], N)
      if(param.name == "t.ij"){
        lower.limit[rj.obj$dat$obs$censored == 1] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]
        upper.limit[rj.obj$dat$obs$censored == -1] <- rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]
      }
    }
    
    if(param.name == "t.ij" & rj.obj$config$function.select){
      
      pval <- rtnorm(n = N, 
                     location = m, 
                     scale = rj.obj$config$prop$mh[[param.name]], 
                     L = p.bounds[, 1],
                     U = p.bounds[, 2])
      
      # pval <- rtnorm(n = N, 
      #                location = m, 
      #                scale = rj.obj$config$prop$mh[[param.name]], 
      #                L = cbind(rj.obj$dat$param$dose.range[1], L.lower)[cbind(1:rj.obj$dat$trials$n, rj.obj$k.ij[rj.obj$iter["k.ij"],])],
      #                U = cbind(L.upper, rj.obj$dat$param$dose.range[2])[cbind(1:rj.obj$dat$trials$n, rj.obj$k.ij[rj.obj$iter["k.ij"],])])
    
    } else if(param.name == "t.ij" & !rj.obj$config$function.select | param.name == "mu.i"){
      
      pval <- rtnorm(n = N, 
                     location = m, 
                     scale = rj.obj$config$prop$mh[[param.name]], 
                     L = lower.limit,
                     U = upper.limit)
      
    } else if(param.name == "mu.ij"){

      pval <- unname(cbind(rtnorm(n = N,
                                  location = m[, 1],
                                  scale = rj.obj$config$prop$mh[[param.name]],
                                  L = p.bounds[, 1], U = p.bounds[, 2]),
                           
                           rtnorm(n = N,
                                  location = m[, 2],
                                  scale = rj.obj$config$prop$mh[[param.name]],
                                  L = p.bounds[, 3], U = p.bounds[, 4])))
      
      # pval <- unname(cbind(rtnorm(n = N,
      #                             location = m[, 1],
      #                             scale = rj.obj$config$prop$mh[[param.name]],
      #                             L = L.1, U = U.1),
      #                      
      #                      rtnorm(n = N,
      #                             location = m[, 2],
      #                             scale = rj.obj$config$prop$mh[[param.name]],
      #                             L = L.2, U = U.2)))
      
    } else if(param.name == "nu"){
      
      if(is.null(dim(m))){
        
        m1 <- m[1] 
        m2 <- m[2]
        pval <- rbind(rnorm(n = N, mean = unique(m1), sd = rj.obj$config$prop$mh[[param.name]]),
              rnorm(n = N, mean = unique(m2), sd = rj.obj$config$prop$mh[[param.name]]))
        
      } else {
        
        m1 <- m[,1]
        m2 <- m[,2]
      
        pval <- rbind(rnorm(n = N, mean = unique(m1), sd = rj.obj$config$prop$mh[[param.name]]),
                      rnorm(n = N, mean = unique(m2), sd = rj.obj$config$prop$mh[[param.name]]))
      }
      
    } else {
      
      pval <- rnorm(n = N, mean = unique(m), sd = rj.obj$config$prop$mh[[param.name]])
      
    }
    
    # Additions
    index <- cbind(ng, ord = 1:length(ng))
    
    if (param.name == "mu") {
      pval <- pval[sapply(
        X = rj.obj$mlist[[rj.obj$current.model]],
        FUN = function(x) index[ng == x, 2]
      )]
    }
    
    if(param.name %in% rj.obj$dat$covariates$names){
      if(rj.obj$dat$covariates$fL[[param.name]]$nL >= 2) 
        pval <- c(0, pval)
    }
  } # end if k.ij
  
  if(param.name == "alpha") pval <- pval[rj.obj$mlist[[rj.obj$current.model]]]
  if(param.name == "nu") pval <- pval[, rj.obj$mlist[[rj.obj$current.model]], drop = FALSE]
  
  if(any(is.infinite(pval))) stop("Infinite values generated")
  res <- list(pval)
  names(res) <- param.name
  res <- append(res, list(p.bounds = p.bounds))
  return(res)
}

# Proposal densities ----------------------------------------------------------------

#' Proposal densities for functional form jumps.
#'
#'@param rj.obj Input object.
#' @param param Parameter values.

propdens_ff <- function(rj.obj, param){

  logprop.mono <- 0
 
  logprop.bi <- sum(dtnorm(x = unique(param$alpha),
               location = unique(sapply(X = param$inits.bi, FUN = function(a){min(a) + (max(a)-min(a))/2})[rj.obj$mlist[[rj.obj$current.model]]]),
               scale = param$prop.scale[["alpha"]], 
               L = rj.obj$dat$param$dose.range[1],
               U = rj.obj$dat$param$dose.range[2],
               # L = unique(param$L.bound[rj.obj$mlist[[rj.obj$current.model]]]),
               # U = unique(param$U.bound[rj.obj$mlist[[rj.obj$current.model]]]),
               log = TRUE)) +
    
    dtnorm(x = param$omega, location = 1,
           scale = param$prop.scale[["omega"]],
           L = rj.obj$config$priors["omega", 1],
           U = rj.obj$config$priors["omega", 2],
           log = TRUE) +
    
    sum(dnorm(x = param$psi.i, mean = param$psi.i_mean, sd = param$prop.scale[["psi.i"]], log = TRUE)) +
    
    sum(dunif(x = param$mu.ij[param$na_1, 1],
              min = rj.obj$dat$param$dose.range[1],
              max = param$L.alpha[param$na_1], log = TRUE)) +

    sum(dunif(x = param$mu.ij[param$na_2, 2],
              min = param$U.alpha[param$na_2],
              max = rj.obj$dat$param$dose.range[2], log = TRUE))

  return(list(mono = logprop.mono, bi = logprop.bi))
  
}

#' Proposal densities for between-model jumps
#'
#' @param rj.obj Input object.
#' @param param Parameter values.
#' @param jump Proposed between-model jump.
#' @param phase Monophasic (1) or biphasic (2).

propdens_rj <- function(rj.obj, param, jump, phase) {
  
  # Adapted to work with Uniform proposal distribution as per Huelsenbeck et al.
  # aux is the value of the auxiliary variable 'u' drawn by proposal_rj
  # during a proposed split/merge move. u does not exist when the independent
  # sampler is used.
  
  if(phase == 1) {
    
  # Cluster and Random moves
  if(jump$type %in% c(2, 3)){
    
    loglik.unif <- unlist(purrr::map(.x = list(param$out$mu, rj.obj$mu[rj.obj$iter["mu"], ]),
                              .f = ~sum(dunif(x = unique(.x), 
                         min = rj.obj$dat$param$dose.range[1],
                         max = rj.obj$dat$param$dose.range[2],
                         log = TRUE))))
    loglik.unif[is.infinite(loglik.unif)] <- -100000
    return(loglik.unif)
    
    # Bootstrap moves
  } else if (jump$type == 1) {
    
    loglik.norm.forward <- dtnorm(x = unique(param$aux$u),
                                 location = 0, 
                                 scale = rj.obj$config$prop$dd, 
                                 L = rj.obj$dat$param$dose.range[1] - unique(param$out$mu),
                                 U = rj.obj$dat$param$dose.range[2] - unique(param$out$mu),
                                 log = TRUE)
    
    if (any(abs(loglik.norm.forward) == Inf)) loglik.norm.forward <- -100000
    loglik.norm.forward <- sum(loglik.norm.forward)
    
    loglik.norm.backward <- dtnorm(x = unique(median_data(rj.obj = rj.obj,
                                                model.id = rj.obj$mlist[[rj.obj$current.model]]) - 
                                    rj.obj$mu[rj.obj$iter["mu"], ]),
                                  location = 0, 
                                  scale = rj.obj$config$prop$dd,
                                  L = rj.obj$dat$param$dose.range[1] - unique(rj.obj$mu[rj.obj$iter["mu"], ]),
                                  U = rj.obj$dat$param$dose.range[2] - unique(rj.obj$mu[rj.obj$iter["mu"], ]),
                                  log = TRUE)
    
    if (any(abs(loglik.norm.backward) == Inf)) loglik.norm.backward <- -100000
    loglik.norm.backward <- sum(loglik.norm.backward)
    
    return(c(loglik.norm.backward, loglik.norm.forward))
    
    # Split / merge moves
  } else if (jump$type == 0) {
    
    if(param$huelsenbeck){
    loglik <- dunif(x = unique(param$aux$u), 
                    min = param$aux$g_u[1], 
                    max = param$aux$g_u[2], 
                    log = TRUE)
    } else {
    loglik <- dtnorm(x = unique(param$aux$u), 
                     location = 0,
                     scale = rj.obj$config$prop$rj, 
                     L = param$aux$g_u[1],
                     U = param$aux$g_u[2],
                     log = TRUE)
    }
    
    if (any(abs(loglik) == Inf)) loglik <- -100000
    # Return two values for the ratio of proposal densities
    # Depending on the type of move, either the numerator or
    # denominator is 0 (as we work in logs)
    if (jump$move == 1) return(c(0, loglik)) else return(c(loglik, 0)) 
  }
    
  } else if(phase ==2) {
    
    if(jump$type %in% c(2, 3)){
      
      loglik.unif <- unlist(purrr::map(.x = list(param$out$alpha, rj.obj$alpha[rj.obj$iter["alpha"], ]),
                                       .f = ~sum(dunif(x = unique(.x), 
                                                   min = rj.obj$config$priors["alpha", 1],
                                                   max = rj.obj$config$priors["alpha", 2],
                                                   log = TRUE)))) +
        
        unlist(purrr::map2(.x = list(param$out$nu[1, ], rj.obj$nu[rj.obj$iter["nu"], , 1]),
                                        .y = list(param$out$alpha, rj.obj$alpha[rj.obj$iter["alpha"], ]),
                                .f = ~sum(dunif(x = unique(.x), 
                                min = rj.obj$dat$param$dose.range[1],
                                max = unique(.y),
                                log = TRUE)))) +
        
        unlist(purrr::map2(.x = list(param$out$nu[2, ], rj.obj$nu[rj.obj$iter["nu"], , 2]),
                           .y = list(param$out$alpha, rj.obj$alpha[rj.obj$iter["alpha"], ]),
                           .f = ~sum(dunif(x = unique(.x), 
                                           min = unique(.y),
                                           max = rj.obj$dat$param$dose.range[2],
                                           log = TRUE))))
      
      loglik.unif[is.infinite(loglik.unif)] <- -100000
      return(loglik.unif)

    } else if(jump$type == 1){
      
      loglik.norm.forward <- dtnorm(x = param$aux$z,
                                   location = 0, 
                                   scale = rj.obj$config$prop$dd, 
                                   L = rj.obj$config$priors["alpha", 1] - unique(param$out$alpha),
                                   U = rj.obj$config$priors["alpha", 2] - unique(param$out$alpha),
                                   log = TRUE) + 
        dtnorm(x = param$aux$u, 
               location = 0, 
               scale = rj.obj$config$prop$dd, 
               L = rj.obj$config$priors["nu", 1] - unique(param$out$alpha), 
               U = 0, 
               log = TRUE) + 
        
        dtnorm(x = param$aux$w, 
               location = 0, 
               scale = rj.obj$config$prop$dd, 
               L = 0,
               U = rj.obj$config$priors["nu", 2] - unique(param$out$alpha),
               log = TRUE)
      
      if (any(abs(loglik.norm.forward) == Inf)) loglik.norm.forward <- -100000
      loglik.norm.forward <- sum(loglik.norm.forward)
      
      z <- rj.obj$alpha[rj.obj$iter["alpha"], ] - median_data(rj.obj = rj.obj, model.id = rj.obj$mlist[[rj.obj$current.model]])
      u <- rj.obj$nu[rj.obj$iter["nu"], , 1] - rj.obj$alpha[rj.obj$iter["alpha"], ]
      w <- rj.obj$nu[rj.obj$iter["nu"], , 2] - rj.obj$alpha[rj.obj$iter["alpha"], ]
      
      loglik.norm.backward <- 
        dtnorm(x = unique(z), 
               location = 0, 
               scale = rj.obj$config$prop$dd, 
               L = rj.obj$config$priors["alpha", 1] - unique(rj.obj$alpha[rj.obj$iter["alpha"], ]),
               U = rj.obj$config$priors["alpha", 2] - unique(rj.obj$alpha[rj.obj$iter["alpha"], ]),
               log = TRUE) + 
        
        dtnorm(x = unique(u), location = 0, scale = rj.obj$config$prop$dd, 
               L = rj.obj$config$priors["nu", 1] - unique(rj.obj$alpha[rj.obj$iter["alpha"], ]),
               U = 0, log = TRUE) + 
        
        dtnorm(x = unique(w), location = 0, scale = rj.obj$config$prop$dd, 
               L = 0, U = rj.obj$config$priors["nu", 2] - unique(rj.obj$alpha[rj.obj$iter["alpha"], ]),
               log = TRUE)
      
      if (any(abs(loglik.norm.backward) == Inf)) loglik.norm.backward <- -100000
      loglik.norm.backward <- sum(loglik.norm.backward)
      
      return(c(loglik.norm.backward, loglik.norm.forward))
      
    } else if (jump$type == 0) {
      
      if(param$huelsenbeck){
        
        loglik <- purrr::map2_dbl(.x = unlist(param$aux[c("u", "w", "z")]),
                              .y = param$aux[c("g_u", "g_w", "g_z")],
                              .f = ~dunif(x = unique(.x),
                                          min = .y[1], 
                                          max = .y[2], 
                                          log = TRUE))
          
      } else {
        
        loglik <- dtnorm(x = unlist(param$aux[c("u", "w", "z")]), 
                        location = 0, 
                        scale = rj.obj$config$prop$rj, 
                        L = sapply(param$aux[c("g_u", "g_w", "g_z")], "[[", 1),
                        U = sapply(param$aux[c("g_u", "g_w", "g_z")], "[[", 2),
                        log = TRUE)
      }
      
      if (any(abs(loglik) == Inf)) loglik <- -100000
      if (jump$move == 1) return(c(0, sum(loglik))) else return(c(sum(loglik), 0)) 
    }
    
  }
}

#' Proposal densities for between-model jumps
#'
#' @param param.name Parameter name
#' @param dest Proposed value(s)
#' @param orig Current value(s) 

propdens_mh <- function(rj.obj, param.name, dest, orig, bounds = NULL) {

  # if(param.name == "mu.ij"){
    
    # L.lower <- L.upper <- rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials]
    # 
    # L.lower[rj.obj$dat$obs$Lc < L.lower & rj.obj$dat$obs$censored == -1] <-
    #   rj.obj$dat$obs$Lc[rj.obj$dat$obs$Lc < L.lower & rj.obj$dat$obs$censored == -1]
    # 
    # L.upper[rj.obj$dat$obs$Rc > L.upper & rj.obj$dat$obs$censored == 1] <-
    #   rj.obj$dat$obs$Rc[rj.obj$dat$obs$Rc > L.upper & rj.obj$dat$obs$censored == 1]
    
  # } else {
    
  pn <- ifelse(param.name %in% c("t.ij", "mu.i"), "mu", param.name)
  lower.limit <- rep(rj.obj$config$priors[pn, 1], length(dest))
  upper.limit <- rep(rj.obj$config$priors[pn, 2], length(dest))
  
  if(param.name == "t.ij"){
    lower.limit[rj.obj$dat$obs$censored == 1] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]
    upper.limit[rj.obj$dat$obs$censored == -1] <- rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]
  }
  
  # }

  if(param.name == "mu.ij") {
    
    loglik <- unname(cbind(dtnorm(x = dest[, 1], 
                                  location = orig[, 1],
                                  scale = rj.obj$config$prop$mh[[param.name]], 
                                  L = bounds[, 1], U = bounds[, 2],
                                  log = TRUE),
                           dtnorm(x = dest[, 2], 
                                  location = orig[, 2],
                                  scale = rj.obj$config$prop$mh[[param.name]], 
                                  bounds[, 3], U = bounds[, 4], 
                                  log = TRUE)))
    
   } else if(param.name == "t.ij" & rj.obj$config$function.select) {
     
     loglik <- dtnorm(x = dest, 
                      location = orig, 
                      scale = rj.obj$config$prop$mh[[param.name]], 
                      L = bounds[, 1],
                      U = bounds[, 2],
                      log = TRUE)
    
   } else if(param.name %in% c("t.ij", "mu.i")) {
    
    loglik <- dtnorm(x = dest, 
                     location = orig, 
                     scale = rj.obj$config$prop$mh[[param.name]], 
                     L = lower.limit,
                     U = upper.limit,
                     log = TRUE)
  } else {
    
    loglik <- dnorm(x = dest, mean = orig, sd = rj.obj$config$prop$mh[[param.name]], log = TRUE) 
    
  }
  
  if (any(abs(loglik) == Inf)) loglik <- -100000
  if(param.name %in% c("t.ij", "mu.i", "mu.ij")) return(loglik) else return(sum(loglik))
}

# Trace ----------------------------------------------------------------

#' Print method for objects of class \code{rjtrace}
#' @export
print.rjtrace <- function(rj.obj){
  print(summary(rj.obj$trace))
}

# check <- function(x, ...) { UseMethod("check") }

# Dose-response ----------------------------------------------------------------

#' Print method for objects of class \code{dose_response}
#' @export
print.dose_response <- function(dr.object){
  if(!"dose_response" %in% class(dr.object)) stop("dr.object must be of class <dose_response>")
  print(dr.object$ranks)
}


# Gibbs Variable Selection ------------------------------------------------

# Functions written by Dina Sadykova

is.even <- function(x){ x %% 2 == 0 }

"%!in%" <- function(x, table) { match(x, table, nomatch = 0) == 0 }

# All possible A+B combinations (divide all species/signals into 2 groups)
combAB.fun <- function(n){
  
  # All combinations of two groups
  id1 <- unlist(lapply(2:(n-1), function(x) combn(1:n, x, simplify = FALSE)), recursive = FALSE) 
  
  Ind <- diag(n)
  IndAB <- matrix(NA, ncol = n, nrow = 2*length(id1))
  
  for (i in 1:length(id1)){
    IndAB[(2*i-1),] <- colSums(Ind[id1[[i]],])
    IndAB[(2*i),] <- 1-colSums(Ind[id1[[i]],])
  }
  
  # Remove duplicate models
  IndAB <- IndAB[!duplicated(IndAB), ]
  n.mod2 <- rep(2:(1+dim(IndAB)[1]/2), each = 2)
  return(list(IndAB,n.mod2))
}


# All combinations
# 1 group and the rest is treated each separately 
# (for example, AB+C+D+E or ABC+D+E for 5 species)
# it is working for n>=4
combAB_C_D.fun <- function(n){
  
  # All combinations
  id1 <- unlist(lapply(2:(n-1), function(x)combn(1:n, x, simplify = FALSE)), recursive = FALSE) 
  
  Ind <- diag(n)
  IndAB <- matrix(NA, ncol = n, nrow = 0)
  n.mod3 <- NA
  
  for (i in 1:length(id1)){
    if ((n-length(id1[[i]]))>1){
      idx = which(c(1:n) %!in% id1[[i]])
      IndT <- rbind(colSums(Ind[id1[[i]],]),Ind[idx,]) 
      IndAB <- rbind(IndAB,IndT)
      n.mod3 <- c(n.mod3,rep(i,dim(IndT)[1]))
    }
  }
  
  n.mod3 <- n.mod3[-1]
  return(list(IndAB,n.mod3))
}


# Function to find all possible ways to split a list of elements
# Split into a given number of groups of the SAME size (x=1:(even number))
# without duplicate models
comb.groups  <-  function(x, gr){
  
  nx  <-  length(x)
  ning  <-  nx/gr
  
  group1  <-  rbind(matrix(rep(x[1],choose(nx-1,ning-1)),nrow=1),combn(x[-1],ning-1))
  ng  <-  ncol(group1)
  
  if(gr > 2){
    out  <-  vector('list',ng)    
    for(i in seq_len(ng)){
      other  <-  comb.groups(setdiff(x,group1[,i]),gr=gr-1)
      out[[i]]  <-  lapply(seq_along(other), function(j) cbind(group1[,i],other[[j]]))
    }
    out  <-  unlist(out,recursive=FALSE)
  } else {
    other  <-  lapply(seq_len(ng),function(i) matrix(setdiff(x,group1[,i]),ncol=1))
    out  <-  lapply(seq_len(ng),function(i) cbind(group1[,i],other[[i]])
    )
  }
  out    
}


# Index matrix for comb.groups function (few groups)
combAB.ind <- function(n,id){
  Ind <- diag(n)
  IndAB <- matrix(NA,ncol=n,nrow=dim(id)[2])
  for (i in 1:dim(IndAB)[1]){
    if (is.null(dim(Ind[id[,i],]))) IndAB[i,] <- Ind[id[,i],] else IndAB[i,] <- colSums(Ind[id[,i],])
  }  
  return(IndAB)
}


# Index matrix2 for comb.groups function (few groups+rest is treated each separately)
combAB_C_D.ind <- function(n,id){
  Ind <- diag(n)
  IndAB <- IndT <- matrix(NA,ncol=n,nrow=0)
  idx = which(c(1:n) %!in% id)
  for (i in 1:dim(id)[2]){
    if (is.null(dim(Ind[id[,i],]))){
      IndT <- rbind(IndT,Ind[id[,i],])} else {
        IndT <- rbind(IndT,colSums(Ind[id[,i],])) 
      }
  }
  IndT <- rbind(IndT,Ind[idx,])  
  IndAB <- rbind(IndAB,IndT)  
  
  return(IndAB)
}


# Function to find all possible combinations for n>=5 using the previous functions 
all.comb <- function(n){
  
  id1 <- unlist(lapply(2:(n-1),function(x)combn(1:n,x,simplify=F)),recursive=F) 
  id.not1 <- list()
  for (i in 1:length(id1)) id.not1[[i]] <- which(c(1:n) %!in% id1[[i]])
  id1.length <- id.not1.length <- rep(0,length(id1))
  for (i in 1:length(id1.length)){
    id1.length[i] <- length(id1[[i]])
    id.not1.length[i] <- length(id.not1[[i]])
  }
  
  id1 <- id1[id.not1.length>2]
  id.not1 <- id.not1[id.not1.length>2]
  
  Ind <- diag(n)
  IndALL2 <- matrix(NA,ncol=n,nrow=0)
  c.l.arch <- mat.cg.l <- mat.c.list <- list()
  n.mod4 <- 0
  c.ind <- 0
  n.comb.run <- floor(n/2)
  for (n.mod in 2:min(3,max(2,(n.comb.run-1)))){
    for (i in 1:length(id1)){    
      idx <- which(c(1:n) %!in% id1[[i]])  
      if (is.even(length(idx))) {      
        len.comb <- length(comb.groups(idx,n.mod))
        c.len <- matrix(0,nrow=len.comb,ncol=dim(comb.groups(idx,n.mod)[[1]])[2]+1)
        c.len[,1] <- id1.length[i]
        for (hh in 1:len.comb){
          for (hv in 2:(dim(comb.groups(idx,n.mod)[[1]])[2]+1)){
            c.len[hh,hv] <- length((comb.groups(idx,n.mod)[[hh]][,hv-1][comb.groups(idx,n.mod)[[hh]][,hv-1]!=0]))
          }
        } 
        c.len <- unique(t(apply(c.len,1,sort))) 
        c.l.arch <- c(c.l.arch,list(c.len))
        ind <- any(sapply(c.l.arch, function(x, want) isTRUE(all.equal(x, want)), c.len))        
      } else {
        len.comb <- length(comb.groups(c(idx,0),n.mod))
        c.len <- matrix(0,nrow=len.comb,ncol=dim(comb.groups(c(idx,0),n.mod)[[1]])[2]+1)
        c.len[,1] <- id1.length[i]
        for (hh in 1:len.comb){
          for (hv in 2:(dim(comb.groups(c(idx,0),n.mod)[[1]])[2]+1)){
            c.len[hh,hv] <- length((comb.groups(c(idx,0),n.mod)[[hh]][,hv-1][comb.groups(c(idx,0),n.mod)[[hh]][,hv-1]!=0]))
          }
        }        
        c.len <- unique(t(apply(c.len,1,sort))) 
        ind <- any(sapply(c.l.arch, function(x, want) isTRUE(all.equal(x, want)), c.len))
        if (i>=2) if (id1.length[i]!=id1.length[i-1]) c.l.arch <- c(c.l.arch,list(c.len))
      }
      for (j in 1:len.comb){
        #if (ind==FALSE){
        if (is.even(length(idx))){
          comb.group <- rbind(colSums(Ind[id1[[i]],]),combAB.ind(n,comb.groups(idx,n.mod)[[j]]))
          mat.cg <- comb.groups(idx,n.mod)[[j]]
        } else {
          comb.group <- rbind(colSums(Ind[id1[[i]],]),combAB.ind(n,comb.groups(c(idx,0),n.mod)[[j]]))
          mat.cg <- comb.groups(c(idx,0),n.mod)[[j]]
        }
        if (length(mat.cg[,1])==length(id1[[i]])){
          mat.cg <- cbind(mat.cg,id1[[i]])          
          mat.cg <- apply(mat.cg,2,sort)
          mat.cg <- mat.cg[,order(mat.cg[1,])]
        } else {
          a <- matrix(0,ncol=length(mat.cg[1,]),nrow=abs(length(id1[[i]])-length(mat.cg[,1])))
          if (length(id1[[i]])>=length(mat.cg[,1])){
            mat.cg <- rbind(mat.cg,a)
            mat.cg <- cbind(mat.cg,id1[[i]])          
            mat.cg <- apply(mat.cg,2,sort)
            mat.cg <- mat.cg[,order(mat.cg[1,])]
          } else {
            t.add <- c(id1[[i]],rep(0,abs(length(id1[[i]])-length(mat.cg[,1]))))
            mat.cg <- cbind(mat.cg,t.add)          
            mat.cg <- apply(mat.cg,2,sort)
            mat.cg <- mat.cg[,order(mat.cg[1,])]
          }
        }        
        c.ind <- rep(0,length(mat.c.list))
        c.ind.b <- 0
        if (length(mat.c.list)>=1) for (ic in 1:length(mat.c.list)){
          if (all(dim(mat.cg)==dim(mat.c.list[[ic]]))) c.ind[ic] <- all(mat.cg==mat.c.list[[ic]])
          if (c.ind[ic]==TRUE) c.ind.b <- TRUE
        }          
        mat.c.list <- c(mat.c.list,list(mat.cg))
        
        if (c.ind.b==0) IndALL2 <- rbind(IndALL2,comb.group)
        if (c.ind.b==0) n.mod4 <- c(n.mod4,rep(max(1,max(n.mod4)+1),dim(comb.group)[1]))
        
        c.ind.b <- 0
        #}
      }  
    }   
  }
  n.mod4 <- n.mod4[-1]
  return(list(IndALL2,n.mod4))
}


# Function to return the index matrix
ind.m.fun <- function(n){
  
  if (n < 2) print("wrong model: there should be at least two species/signals to compare")
  
  # If only two species or signals
  
  if (n == 2) {
    Ind <- diag(n)
    n.mod <- c(rep(1,n))
    n.mod <- c(n.mod, max(n.mod)+1)
  }
  
  # If more than 2 species or signals
  
  if (n >= 3) {
    I.AB <- rbind(diag(n), combAB.fun(n)[[1]])
    n.mod2 <- c(rep(1,n), combAB.fun(n)[[2]])
    if (n==3) {
      Ind <- I.AB
      n.mod <- c(n.mod2, max(n.mod2)+1)
    }
  }
  
  if (n >= 4) {
    I.AB_C_D <- rbind(I.AB, combAB_C_D.fun(n)[[1]])
    n.mod3 <- c(n.mod2, combAB_C_D.fun(n)[[2]]+max(n.mod2))
    if (n==4) {
      Ind <- I.AB_C_D
      n.mod <- c(n.mod3, max(n.mod3)+1)
    }
  }
  
  if (n >= 5) {
    Ind <- rbind(I.AB_C_D, all.comb(n)[[1]])
    n.mod <- c(n.mod3, all.comb(n)[[2]]+max(n.mod3))
    n.mod <- c(n.mod,max(n.mod)+1)
  }
  
  Ind <- rbind(Ind,rep(1,n))
  return(list(Ind,n.mod))
}


# Get species groupings
# Written by Phil Bouchet
get_groupings <- function(species.matrix, 
                          model.index, 
                          species.names, 
                          latin = TRUE, 
                          simulation){
  
  if(latin){
    for(i in 1:length(species.names)){
      if(species.names[i]== "killer") species.names[i] <- "Oo"
      if(species.names[i]== "lf_pilot") species.names[i] <- "Gm"
      if(species.names[i]== "sperm") species.names[i] <- "Pm"
      if(species.names[i]== "humpback") species.names[i] <- "Mn"
      if(species.names[i]== "minke") species.names[i] <- "Ba"
      if(species.names[i]== "blue") species.names[i] <- "Bm"
      if(species.names[i]== "nbottle") species.names[i] <- "Ha"
      if(species.names[i]== "blain") species.names[i] <- "Md"
      if(species.names[i]== "cuvier") species.names[i] <- "Zc"
      if(species.names[i]== "baird") species.names[i] <- "Bb"
    }
  }
  
  spn <- apply(X = species.matrix, MARGIN = 1, FUN = function(x) {
    species.names[which(x == 1)]
  }) 
  
  if(simulation) spn <- purrr::map(.x = spn, .f = ~gsub(pattern = "Sp", replacement = "", x = .x))
  
  spn <- purrr::map(.x = spn, .f = ~ {
    if (length(.x) > 1) paste0("(", paste(.x, collapse = ","), ")") else paste0("(", .x, ")")
  })
  
  purrr::map_chr(.x = unique(model.index[,1]),
                 .f = ~paste(spn[which(model.index[,1]==.x)], collapse = "+"))
}

#' Print method for objects of class \code{gvs}
#' @export
print.gvs <- function(gvs.dat){
 print(summary(gvs.dat$trace))
}

# Convenience ----------------------------------------------------------------

# sense_check <- function(print.text){
#   
#   out <- 
#     
#     c(sum(rj$t.ij[rj$iter["t.ij"], rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 2] < rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials][rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 2]),
#       
#       sum(rj$t.ij[rj$iter["t.ij"], rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 1] > rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials][rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 1]),
#       
#       sum(rj$t.ij[rj$iter["t.ij"], rj$k.ij[rj$iter["k.ij"], ] == 2] < rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials][rj$k.ij[rj$iter["k.ij"], ] == 2]),
#       
#       sum(rj$t.ij[rj$iter["t.ij"], rj$k.ij[rj$iter["k.ij"], ] == 1] > rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials][rj$k.ij[rj$iter["k.ij"], ] == 1]),
#       
#       sum(rj$mu.ij[rj$iter["mu.ij"], , 2] < rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials]),
#       sum(rj$mu.ij[rj$iter["mu.ij"], , 1] > rj$alpha[rj$iter["alpha"], rj.obj$dat$species$trials]),
#       
#       sum(rj$t.ij[rj$iter["t.ij"], rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 2] < rj.obj$dat$obs$Rc[rj$dat$obs$censored == 1 & rj$k.ij[rj$iter["k.ij"], ] == 2]),
#       
#       sum(rj$t.ij[rj$iter["t.ij"], rj$dat$obs$censored == -1 & rj$k.ij[rj$iter["k.ij"], ] == 1] > rj.obj$dat$obs$Lc[rj$dat$obs$censored == -1 & rj$k.ij[rj$iter["k.ij"], ] == 1]))
#   
#   cat(print.text, sum(out), '\n')
#   return(out)
# }

sense_check <- function(rj.obj, print.text){
  
  out <- 
    c(sum(rj.obj$mu.ij[rj.obj$iter["mu.ij"], , 1] > rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials]),
      sum(rj.obj$mu.ij[rj.obj$iter["mu.ij"], , 2] < rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials]),
      sum(rj.obj$nu[rj.obj$iter["nu"], , 1] > rj.obj$alpha[rj.obj$iter["alpha"],]),
      sum(rj.obj$nu[rj.obj$iter["nu"], , 2] < rj.obj$alpha[rj.obj$iter["alpha"],]),
      sum(rj.obj$t.ij[rj.obj$iter["t.ij"], rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1] > rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials][rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 1]),
      sum(rj.obj$t.ij[rj.obj$iter["t.ij"], rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2] < rj.obj$alpha[rj.obj$iter["alpha"], rj.obj$dat$species$trials][rj.obj$k.ij[rj.obj$iter["k.ij"], ] == 2]),
      sum(rj.obj$t.ij[rj.obj$iter["t.ij"], rj.obj$dat$obs$censored == -1] > rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]),
      sum(rj.obj$t.ij[rj.obj$iter["t.ij"], rj.obj$dat$obs$censored == 1] < rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]))
    
  # cat(print.text, sum(out), '\n')
  return(out)
}


#' Transmission loss

#' @param r Range in km.
#' @param a Sound absorption coefficient, in dB per km. This is frequency-dependent, and takes a value of 0.185 for a 3 kHz signal under normal sea conditions.

TL <- function(rge, a = 0.185){ 
  loss <- 20*log10(rge*1000)
  loss[loss<0] <- 0
  loss <- loss+a*rge
  return(loss)}

# Function to calculate the range corresponding to a given RL
# Return squared difference between target.TL and actual TL

#' Range finder

#' @param rge Range in km.
#' @param SL Level of the noise source.
#' @param target.L Target noise level.

range_finder <- function(rge, SL, target.L){
  return((SL-TL(rge)-target.L)^2)}

#' Mimick transparent colours
#' 
#' Function to find the HEX colour code corresponding to an input colour 
#' with a set opacity level (i.e. emulate transparency)

#' @param input.colour Initial colour.
#' @param opacity Desired level of transparency (number between 0 and 1).
#' @param bg.colour Colour of the background. Defaults to 'white'.

hexa2hex <- function(input.colour, 
                     opacity, 
                     bg.colour = "white"){
  
  # White background
  bg <- grDevices::col2rgb(bg.colour, alpha = FALSE)
  
  # Convert input colour to RGB
  rgbcol <- grDevices::col2rgb(input.colour, alpha = FALSE)
  
  # Calculate red, green, blue values corresponding to input colour at chosen transparency level
  rc <- (1 - opacity) * bg[1,] + opacity * rgbcol[1,]
  gc <- (1 - opacity) * bg[2,] + opacity * rgbcol[2,]
  bc <- (1 - opacity) * bg[3,] + opacity * rgbcol[3,]
  
  # Convert back to hex
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  return(rgb2hex(r = rc, g = gc, b = bc))
}

# Covariate effects (on probit scale, biphasic model)
cov_effects <- function(rj.obj){
  
  included.cov <- rj.obj$include.covariates[rj.obj$iter["covariates"], ]
  
  if (sum(included.cov) == 0) {
    covnames <- "nil"
    cov.effects <- 0
  } else {
    included.cov <- included.cov[included.cov == 1]
    covvalues <- do.call(c, purrr::map(
      .x = names(included.cov),
      .f = ~ {qnorm(pnorm(q = rj.obj[[.x]][rj.obj$iter[.x], ], 
                          mean = rj.obj$config$priors[.x, 1], 
                          sd = rj.obj$config$priors[.x, 2]))}
    ))
    cov.effects <- colSums(t(do.call(cbind, rj.obj$dat$covariates$dummy[names(included.cov)])) * covvalues)
  }
  return(cov.effects)
}

# Rescale probability vector
rescale_p <- function(p, default.value = 0.05){
  if(length(p) == 1){
    p.out <- 1
  } else {
  min.p <- min(p[p > 0])
  if(min.p == 1) min.p <- default.value
  if(length(p) > 1) p[p == 0] <- min.p
  p.out <- p / sum(p)
  }
  return(p.out)
}

# Posterior model probabilities
prob_models <- function(input.obj, 
                        n.top = NULL,
                        mlist = NULL, 
                        select = NULL, 
                        do.combine, 
                        gvs){
  
  input.trace <- input.obj$trace
  if(gvs) species.Groups <- input.obj$species.Groups else species.Groups <- NULL
  
  if(do.combine){
    ptrace <- do.call(rbind, input.trace) %>% coda::as.mcmc()
    p_models(input.trace = ptrace,
             n.top = n.top,
             mlist = mlist, 
             select = select, 
             gvs = gvs, 
             species.Groups = species.Groups)
  } else {
    purrr::map(.x = input.trace, 
               .f = ~p_models(input.trace = .x, 
                              n.top = n.top,
                              mlist = mlist, 
                              select = select, 
                              gvs = gvs, 
                              species.Groups = species.Groups))
  }
}

p_form <- function(input.trace){

  res <- list()
  
  # Model trace
  phase.trace <- input.trace[, "phase"]
  
  # Tabulate
  res$est <- table(phase.trace)/nrow(input.trace)
  
  # Generate tibble output
  res$est <- res$est %>% 
    tibble::enframe() %>% 
    dplyr::rename(phase = name, p = value) %>% 
    dplyr::mutate(phase = dplyr::case_when(phase == "2" ~ "biphasic",
                                   TRUE ~ "monophasic")) %>%
    dplyr::mutate(p = as.numeric(p))
  
  # Compute Monte Carlo error
  mce <- purrr::map(.x = 1:2, .f = ~as.numeric(as.numeric(phase.trace) == .x))
  
  res$mc.error <- purrr::map_dbl(.x = mce, .f = ~mcmcse::mcse(.x)$se) %>% 
    purrr::set_names(x = ., nm = c("monophasic", "biphasic")) %>% 
    tibble::enframe() %>% 
    dplyr::rename(phase = name, monte_carlo = value)
  
  # Compute lower and upper interval bounds
  res$est <- res$est %>% 
    dplyr::left_join(x = ., y = res$mc.error, by = "phase") %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(lower = max(0, round(p - 1.98 * monte_carlo, 4)),
                  upper = min(1, round(p + 1.98 * monte_carlo, 4))) %>% 
    dplyr::ungroup()
  
  return(res)
  
}

p_models <- function(input.trace, 
                     n.top,
                     mlist, 
                     select, 
                     gvs, 
                     species.Groups){
  
  res <- list()
  
  if(gvs){
    
    # Extract posterior samples

    mcmc.values <- as.matrix(input.trace) %>%
      tibble::as_tibble() %>%
      janitor::clean_names(.)
    
    # Calculate posterior model probabilities
    mselect.species <- table(mcmc.values$theta)/nrow(mcmc.values)
    model.ID.posterior <- as.numeric(names(sort(mselect.species, decreasing = TRUE)))
    names(mselect.species) <- species.Groups[as.numeric(names(mselect.species))]
    mselect.species <- sort(mselect.species, decreasing = TRUE)
    m_prob <- tibble::tibble(as.data.frame(mselect.species)) %>% 
      dplyr::arrange(-value) %>% 
      dplyr::slice(1:n.top)
    
    if(ncol(m_prob) == 1) m_prob <- m_prob %>%
      dplyr::rename(p = mselect.species) %>%
      dplyr::mutate(model = names(mselect.species)) %>%
      dplyr::select(model, p) else m_prob <- m_prob %>%
      dplyr::rename(model = Var1, p = Freq)
    
    # Monte Carlo error
    mtrace <- mcmc.values[["theta"]]
    
    # Compute Monte Carlo error
    mce <- sapply(X = unique(mtrace), FUN = function(x) as.numeric(mtrace == x))
    mc.error <- apply(X = mce, MARGIN = 2, FUN = mcmcse::mcse) %>% 
      purrr::map_dbl(.x = ., .f = "se")
    
    mc.error <- tibble::as_tibble(mc.error, rownames = "model") %>% 
      dplyr::mutate(model = m_prob$model) %>% 
      dplyr::rename(monte_carlo = value)
    
    # Compute lower and upper interval bounds
    m_prob <- m_prob %>% 
      dplyr::left_join(x = ., y = mc.error, by = "model") %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(lower = max(0, round(p - 1.98 * monte_carlo[dplyr::row_number()], 4)),
                    upper = min(1, round(p + 1.98 * monte_carlo[dplyr::row_number()], 4))) %>% 
      dplyr::ungroup()
    
    bestmod <- as.character(m_prob$model[(which.max(m_prob$p))])
    
    res$model$m_prob <- m_prob
    res$model$bestmod <- bestmod
    res$mc.error <- mc.error
    
    
  } else {
    
    # Model trace
    mtrace <- input.trace[, "model_ID"]
    
    # Tabulate
    res$model$m_prob <- table(mtrace)/nrow(input.trace)
    
    # Generate tibble output
    res$model$m_prob <- res$model$m_prob %>% 
      tibble::enframe() %>% 
      dplyr::arrange(-value) %>% 
      dplyr::rename(ID = name, p = value) %>% 
      dplyr::mutate(ID = as.numeric(ID)) %>% 
      dplyr::left_join(., mlist[, c("ID", "model")], by = "ID") %>% 
      dplyr::mutate(p = as.numeric(p)) %>% 
      dplyr::slice(1:n.top)
    
    # Identify top-ranking model
    res$model$bestmod <- res$model$m_prob$model[1]
    
    # Compute Monte Carlo error
    mce <- as.numeric(mtrace)
    mce <- purrr::map(.x = res$model$m_prob$ID, .f = ~as.numeric(mce == .x))

    mce <- lapply(X = res$model$m_prob$ID, FUN = function(x) as.numeric(mtrace == x))
    
    # There are some differences in the ESS values calculated by coda and mcmcse
    # The original code (using mean * 1 - mean etc.) gives exactly the same Monte Carlo errors
    # as the mcse function from <mcmcse> when neff is calculated using the mcmcse package,
    # neff <- apply(X = mce, MARGIN = 2, FUN = coda::effectiveSize)
    # neff <- apply(X = mce, MARGIN = 2, FUN = mcmcse::ess)
    
    res$mc.error <- purrr::map_dbl(.x = mce, .f = ~mcmcse::mcse(.x)$se) %>% 
      purrr::set_names(x = ., nm = res$model$m_prob$model) %>% 
      tibble::enframe() %>% 
      dplyr::rename(model = name, monte_carlo = value)
    
    # res$mc.error <- sapply(X = seq_along(unique(mtrace)),
    #                        FUN = function(x) sqrt(mean(mce[, x]) * (1 - mean(mce[, x])) / neff[x]))
    
    # names(res$mc.error) <- purrr::map_chr(.x = unique(mtrace), 
    #                                       .f = ~mlist %>% 
    #                                         dplyr::filter(ID == .x) %>% 
    #                                         dplyr::pull(model))
    
    # res$mc.error <- tibble::as_tibble(res$mc.error, rownames = "model") %>% 
    #   dplyr::rename(monte_carlo = value)
    
    # Compute lower and upper interval bounds
    res$model$m_prob <- res$model$m_prob %>% 
      dplyr::left_join(x = ., y = res$mc.error, by = "model") %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(lower = max(0, round(p - 1.98 * monte_carlo, 4)),
                    upper = min(1, round(p + 1.98 * monte_carlo, 4))) %>% 
      dplyr::ungroup()
    
  }
  
  return(res)
  
}

# Function to create a tiled representation of model rankings using ggplot
gg_model <- function(dat, 
                     post.probs, 
                     rj.obj, 
                     colours, 
                     n.top, 
                     combine, 
                     no = 0, 
                     p.size = 4, 
                     x.offset = 0.5, 
                     x.margin = 2){
  
  rotate.x <- min(nchar(unique(dat$species))) > 5
  
  # Probability values
  if(combine){
    p1 <- tibble::tibble(x = rep(max(dat$x) + x.offset, n.top),
                         y = unique(dat$y), 
                         label = post.probs)
  } else {
    p1 <- tibble::tibble(x = rep(max(dat$x) + x.offset, length(unique(dat$y))), 
                         y = unique(dat$y), 
                         label = post.probs)
  }
  
  # Extra blank space
  p2 <- p1 %>% dplyr::mutate(x = x + x.margin, label = "")
  
  out.plot <- ggplot2::ggplot(data = dat, aes(x = x, y = y)) + 
    ggplot2::geom_tile(aes(fill = as.factor(grouping)), col = "white", size = 0.25) +
    ggplot2::scale_fill_manual(values = colours) + 
    ggplot2::xlab("") +
    {if(!combine) ggplot2::ylab("Chain")} +
    {if(combine) ggplot2::ylab("Rank")} +
    ggplot2::scale_y_continuous(breaks = rev(seq_len(n.top)), 
                                labels = 1:n.top,
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous(breaks = seq_len(rj.obj$dat$species$n), 
                                labels = rj.obj$dat$species$names, 
                                expand = c(0, 0)) +
    ggplot2::theme(axis.text = element_text(size = 12, colour = "black"),
                   axis.title = element_text(size = 12, colour = "black"),
                   axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                   axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 00, l = 0), 
                                               angle = 90, vjust = 0.5, hjust = 0.5),
                   plot.margin = margin(t = 1, r = 1, b = 0.25, l = 1, "cm"),
                   legend.position = "top",
                   # legend.title = element_blank(),
                   legend.text = element_text(size = 12),
                   panel.background = element_rect(fill = "white")) + 
    {if(rotate.x) ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))} +
    
    # Annotation for posterior probabilities
    ggplot2::geom_text(data = p1, aes(x = x , y = y, label = label), size = p.size, hjust = 0) +
    ggplot2::geom_text(data = p2, aes(x = x , y = y, label = label), size = p.size) +
    
    labs(fill = "Groups") +
    
    ggplot2::guides(fill = guide_legend(nrow = 1)) +
    
    {if(!combine) ggplot2::theme(legend.position = "none")} +
    {if(!combine) ggtitle(paste0("Rank: ", no)) }
  
  return(out.plot)
}

prob_form <- function(obj, do.combine){
  
  if(do.combine){
    ptrace <- do.call(rbind, obj$trace) %>% coda::as.mcmc()
    p_form(input.trace = ptrace)
  } else {
    purrr::map(.x = obj$trace, .f = ~p_form(input.trace = .x))
  }
}

#' Posterior inclusion probabilities for contextual covariates
prob_covariates <- function(obj, do.combine){

  if(do.combine){
    ptrace <- do.call(rbind, obj$trace) %>% coda::as.mcmc()
    p_covariates(input.trace = ptrace, ddf = obj$dat)
  } else {
    purrr::map(.x = obj$trace, .f = ~p_covariates(input.trace = .x, ddf = obj$dat))
  }
}

p_covariates <- function(input.trace, ddf){
  
  # Identify the covariate columns
  col.indices <- which(grepl(pattern = "incl.", x = colnames(input.trace)))
  
  # Extract the data and update column names
  covtrace <- input.trace[, col.indices, drop = FALSE] %>% 
    tibble::as_tibble(.)
  colnames(covtrace) <- gsub(x = colnames(covtrace), pattern = "incl.", replacement = "")
  
  # Create a duplicate tibble used to generate labels for each combination of covariates
  covtrace.labels <- covtrace
  for(j in ddf$covariates$names) covtrace.labels[, j] <- ifelse(covtrace.labels[, j] == 0, NA, j)
  covtrace.labels <- covtrace.labels %>% 
    tidyr::unite("comb", ddf$covariates$names, na.rm = TRUE, remove = TRUE, sep = " + ")
  
  if(ddf$covariates$n > 1){
  covtrace.labels$comb <- ifelse(stringr::str_detect(pattern = "\\+", string = covtrace.labels$comb),
                                 covtrace.labels$comb, paste0(covtrace.labels$comb, " (only)"))
  covtrace.labels$comb <- ifelse(covtrace.labels$comb == " (only)", "No covariates", covtrace.labels$comb)
  }
  
  # Dummy coding for each covariate combination
  dummy.covariates <- fastDummies::dummy_cols(covtrace.labels)
  names(dummy.covariates) <- gsub(x = names(dummy.covariates), pattern = "comb_", replacement = "")
  
  # Calculate overall posterior probabilities for each covariate (i.e., inclusive of all combinations
  # in which the covariate appears)
  tab.overall <- purrr::map_dbl(.x = ddf$covariates$names, 
                                .f = ~{
                                  tmp <- table(covtrace[, .x]) / nrow(covtrace)
                                  tmp[names(tmp)=="1"]}) %>% purrr::set_names(x = ., nm = ddf$covariates$names) %>% 
    tibble::enframe(.) %>% 
    dplyr::rename(covariate = name, p = value) %>% 
    dplyr::rowwise()  %>% 
    dplyr::mutate(monte_carlo = mcmcse::mcse(dplyr::pull(covtrace[, covariate]))$se) %>% 
    dplyr::mutate(lower = max(0, round(p - 1.98 * monte_carlo, 4)),
                  upper = min(1, round(p + 1.98 * monte_carlo, 4))) %>% 
    dplyr::ungroup()
  
  if(ddf$covariates$n > 1){
    tab.overall <- tab.overall %>% 
      dplyr::mutate(covariate = paste0(covariate, " (overall)"))
  }
  
  if(ddf$covariates$n > 1){
    
    # Calculate posterior probabilities for specific covariate combinations
    tab.combinations <- table(covtrace.labels) / nrow(covtrace.labels)
    tab.combinations <- tibble::as_tibble(tab.combinations) %>% 
      dplyr::rename(covariate = covtrace.labels, p = n) %>% 
      dplyr::filter(!covariate %in% c("No covariates", "")) %>% 
      dplyr::arrange(covariate) %>% 
      dplyr::rowwise()  %>% 
      dplyr::mutate(monte_carlo = mcmcse::mcse(dplyr::pull(dummy.covariates[, covariate]))$se) %>% 
      dplyr::mutate(lower = max(0, round(p - 1.98 * monte_carlo, 4)),
                    upper = min(1, round(p + 1.98 * monte_carlo, 4))) %>% 
      dplyr::mutate(ind_order = ifelse(stringr::str_detect(string = covariate, pattern = "(only)"), 1, 99)) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(ind_order, covariate) %>% 
      dplyr::select(-ind_order)
    
    tb <- dplyr::bind_rows(tab.overall, tab.combinations)
    
  } else {
    
    tb <- tab.overall
  }
  
  return(tb)
  
}

gg_color_hue <- function(n) {
  
  # Custom palette - to match <bayesplot>
  if(n == 1) col <-"#6497b1"
  if(n == 2) col <-c("#6497b1", "#005b96")
  if(n == 3) col <-c("#b3cde0", "#6497b1", "#005b96")
  if(n > 3) col <- c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")[1:n]
  
  # Original palette
  hues <- seq(50, 300, length = n + 1)
  col <- grDevices::hcl(h = hues, l = 70, c = 100)[1:n]
  return(col)
}

MCMC_trace <- function (object, params = "all", excl = NULL, ISB = TRUE, iter = 5000, 
            gvals = NULL, priors = NULL, post_zm = TRUE, PPO_out = FALSE, 
            Rhat = FALSE, n.eff = FALSE, ind = FALSE, pdf = TRUE, plot = TRUE, 
            open_pdf = TRUE, filename, wd = getwd(), type = "both", ylim = NULL, 
            xlim = NULL, xlab_tr, ylab_tr, xlab_den, ylab_den, main_den = NULL, 
            main_tr = NULL, lwd_den = 1, lwd_pr = 1, lty_den = 1, lty_pr = 1, 
            col_den, col_pr, col_txt, sz_txt = 1.2, sz_ax = 1, sz_ax_txt = 1, 
            sz_tick_txt = 1, sz_main_txt = 1.2, pos_tick_x_tr = NULL, 
            pos_tick_y_tr = NULL, pos_tick_x_den = NULL, pos_tick_y_den = NULL) 
  {
  
  
  # Adapted from MCMCvis::MCMCtrace
  # Only change made is to gg_color_hue to allow for different colour palette
  
    .pardefault <- graphics::par(no.readonly = T)
    if (methods::is(object, "matrix")) {
      warning("Input type matrix - assuming only one chain for each parameter.")
      object1 <- coda::as.mcmc.list(coda::as.mcmc(object))
      object2 <- MCMCvis::MCMCchains(object1, params, excl, ISB, mcmc.list = TRUE)
    }
    else {
      object2 <- MCMCvis::MCMCchains(object, params, excl, ISB, mcmc.list = TRUE)
    }
    np <- colnames(object2[[1]])
    n_chains <- length(object2)
    if (nrow(object2[[1]]) > iter) {
      it <- (nrow(object2[[1]]) - iter + 1):nrow(object2[[1]])
    }
    else {
      it <- 1:nrow(object2[[1]])
    }
    if (!is.null(priors)) {
      if (NCOL(priors) == 1 & length(np) > 1) {
        warning("Only one prior specified for > 1 parameter. Using a single prior for all parameters.")
      }
      if ((NCOL(priors) > 1 & NCOL(priors) != length(np))) {
        stop("Number of priors does not equal number of specified parameters.")
      }
      if (NROW(priors) > length(it) * n_chains) {
        warning(paste0("Number of samples in prior is greater than number of total or specified iterations (for all chains) for specified parameter. Only last ", 
                       length(it) * n_chains, " iterations will be used."))
      }
      if (NROW(priors) < length(it) * n_chains) {
        warning(paste0("Number of samples in prior is less than number of total or specified iterations (for all chains) for specified parameter. Resampling from prior to generate ", 
                       length(it) * n_chains, " total iterations."))
      }
      if (type == "trace") {
        warning("Prior posterior overlap (PPO) cannot be plotting without density plots. Use type = 'both' or type = 'density'.")
      }
    }
    if (!is.null(gvals)) {
      if (length(gvals) == 1 & length(np) > 1) {
        warning("Only one generating value specified for > 1 parameter. Using a single generating value for all parameters.")
      }
      if (length(gvals) > 1 & length(gvals) != length(np)) {
        stop("Number of generating values does not equal number of specified parameters.")
      }
    }
    if (PPO_out == TRUE) {
      PPO_df <- data.frame(param = rep(NA, length(np)), percent_PPO = rep(NA, length(np)))
    }
    if (plot == TRUE) {
      if (pdf == TRUE) {
        if (missing(filename)) {
          file_out <- paste0(wd, "/MCMCtrace.pdf")
        }
        else {
          if (grepl(".pdf", filename, fixed = TRUE)) {
            file_out <- paste0(wd, "/", filename)
          }
          else {
            file_out <- paste0(wd, "/", filename, ".pdf")
          }
        }
        pdf(file = file_out)
      }
      ref_col <- "red"
      A_VAL <- 0.5
      if (type == "both") {
        if (length(np) >= 3) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  3, 2, byrow = TRUE))
          graphics::par(mar = c(4.1, 4.1, 2.1, 1.1))
          MN_LINE <- NULL
        }
        if (length(np) == 2) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  2, 2, byrow = TRUE))
          graphics::par(mar = c(4.1, 4.1, 2.1, 1.1))
          MN_LINE <- NULL
        }
        if (length(np) == 1) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  1, 2, byrow = TRUE))
          graphics::par(mar = c(8.1, 4.1, 7.1, 1.1))
          MN_LINE <- 1.1
        }
      }
      else {
        if (length(np) >= 5) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  3, 2, byrow = TRUE))
          graphics::par(mar = c(4.1, 4.1, 2.1, 1.1))
          MN_LINE <- NULL
        }
        if (length(np) == 3 | length(np) == 4) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  2, 2, byrow = TRUE))
          graphics::par(mar = c(4.1, 4.1, 2.1, 1.1))
          MN_LINE <- NULL
        }
        if (length(np) == 2) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  1, 2, byrow = TRUE))
          graphics::par(mar = c(8.1, 4.1, 7.1, 1.1))
          MN_LINE <- 1.1
        }
        if (length(np) == 1) {
          graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 
                                  1, 1, byrow = TRUE))
          graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
          MN_LINE <- NULL
        }
      }
      graphics::par(mgp = c(2.5, 1, 0))
      
      colors <- gg_color_hue(n_chains)
      gg_cols <- grDevices::col2rgb(colors)/255
      YLIM <- ylim
      XLIM <- xlim
      if (!missing(xlab_tr)) {
        xlab_tr <- xlab_tr
      }
      else {
        xlab_tr <- "Iteration"
      }
      if (!missing(ylab_tr)) {
        ylab_tr <- ylab_tr
      }
      else {
        ylab_tr <- "Value"
      }
      if (!missing(xlab_den)) {
        xlab_den <- xlab_den
      }
      else {
        xlab_den <- "Parameter estimate"
      }
      if (!missing(ylab_den)) {
        ylab_den <- ylab_den
      }
      else {
        ylab_den <- "Density"
      }
      if (!missing(main_den)) {
        if (length(main_den) != 1 & length(main_den) != length(np)) {
          stop("Number of elements for 'main_den' does not equal number of specified parameters.")
        }
      }
      if (!missing(main_tr)) {
        if (length(main_tr) != 1 & length(main_tr) != length(np)) {
          stop("Number of elements for 'main_tr' does not equal number of specified parameters.")
        }
      }
      MAIN_DEN <- function(x, md = main_den, idx) {
        if (!is.null(md)) {
          if (length(md) == 1) {
            return(md)
          }
          else {
            return(md[idx])
          }
        }
        else {
          # return(paste0("Density - ", x))
          return(paste0(x))
        }
      }
      MAIN_TR <- function(x, mtr = main_tr, idx) {
        if (!is.null(mtr)) {
          if (length(mtr) == 1) {
            return(mtr)
          }
          else {
            return(mtr[idx])
          }
        }
        else {
          # return(paste0("Trace - ", x))
          return(paste0(x))
        }
      }
      if (missing(col_den)) {
        COL_DEN <- "black"
      }
      else {
        COL_DEN <- col_den
      }
      if (missing(col_pr)) {
        COL_PR <- "red"
      }
      else {
        COL_PR <- col_pr
      }
      if (missing(col_txt)) {
        COL_TXT <- "red"
      }
      else {
        COL_TXT <- col_txt
      }
      if (is.null(pos_tick_x_den)) {
        XAXT_DEN <- "s"
      }
      else {
        XAXT_DEN <- "n"
      }
      if (is.null(pos_tick_y_den)) {
        YAXT_DEN <- "s"
      }
      else {
        YAXT_DEN <- "n"
      }
      if (is.null(pos_tick_x_tr)) {
        XAXT_TR <- "s"
      }
      else {
        XAXT_TR <- "n"
      }
      if (is.null(pos_tick_y_tr)) {
        YAXT_TR <- "s"
      }
      else {
        YAXT_TR <- "n"
      }
      if (Rhat == TRUE | n.eff == TRUE) {
        summ <- MCMCsummary(object, params = params, excl = excl, 
                            ISB = ISB, Rhat = Rhat, n.eff = n.eff)
        if (Rhat == TRUE) {
          rhat <- summ[, grep("Rhat", colnames(summ))]
        }
        if (n.eff == TRUE) {
          neff <- summ[, grep("n.eff", colnames(summ))]
        }
      }
      if (type == "both") {
        for (j in 1:length(np)) {
          tmlt <- do.call("cbind", object2[it, np[j]])
          graphics::matplot(it, tmlt, lwd = 1, lty = 1, 
                            type = "l", main = NULL, col = grDevices::rgb(red = gg_cols[1, 
                            ], green = gg_cols[2, ], blue = gg_cols[3, 
                            ], alpha = A_VAL), xlab = xlab_tr, ylab = ylab_tr, 
                            cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, 
                            xaxt = XAXT_TR, yaxt = YAXT_TR)
          graphics::title(main = MAIN_TR(np[j], main_tr, 
                                         j), line = MN_LINE, cex.main = sz_main_txt)
          graphics::axis(side = 1, at = pos_tick_x_tr, 
                         cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_tr, 
                         cex.axis = sz_tick_txt)
          if (!missing(sz_ax)) {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                           at = pos_tick_x_tr)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                           at = pos_tick_y_tr)
          }
          if (!is.null(priors)) {
            if (NCOL(priors) == 1) {
              wp <- priors
            }
            else {
              wp <- priors[, j]
            }
            lwp <- length(wp)
            if (lwp > length(it) * n_chains) {
              pit <- (lwp - (length(it) * n_chains) + 1):lwp
              wp2 <- wp[pit]
            }
            if (lwp < length(it) * n_chains) {
              samps <- sample(wp, size = ((length(it) * 
                                             n_chains) - lwp), replace = TRUE)
              wp2 <- c(wp, samps)
            }
            if (lwp == length(it) * n_chains) {
              wp2 <- wp
            }
            dpr <- stats::density(wp2)
            PPO_x_rng <- range(dpr$x)
            PPO_y_rng <- range(dpr$y)
            tmlt_1c <- matrix(tmlt, ncol = 1)
            pp <- list(wp2, tmlt_1c)
            ovr_v <- round((overlapping::overlap(pp)$OV[[1]]) * 
                             100, digits = 1)
            ovrlap <- paste0(ovr_v, "% overlap")
            if (PPO_out == TRUE) {
              PPO_df$param[j] <- np[j]
              PPO_df$percent_PPO[j] <- ovr_v
            }
          }
          if (ind == TRUE & n_chains > 1) {
            dens <- apply(tmlt, 2, stats::density)
            max_den_y <- c()
            rng_den_x <- c()
            for (k in 1:NCOL(tmlt)) {
              max_den_y <- c(max_den_y, max(dens[[k]]$y))
              rng_den_x <- c(rng_den_x, range(dens[[k]]$x))
            }
            if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & 
                is.null(xlim)) {
              ylim <- range(c(0, max(max_den_y), PPO_y_rng))
              xlim <- range(c(range(rng_den_x), PPO_x_rng))
            }
            else {
              if (!is.null(ylim)) {
                ylim <- YLIM
                xlim <- XLIM
              }
              if (is.null(ylim) & is.null(xlim)) {
                ylim <- c(0, max(max_den_y))
                xlim <- NULL
              }
            }
            graphics::plot(dens[[1]], xlab = xlab_den, 
                           ylab = ylab_den, ylim = ylim, xlim = xlim, 
                           lty = lty_den, lwd = lwd_den, main = "", 
                           col = grDevices::rgb(red = gg_cols[1, 1], 
                                                green = gg_cols[2, 1], blue = gg_cols[3, 
                                                                                      1]), cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, 
                           xaxt = XAXT_DEN, yaxt = YAXT_DEN)
            graphics::title(main = MAIN_DEN(np[j], main_den, 
                                            j), line = MN_LINE, cex.main = sz_main_txt)
            graphics::axis(side = 1, at = pos_tick_x_den, 
                           cex.axis = sz_tick_txt)
            graphics::axis(side = 2, at = pos_tick_y_den, 
                           cex.axis = sz_tick_txt)
            for (l in 2:NCOL(tmlt)) {
              graphics::lines(dens[[l]], lty = lty_den, 
                              lwd = lwd_den, col = grDevices::rgb(red = gg_cols[1, 
                                                                                l], green = gg_cols[2, l], blue = gg_cols[3, 
                                                                                                                          l]))
            }
            if (!missing(sz_ax)) {
              graphics::box(lwd = sz_ax)
              graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_x_den)
              graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_y_den)
            }
          }
          else {
            dens <- stats::density(rbind(tmlt))
            rng_den_x <- range(dens$x)
            if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & 
                is.null(xlim)) {
              ylim <- range(c(range(dens$y), PPO_y_rng))
              xlim <- range(c(range(dens$x), PPO_x_rng))
            }
            else {
              if (!is.null(ylim)) {
                ylim <- YLIM
                xlim <- XLIM
              }
              if (is.null(ylim) & is.null(xlim)) {
                ylim <- NULL
                xlim <- NULL
              }
            }
            graphics::plot(dens, xlab = xlab_den, ylab = ylab_den, 
                           ylim = ylim, main = "", col = COL_DEN, xlim = xlim, 
                           lty = lty_den, lwd = lwd_den, cex.axis = sz_tick_txt, 
                           cex.lab = sz_ax_txt, xaxt = XAXT_DEN, yaxt = YAXT_DEN)
            graphics::title(main = MAIN_DEN(np[j], main_den, 
                                            j), line = MN_LINE, cex.main = sz_main_txt)
            graphics::axis(side = 1, at = pos_tick_x_den, 
                           cex.axis = sz_tick_txt)
            graphics::axis(side = 2, at = pos_tick_y_den, 
                           cex.axis = sz_tick_txt)
            if (!missing(sz_ax)) {
              graphics::box(lwd = sz_ax)
              graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_x_den)
              graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_y_den)
            }
          }
          if (!is.null(priors)) {
            graphics::lines(dpr, col = COL_PR, lwd = lwd_pr, 
                            lty = lty_pr)
            if (!is.null(sz_txt) & !is.null(COL_TXT)) {
              graphics::legend("topright", legend = ovrlap, 
                               bty = "n", pch = NA, text.col = COL_TXT, 
                               cex = sz_txt)
            }
          }
          if (Rhat == TRUE & n.eff == TRUE) {
            diag_txt <- list(paste0("Rhat: ", rhat[j]), 
                             paste0("n.eff: ", neff[j]))
          }
          if (Rhat == TRUE & n.eff == FALSE) {
            diag_txt <- paste0("Rhat: ", rhat[j])
          }
          if (Rhat == FALSE & n.eff == TRUE) {
            diag_txt <- paste0("n.eff: ", neff[j])
          }
          if (Rhat == TRUE | n.eff == TRUE) {
            if (!is.null(sz_txt) & !is.null(COL_TXT)) {
              graphics::legend("topleft", x.intersp = -0.5, 
                               legend = diag_txt, bty = "n", pch = NA, 
                               text.col = COL_TXT, cex = sz_txt)
            }
          }
          if (!is.null(gvals)) {
            if (length(gvals) == 1) {
              gv <- gvals
            }
            else {
              gv <- gvals[j]
            }
            graphics::abline(v = gv, lty = 2, lwd = 3, 
                             col = ref_col)
          }
        }
      }
      if (type == "trace") {
        for (j in 1:length(np)) {
          tmlt <- do.call("cbind", object2[it, np[j]])
          graphics::matplot(it, tmlt, lwd = 1, lty = 1, 
                            type = "l", main = NULL, col = grDevices::rgb(red = gg_cols[1, 
                            ], green = gg_cols[2, ], blue = gg_cols[3, 
                            ], alpha = A_VAL), xlab = xlab_tr, ylab = ylab_tr, 
                            cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, 
                            xaxt = XAXT_TR, yaxt = YAXT_TR)
          graphics::title(main = MAIN_TR(np[j], main_tr, 
                                         j), line = MN_LINE, cex.main = sz_main_txt)
          graphics::axis(side = 1, at = pos_tick_x_tr, 
                         cex.axis = sz_tick_txt)
          graphics::axis(side = 2, at = pos_tick_y_tr, 
                         cex.axis = sz_tick_txt)
          if (!missing(sz_ax)) {
            graphics::box(lwd = sz_ax)
            graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                           at = pos_tick_x_tr)
            graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                           at = pos_tick_y_tr)
          }
        }
      }
      if (type == "density") {
        for (j in 1:length(np)) {
          tmlt <- do.call("cbind", object2[it, np[j]])
          if (!is.null(priors)) {
            if (NCOL(priors) == 1) {
              wp <- priors
            }
            else {
              wp <- priors[, j]
            }
            lwp <- length(wp)
            if (lwp > length(it) * n_chains) {
              pit <- (lwp - (length(it) * n_chains) + 1):lwp
              wp2 <- wp[pit]
            }
            if (lwp < length(it) * n_chains) {
              samps <- sample(wp, size = ((length(it) * 
                                             n_chains) - lwp), replace = TRUE)
              wp2 <- c(wp, samps)
            }
            if (lwp == length(it) * n_chains) {
              wp2 <- wp
            }
            dpr <- stats::density(wp2)
            PPO_x_rng <- range(dpr$x)
            PPO_y_rng <- range(dpr$y)
            tmlt_1c <- matrix(tmlt, ncol = 1)
            pp <- list(wp2, tmlt_1c)
            ovr_v <- round((overlapping::overlap(pp)$OV[[1]]) * 
                             100, digits = 1)
            ovrlap <- paste0(ovr_v, "% overlap")
            if (PPO_out == TRUE) {
              PPO_df$param[j] <- np[j]
              PPO_df$percent_PPO[j] <- ovr_v
            }
          }
          if (ind == TRUE & n_chains > 1) {
            dens <- apply(tmlt, 2, stats::density)
            max_den_y <- c()
            rng_den_x <- c()
            for (k in 1:NCOL(tmlt)) {
              max_den_y <- c(max_den_y, max(dens[[k]]$y))
              rng_den_x <- c(rng_den_x, range(dens[[k]]$x))
            }
            if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & 
                is.null(xlim)) {
              ylim <- range(c(0, max(max_den_y), PPO_y_rng))
              xlim <- range(c(range(rng_den_x), PPO_x_rng))
            }
            else {
              if (!is.null(ylim)) {
                ylim <- YLIM
                xlim <- XLIM
              }
              if (is.null(ylim) & is.null(xlim)) {
                ylim <- c(0, max(max_den_y))
                xlim <- NULL
              }
            }
            graphics::plot(dens[[1]], xlab = xlab_den, 
                           ylab = ylab_den, ylim = ylim, xlim = xlim, 
                           lty = lty_den, lwd = lwd_den, main = "", 
                           col = grDevices::rgb(red = gg_cols[1, 1], 
                                                green = gg_cols[2, 1], blue = gg_cols[3, 
                                                                                      1]), cex.axis = sz_tick_txt, cex.lab = sz_ax_txt, 
                           xaxt = XAXT_DEN, yaxt = YAXT_DEN)
            graphics::title(main = MAIN_DEN(np[j], main_den, 
                                            j), line = MN_LINE, cex.main = sz_main_txt)
            graphics::axis(side = 1, at = pos_tick_x_den, 
                           cex.axis = sz_tick_txt)
            graphics::axis(side = 2, at = pos_tick_y_den, 
                           cex.axis = sz_tick_txt)
            for (l in 2:NCOL(tmlt)) {
              graphics::lines(dens[[l]], lty = lty_den, 
                              lwd = lwd_den, col = grDevices::rgb(red = gg_cols[1, 
                                                                                l], green = gg_cols[2, l], blue = gg_cols[3, 
                                                                                                                          l]))
            }
            if (!missing(sz_ax)) {
              graphics::box(lwd = sz_ax)
              graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_x_den)
              graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_y_den)
            }
          }
          else {
            dens <- stats::density(rbind(tmlt))
            if (!is.null(priors) & post_zm == FALSE & is.null(ylim) & 
                is.null(xlim)) {
              ylim <- range(c(range(dens$y), PPO_y_rng))
              xlim <- range(c(range(dens$x), PPO_x_rng))
            }
            else {
              if (!is.null(ylim)) {
                ylim <- YLIM
                xlim <- XLIM
              }
              if (is.null(ylim) & is.null(xlim)) {
                ylim <- NULL
                xlim <- NULL
              }
            }
            graphics::plot(stats::density(rbind(tmlt)), 
                           xlab = xlab_den, ylab = ylab_den, ylim = ylim, 
                           col = COL_DEN, xlim = xlim, lty = lty_den, 
                           lwd = lwd_den, main = "", cex.axis = sz_tick_txt, 
                           cex.lab = sz_ax_txt, xaxt = XAXT_DEN, yaxt = YAXT_DEN)
            graphics::title(main = MAIN_DEN(np[j], main_den, 
                                            j), line = MN_LINE, cex.main = sz_main_txt)
            graphics::axis(side = 1, at = pos_tick_x_den, 
                           cex.axis = sz_tick_txt)
            graphics::axis(side = 2, at = pos_tick_y_den, 
                           cex.axis = sz_tick_txt)
            if (!missing(sz_ax)) {
              graphics::box(lwd = sz_ax)
              graphics::axis(1, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_x_den)
              graphics::axis(2, lwd.ticks = sz_ax, labels = FALSE, 
                             at = pos_tick_y_den)
            }
          }
          if (!is.null(priors)) {
            graphics::lines(dpr, col = COL_PR, lwd = lwd_pr, 
                            lty = lty_pr)
            if (!is.null(sz_txt) & !is.null(COL_TXT)) {
              graphics::legend("topright", legend = ovrlap, 
                               bty = "n", pch = NA, text.col = COL_TXT, 
                               cex = sz_txt)
            }
          }
          if (Rhat == TRUE & n.eff == TRUE) {
            diag_txt <- list(paste0("Rhat: ", rhat[j]), 
                             paste0("n.eff: ", neff[j]))
          }
          if (Rhat == TRUE & n.eff == FALSE) {
            diag_txt <- paste0("Rhat: ", rhat[j])
          }
          if (Rhat == FALSE & n.eff == TRUE) {
            diag_txt <- paste0("n.eff: ", neff[j])
          }
          if (Rhat == TRUE | n.eff == TRUE) {
            if (!is.null(sz_txt) & !is.null(COL_TXT)) {
              graphics::legend("topleft", x.intersp = -0.5, 
                               legend = diag_txt, bty = "n", pch = NA, 
                               text.col = COL_TXT, cex = sz_txt)
            }
          }
          if (!is.null(gvals)) {
            if (length(gvals) == 1) {
              gv <- gvals
            }
            else {
              gv <- gvals[j]
            }
            graphics::abline(v = gv, lty = 2, lwd = 3, 
                             col = ref_col)
          }
        }
      }
      if (type != "both" & type != "density" & type != "trace") {
        stop("Invalid argument for \"type\". Valid inputs are \"both\", \"trace\", and \"density\".")
      }
      if (pdf == TRUE) {
        invisible(grDevices::dev.off())
        if (open_pdf == TRUE) {
          system(paste0("open ", file_out))
        }
      }
      else {
        graphics::par(.pardefault)
      }
    }
    if (plot == FALSE) {
      for (j in 1:length(np)) {
        tmlt <- do.call("cbind", object2[it, np[j]])
        if (!is.null(priors)) {
          if (NCOL(priors) == 1) {
            wp <- priors
          }
          else {
            wp <- priors[, j]
          }
          lwp <- length(wp)
          if (lwp > length(it) * n_chains) {
            pit <- (lwp - (length(it) * n_chains) + 1):lwp
            wp2 <- wp[pit]
          }
          if (lwp < length(it) * n_chains) {
            samps <- sample(wp, size = ((length(it) * n_chains) - 
                                          lwp), replace = TRUE)
            wp2 <- c(wp, samps)
          }
          if (lwp == length(it) * n_chains) {
            wp2 <- wp
          }
          dpr <- stats::density(wp2)
          PPO_x_rng <- range(dpr$x)
          PPO_y_rng <- range(dpr$y)
          tmlt_1c <- matrix(tmlt, ncol = 1)
          pp <- list(wp2, tmlt_1c)
          ovr_v <- round((overlapping::overlap(pp)$OV[[1]]) * 
                           100, digits = 1)
          ovrlap <- paste0(ovr_v, "% overlap")
          if (PPO_out == TRUE) {
            PPO_df$param[j] <- np[j]
            PPO_df$percent_PPO[j] <- ovr_v
          }
        }
      }
    }
    if (PPO_out == TRUE) {
      if (is.null(priors)) {
        warning("NAs produced for PPO dataframe as priors not specified.")
      }
      return(PPO_df)
    }
  }

median_data <- function(rj.obj, model.id){
  
  purrr::map_dbl(.x = model.id, 
                 .f = ~{
                   
                   input.data <- tibble::tibble(species = rj.obj$dat$species$trials,
                                                ID = model.id[rj.obj$dat$species$trials],
                                                y = rj.obj$dat$obs$y_ij,
                                                censored = rj.obj$dat$obs$censored,
                                                rc = rj.obj$dat$obs$Rc,
                                                lc = rj.obj$dat$obs$Lc)
                   
                   sp.data <- input.data %>% 
                     dplyr::filter(ID == .x)
                   
                   if(all(is.na(sp.data$y))){

                       sp.data %>% 
                         dplyr::rowwise() %>% 
                         dplyr::mutate(y = dplyr::case_when(
                           censored == 1 ~ runif(n = 1, min = rc, max = rj.obj$dat$param$dose.range[2]),
                           censored == -1 ~ runif(n = 1, min = rj.obj$dat$param$dose.range[1], max = lc),
                           TRUE ~ runif(n = 1, min = rj.obj$dat$param$dose.range[1], max = rj.obj$dat$param$dose.range[2]))) %>% 
                       dplyr::ungroup() %>% 
                       dplyr::pull(y) %>% 
                       median(., na.rm = TRUE)
                   
                   } else {
                   
                     median(sp.data$y, na.rm = TRUE) }})
}

relabel <- function(vec){
  vec.list <- partitions::vec_to_eq(vec)[]
  vec.df <- tibble::tibble(num = seq_len(length(vec.list)))
  vec.df$len <- sapply(vec.list, length)
  vec.df$first <- purrr::map_dbl(.x = vec.list, .f = ~dplyr::first(sort(.x)))
  new.order <- vec.df %>% dplyr::arrange(-len, first) %>% dplyr::pull(num)
  vec.list <- vec.list[new.order]
  res <- format_group(vec.list)
  return(res$group)
}

vec_to_model <- function(input.vector, sp.names){
  partitions::vec_to_eq(vec = input.vector) %>% 
    list() %>% 
    purrr::map_depth(.x = ., .depth = 2, .f = ~sp.names[.x]) %>%
    purrr::map(.x = ., .f = ~ lapply(X = .x, FUN = function(x) paste(x, collapse = ","))) %>%
    purrr::map(.x = ., .f = ~ paste0("(", .x, ")")) %>% 
    purrr::map(.x = ., .f = ~paste0(.x, collapse = "+")) %>% do.call(c, .)
}

nb_groups <- function(vec){
  length(unique(vec))
}

l_groups <- function(vec){
  purrr::map_dbl(.x = unique(sort(vec)), .f = ~length(vec[vec == .x]))}

N_groups <- function(rj.obj, vec){
  
  purrr::map_dbl(.x = unique(sort(vec)),
                 .f = ~sum(rj.obj$dat$species$nper[which(vec == .x)]))
  
}

# Inverse of intersect
outersect <- function(x, y) {
  sort(c(setdiff(x, y), setdiff(y, x)))
}

# Parallel processing

start_cluster <- function(n.cores){
  cl <<- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)}

stop_cluster <- function(worker = cl){
  parallel::stopCluster(cl = worker) # Stop the cluster
  rm(cl, envir = .GlobalEnv)
  gc()
}


split_count <- function(speciesFrom){
  
  comb.list <- lapply(X = seq_len(length(speciesFrom) - 1),
                      FUN = function(N){
                        combn(x = speciesFrom, m = N) })
                        # gRbase::combn_prim(x = speciesFrom, m = N) })
  
  res <- purrr::map(.x = comb.list, 
                    .f = ~{
                      lapply(X = seq_len(ncol(.x)), 
                             FUN = function(x){
                               sort(c(paste0(.x[, x], collapse = ","),
                                      paste0(speciesFrom[which(!speciesFrom %in% .x[,x])], collapse = ",") ))})}) %>% purrr::flatten() 
  
  res.char <- purrr::map(.x = res, .f = ~paste0(.x, collapse = "+")) %>% 
    unlist()
  
  dd <- duplicated(res.char)
  
  res.char <- res.char[dd]
  res.num <- purrr::map(.x = res, .f = ~stringr::str_split(.x, ",")) %>% 
    purrr::map_depth(.x = ., .depth = 2, .f = ~as.numeric(.x))
  res.num <- res.num[dd]
  
  return(res.num)
}

merge_count <- function(speciesFrom){
  
  # Use simplify = FALSE to return a list, even if only one single pairing is possible
  comb.list <- purrr::map(.x = speciesFrom, .f = ~paste0(.x, collapse = ",")) %>%
    unlist() %>% 
    combn(x = ., m = 2, simplify = FALSE)
    # gRbase::combn_prim(x = ., m = 2, simplify = FALSE)
  return(comb.list)
}

acceptance_rate <- function(AR.obj, mp, cov.n = NULL){
  
  ind <- which(grepl(pattern = "move.", x = names(AR.obj$accept)))
  if(!is.null(cov.n)) ind <- c(ind, which(names(AR.obj$accept) %in% names(cov.n)))
  
  AR.obj$accept[-ind] <- 
    purrr::map(.x = AR.obj$accept[-ind],
               .f = ~round( .x/ mp$n.iter, 3))
  
  if(!is.null(cov.n)){
    AR.obj$accept[names(cov.n)] <- purrr::map(.x = names(cov.n),
               .f = ~unname(round(AR.obj$accept[[.x]]/ cov.n[.x], 3)))
    }
  
  if("move.0" %in% names(mp$move$m)){
    if(mp$move$m["move.0"] > 0){
      AR.obj$accept[["move.0"]] <- round(AR.obj$accept[["move.0"]] / mp$move$m["move.0"], 3)}}
  
  if("move.1" %in% names(mp$move$m)){
    if(mp$move$m["move.1"] > 0){
    AR.obj$accept[["move.1"]] <- round(AR.obj$accept[["move.1"]] / mp$move$m["move.1"], 3)}}
  
  if("move.2" %in% names(mp$move$m)){
    if(mp$move$m["move.2"] > 0){
    AR.obj$accept[["move.2"]] <- round(AR.obj$accept[["move.2"]] / mp$move$m["move.2"], 3)}}
  
  if("move.3" %in% names(mp$move$m)){
    if(mp$move$m["move.3"] > 0){
    AR.obj$accept[["move.3"]] <- round(AR.obj$accept[["move.3"]] / mp$move$m["move.3"], 3)}}
  
  if ("move.covariates" %in% names(AR.obj$accept)) {
    AR.obj$accept[["move.covariates"]] <-
      round(AR.obj$accept[["move.covariates"]] / mp$n.iter, 3)
  }
  
  return(AR.obj)
}

format_group <- function(input.list){
  
  # Sort by vector length
  len <- sapply(input.list, length)
  if(!identical(len, sort(len, decreasing = TRUE))) input.list <- input.list[order(-len)]
  
  # Sort by species
  input.list <- lapply(X = input.list, FUN = sort)
  len <- sapply(input.list, length)
  
  input.list.over <- input.list[which(len > 1)]
  input.list.one <- input.list[which(len == 1)]
  
  if(length(input.list.one) > 0) input.list.one <- input.list.one[order(unlist(input.list.one))]
  
  output.list <- append(input.list.over, input.list.one)
  output.list <- purrr::set_names(output.list, nm = seq_len(length(input.list)))
  
  list.order <- sapply(X = seq_len(max(unlist(output.list))), FUN = function(x) which(unlist(output.list)==x))
  
  # list.order <- sapply(X = seq_len(rj$dat$species$n), FUN = function(x) which(unlist(output.list)==x))
  
  groupings <- list()
  for(u in 1:length(output.list)) groupings[[u]] <- rep(u, length = length(output.list[[u]]))
  
  groupings <- unlist(groupings)[list.order]
  # groupings <- unlist(groupings)[unlist(output.list)]
  
  return(list(species = output.list, group = groupings))
}


# Extract every nth element from a vector
nth_element <- function(vector, starting.position, n) {
  vector[seq(starting.position, length(vector), n)] 
}

glance <- function(dat, which.chain = 1, f = "head", start = 1, end = 6, reset = TRUE) {
  options(max.print = 999999)
  dat[[which.chain]]$mcmc$move$m <- NULL
  dat <- dat[[which.chain]]
  if(f == "update"){
    dat.list <- dat[!names(dat) %in% c("abbrev", "mcmc", "run_time", "accept", "config", "dat")]
    res <- purrr::map(.x = dat.list, .f = ~ {
      if(is.null(dim(.x))){
        out <- .x[length(.x)] 
      } else {
        if(length(dim(.x)) == 2) out <- .x[nrow(.x), ] 
        if(length(dim(.x)) == 3) out <- .x[nrow(.x), ,]
      }
    out})
    res$mlist <- dat$mlist
    res$iter <- dat$iter
    return(res)
  } else {
    dat <- dat[!names(dat)%in%c("accept", "mlist", "config", "dat")]
    purrr::map(.x = dat, .f = ~ tail(get(f)(.x, end), (end - start) + 1))}
}

outOfbounds <- function(rj.obj, v, p){
  if(p == "alpha") {
    any(v < rj.obj$config$priors[p, 1]) | 
      any(v > rj.obj$config$priors[p, 2]) | 
      any(v < rj.obj$nu[rj.obj$iter["nu"], , 1]) |
      any(v > rj.obj$nu[rj.obj$iter["nu"], , 2])
  } else if(p == "nu1"){
    any(v < rj.obj$config$priors["nu", 1]) | any(v > rj.obj$alpha[rj.obj$iter["alpha"], ])
  } else if(p == "nu2"){
    any(v < rj.obj$alpha[rj.obj$iter["alpha"], ]) | any(v > rj.obj$config$priors["nu", 2])
  } else if(p == "mu"){
    any(v < rj.obj$config$priors["mu", 1]) | any(v > rj.obj$config$priors["mu", 2])
  } else {
    any(v < rj.obj$config$priors[p, 1]) | any(v > rj.obj$config$priors[p, 2])
  }
}

# Convenience function to list properties associated with a categorical covariate
# nL = Number of levels of the factor
# Lnames = Names of factor levels.
# nparam = Number of levels other than the baseline
# and index = Column indices for those levels.

factor_levels <- function(covname, dat){
  
  # Number of factor levels
  nL <- nlevels(dat[, covname])
  
  # Number of levels other than baseline
  if(nL == 0) nparam <- 1 else nparam <- nL - 1
  
  index <- lapply(X = nL, FUN = function(x){
    if(x == 0) res <- 1
    if(x == 2) res <- 2
    if(x > 2) res <- seq(from = x-1, to = x)
    return(res)})
  
  return(list(nL = nL, Lnames = levels(dat[, covname]),
              nparam = nparam, index = index[[1]]))
  
}