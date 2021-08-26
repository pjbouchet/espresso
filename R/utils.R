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

# Likelihood --------------------------------------------------------------

#' Likelihood
#'
#' Calculate log-likelihoods required by the rjMCMC sampler.
#' 
#' @param rj.obj rjMCMC object.
#' @param iter Current iteration.
#' @param model Model ID number.
#' @param values List of proposed values, if available.
#' @param include.cov Boolean vector indicating which contextual covariates are included in the current iteration of the rjMCMC sampler.
#' @param param Paramter name. Ignored if \code{values} is specified.
#' @param RJ Logical. If \code{TRUE}, returns the likelihoods associated with between-model jumps. If \code{FALSE}, returns the likelihoods associated with the update step of the Metropolis sampler.
#' @param lprod If \code{TRUE}, returns the product of individual likelihoods. 

likelihood <- function(rj.obj,
                       iter,
                       model, 
                       values = list(), 
                       include.cov = NULL, 
                       param.name = NULL, 
                       RJ = FALSE, 
                       lprod = TRUE){
  
  # Define model parameters
  # Must be in the order in which they are updated within the MCMC sampler
  model.parameters <- c("y.ij", "t.ij", "sigma", rj.obj$dat$covariates$names, "mu.i", "mu", "phi")
  not.covariates <- c("y.ij", "t.ij", "sigma", "mu.i", "mu", "phi")
  
  if (length(values) > 0) {
    par.name <- names(values)
  } else {
    if (is.null(param.name)) par.name <- "" else par.name <- param.name
  }
  
  # Check inputs are correct
  if (all(!par.name == "" & !par.name %in% model.parameters)) stop("Parameter(s) not found.")
  if (!RJ & length(values) == 0 & is.null(param.name)) stop("Must provide parameter name.")
  
  # Determine correct indices to retrieve to relevant parameter values 
  iteration <- rep(iter - 1, times = length(model.parameters))
  iter.ind <- which(model.parameters %in% as.character(par.name))
  iter.ind <- ifelse(length(iter.ind) == 0, 1, iter.ind)
  if (!RJ) iteration[(1:iter.ind) - 1] <- iter
  if (!RJ | !is.null(include.cov)) iteration[which(model.parameters == "mu")] <- iter
  
  # When updating covariates in MH step
  if (!RJ & is.null(include.cov) & all(!par.name %in% not.covariates)) 
    iteration[which(model.parameters %in% rj.obj$dat$covariates$names)] <- iter
  
  # Function to extract the right parameter values
  extract_pars <- function(m){
    x <- rj.obj[[m]]
    if("matrix" %in% class(x)) x <- x[iteration[which(model.parameters == m)], ] else 
      x <- x[iteration[which(model.parameters == m)]]}
  
  
  # Extract parameter values
  params <- lapply(X = model.parameters, FUN = extract_pars)
  names(params) <- model.parameters
  
  if (any(!par.name == "")) params[names(values)] <- values[]
  
  # include.cov is used to update covariate values
  if (!is.null(include.cov)){
    if (sum(include.cov) == 0) params$include.covariates <- "None" else 
      params$include.covariates <- names(include.cov[include.cov == 1])
  } else { # if include.cov is NULL
    if(RJ) whichrow <- iter - 1 else whichrow <- iter
    include.cov <- rj.obj$include.covariates[whichrow, ]
    params$include.covariates <- names(include.cov[include.cov == 1])
  }
  
  # Log-likelihoods
  if(!lprod) loglikelihood <- list() else loglikelihood <- numeric(length = 2 + nb_groups(model)) 
  
  if(all(par.name %in% "t.ij") | RJ){
    LL.1 <- dnorm(x = params$y.ij, mean = params$t.ij, sd = rj.obj$dat$obs$sd, log = TRUE)
    if(!lprod) loglikelihood <- append(loglikelihood, list(LL.1)) else loglikelihood[1] <- sum(LL.1)}
  
  if(all(par.name %in% c("t.ij", "sigma", "mu.i", rj.obj$dat$covariates$names)) | RJ){
    
    LL.2 <- dtnorm(
      x = params$t.ij,
      location = Reduce("+", append(
        list(params$mu.i[rj.obj$dat$whales$id]),
        lapply(
          X = params$include.covariates,
          FUN = function(k) {
            if (k == "None") {
              0
            } else {
              apply(t(t(rj.obj$dat$covariates$dummy[[k]]) * params[[k]]), 1, sum)
            }
          }
        )
      )),
      scale = params$sigma,
      L = rj.obj$dat$param$bounds["t.ij", 1],
      U = rj.obj$dat$param$bounds["t.ij", 2],
      log = TRUE
    )
    
    if (!lprod) {
      if (all(par.name == "mu.i")) {
        LL.2 <- purrr::map_dbl(.x = seq_len(rj.obj$dat$whales$n), 
                               .f = ~ sum(LL.2[rj.obj$dat$whales$id == .x]))
      }
      
      loglikelihood <- append(loglikelihood, list(LL.2))
    } else {
      loglikelihood[2] <- sum(LL.2)
    }
  }
  
  
  if(all(par.name %in% c("mu.i", "mu", "phi")) | RJ){

    LL.3 <- sapply(X = 1:nb_groups(model), 
                   FUN = function(p){
                     dtnorm(
                       x = params$mu.i[model[rj.obj$dat$species$id] == p],
                       location = unique(params$mu[model == p]),
                       scale = params$phi,
                       L = rj.obj$dat$param$bounds["mu.i", 1],
                       U = rj.obj$dat$param$bounds["mu.i", 2],
                       log = TRUE)
                   })
    
    if (!lprod) {
      
      if (all(par.name == "mu.i")) {
        mu.i.vec <- sort(as.numeric(gsub(
          pattern = "mu.", replacement = "",
          x = names(unlist(c(setNames(LL.3, rownames(LL.3)))))
        )),
        index.return = TRUE
        )
        loglikelihood <- append(loglikelihood, list(unlist(LL.3)[mu.i.vec$ix]))
      } else {
        loglikelihood <- append(loglikelihood, list(LL.3))
      }
      
    } else {
      
      if(is.list(LL.3)) LL.3 <- sapply(LL.3, sum) else LL.3 <- colSums(LL.3)
      
      if (nb_groups(model) == 1) {
        loglikelihood[2 + seq_len(nb_groups(model))] <- sum(LL.3)
      } else {
        loglikelihood[2 + seq_len(nb_groups(model))] <- LL.3
      }
    }
  }
  
  if(!lprod) return(unname(Reduce("+", loglikelihood))) else return(sum(loglikelihood, na.rm = TRUE))
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
  
  loglik.unif <- dunif(x = unique(param), 
                       min = rj.obj$dat$param$bounds[param.name, 1],
                       max = rj.obj$dat$param$bounds[param.name, 2],
                       log = TRUE)
  
  if (any(abs(loglik.unif) == Inf)) loglik.unif <- -100000
  return(sum(loglik.unif))
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
                       mean = rj.obj$config$prior[[param.name]][1], 
                       sd = rj.obj$config$prior[[param.name]][2], 
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
propose_jump <- function(rj.obj, move.type) {
  
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
      J <- (n[1] + n[2]) / (n[1] * n[2])
      
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
      J <- (n[1] * n[2]) / (n[1] + n[2])
      
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

#' Proposal for between-model jumps
#'
#' Generate proposed values for relevant model parameters
#' @param jump Proposed model jump, as defined by \code{\link{propose_jump}}.

proposal_rj <- function(rj.obj, jump, iter) {

  # Current mu
  rj.means <- rj.obj$mu[iter - 1, ]
  
  # Bootstrap moves
  if(jump$type == 1){
    
    u <- rnorm(n = length(unique(jump$model$id)), mean = 0, sd = rj.obj$config$prop$dd)
    rj.means <- median_data(rj.obj = rj.obj, model.id = jump$model$id) + u[jump$model$id]
    g_u <- NULL
    
    # Cluster and Random moves
  } else if (jump$type %in% c(2, 3)) {
    
    rj.means <- sapply(X = 1:nb_groups(jump$model$id), 
                       FUN = function(m) runif(n = 1, min = rj.obj$dat$param$bounds["mu", 1],
                                               max = rj.obj$dat$param$bounds["mu", 2]))[jump$model$id]
    g_u <- u <- NULL
    
    # Split / merge move
  } else if (jump$type == 0){
    
    # Split move
    if (jump$move == 1) {
      
      g_u <- numeric(2)
      
      g_u[1] <- unname(min(c(jump$n[1] * (rj.obj$dat$param$bounds["mu", 2] - rj.means[jump$species$from[1]]),
                             jump$n[2] * (rj.means[jump$species$from[2]] - rj.obj$dat$param$bounds["mu", 1]))))
      
      g_u[2] <- unname(max(c(-jump$n[2] * (rj.obj$dat$param$bounds["mu", 2] - rj.means[jump$species$from[2]]),
                             -jump$n[1] * (rj.means[jump$species$from[1]] - rj.obj$dat$param$bounds["mu", 1]))))
      
      g_u <- sort(g_u)
      
      u <- runif(n = 1, min = g_u[1], max = g_u[2])
      
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
      
      # Or: u <- (wm - m[2]) * n[2]
      # Or: u = (m[2] - m[1]) * ((n[1]*n[2])/(n[1]+n[2]))
      
    } # End move == 2
  } # End random
  
  return(list(mu = rj.means, u = u, g_u = g_u, n = jump$n))
}

#' Proposal for within-model jumps
#'
#' Generate proposed values for relevant model parameters
#'
#' @param param.name Parameter name.
#' @param seed Random seed. Only used for testing purposes.

proposal_mh <- function(rj.obj, param.name, iter, seed = NULL) {
  
  # Set the random seed
  if(!is.null(seed)) set.seed(seed) 
  
  # Species groups
  ng <- unique(rj.obj$mlist[[rj.obj$current.model]])
  
  # Number of proposal values
  if(param.name == "t.ij") N <- rj.obj$dat$trials$n else 
    if(param.name == "mu.i") N <- rj.obj$dat$whales$n else 
      if(param.name == "mu") N <- length(ng) else 
        if(param.name %in% c("phi", "sigma")) N <- 1 else
          if(param.name %in% rj.obj$dat$covariates$names) N <- rj.obj$dat$covariates$fL[[param.name]]$nparam
  
  # Mean(s) for proposal(s)
  if (param.name %in% c("t.ij", "mu.i")) m <- rj.obj[[param.name]][iter - 1, ]
  if (param.name %in% c("mu")) m <- rj.obj[[param.name]][iter, ]
  if (param.name %in% c("sigma", "phi")) m <- rj.obj[[param.name]][iter - 1]
  if (param.name %in% rj.obj$dat$covariates$names) m <- rj.obj[[param.name]][iter, rj.obj$dat$covariates$fL[[param.name]]$index]
  
  lower.limit <- rep(rj.obj$dat$param$bounds[param.name, 1], N)
  upper.limit <- rep(rj.obj$dat$param$bounds[param.name, 2], N)
  
  # Generate proposal(s)
  if(param.name == "t.ij"){
    
    lower.limit[rj.obj$dat$obs$censored == 1] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]
    upper.limit[rj.obj$dat$obs$censored == -1] <- rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]
   
  }
  
  if(param.name %in% c("t.ij", "mu.i")){
    
    pval <- rtnorm(n = N, 
                   location = m, 
                   scale = rj.obj$config$prop$mh[[param.name]], 
                   L = lower.limit,
                   U = upper.limit)
    
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
  if(any(is.infinite(pval))) stop("Infinite values generated")
  return(pval)
}

# Proposal densities ----------------------------------------------------------------

#' Proposal densities for between-model jumps
#'
#' @param param Parameter values.
#' @param jump Proposed between-model jump.

propdens_rj <- function(rj.obj, param, jump, iter) {
  
  # Adapted to work with Uniform proposal distribution as per Huelsenbeck et al.
  # aux is the value of the auxiliary variable 'u' drawn by proposal_rj
  # during a proposed split/merge move. u does not exist when the independent
  # sampler is used.
  
  # Cluster and Random moves
  if(jump$type %in% c(2, 3)){
    
    loglik.unif <- list()
    
    loglik.unif[[1]] <- dunif(x = unique(param$mu), 
                              min = rj.obj$dat$param$bounds["mu", 1],
                              max = rj.obj$dat$param$bounds["mu", 2],
                              log = TRUE)
    
    loglik.unif[[2]] <- dunif(x = unique(rj.obj$mu[iter - 1, ]), 
                              min = rj.obj$dat$param$bounds["mu", 1],
                              max = rj.obj$dat$param$bounds["mu", 2],
                              log = TRUE)
    
    loglik.unif <- purrr::map_dbl(.x = loglik.unif, 
                                  .f = ~ifelse(any(abs(.x) == Inf), -100000, sum(.x)))
    
    return(unlist(loglik.unif))
    
    # Bootstrap moves
  } else if (jump$type == 1) {
    
    loglik.norm.forward <- dnorm(x = param$u,
                                 # x = unique(param$u[jump$model$id]), 
                                 mean = 0, 
                                 sd = rj.obj$config$prop$dd, 
                                 log = TRUE)
    
    if (any(abs(loglik.norm.forward) == Inf)) loglik.norm.forward <- -100000
    loglik.norm.forward <- sum(loglik.norm.forward)
    
    loglik.norm.backward <- dnorm(x = median_data(rj.obj = rj.obj,
                                                model.id = rj.obj$mlist[[rj.obj$current.model]]) - 
                                    rj.obj$mu[iter - 1, ],
                                  mean = 0, 
                                  sd = rj.obj$config$prop$dd,
                                  log = TRUE)
    
    if (any(abs(loglik.norm.backward) == Inf)) loglik.norm.backward <- -100000
    loglik.norm.backward <- sum(loglik.norm.backward)
    
    return(c(loglik.norm.backward, loglik.norm.forward))
    
    # Split / merge moves
  } else if (jump$type == 0) {
    
    loglik.unif <- dunif(x = unique(param$u), min = param$g_u[1], max = param$g_u[2], log = TRUE)
    if (any(abs(loglik.unif) == Inf)) loglik.unif <- -100000
    
    # Return two values for the ratio of proposal densities
    # Depending on the type of move, either the numerator or
    # denominator is 0 (as we work in logs)
    if (jump$move == 1) return(c(0, loglik.unif)) else return(c(loglik.unif, 0))
  }
}

#' Proposal densities for between-model jumps
#'
#' @param param.name Parameter name
#' @param dest Proposed value(s)
#' @param orig Current value(s) 

propdens_mh <- function(rj.obj, param.name, dest, orig) {

  lower.limit <- rep(rj.obj$dat$param$bounds[param.name, 1], length(dest))
  upper.limit <- rep(rj.obj$dat$param$bounds[param.name, 2], length(dest))
  
  if(param.name == "t.ij") {
    
    lower.limit[rj.obj$dat$obs$censored == 1] <- rj.obj$dat$obs$Rc[rj.obj$dat$obs$censored == 1]
    upper.limit[rj.obj$dat$obs$censored == -1] <- rj.obj$dat$obs$Lc[rj.obj$dat$obs$censored == -1]

    }
  
  if(param.name %in% c("t.ij", "mu.i")) {
    
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
  if(param.name %in% c("t.ij", "mu.i")) return(loglik) else return(sum(loglik))
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
    
    # model.ID.posterior <- as.numeric(names(res$model$m_prob))[order(res$model$m_prob, decreasing = TRUE)]
    # bestmod <- as.numeric(names(which.max(res$model$m_prob))); names(bestmod) <- NULL
    
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
    
    # # Assign names
    # if (select) {
    #   names(res$model$m_prob) <- 
    #     sapply(X = names(res$model$m_prob), FUN = function(nn) mlist$model[as.numeric(nn)])
    # } else {
    #   names(res$model$m_prob) <- mlist$model
    # }
    

    # res$model$m_prob <- sort(res$model$m_prob, decreasing = TRUE)
    # res$model$bestmod <- names(which.max(res$model$m_prob))
    # res$model$m_prob <- tibble::enframe(res$model$m_prob) %>% 
    #   dplyr::rename(model = name, p = value) %>% 
    #   dplyr::mutate(p = as.numeric(p)) %>% 
    #   dplyr::slice(1:n.top)
    
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
                   legend.title = element_blank(),
                   legend.text = element_text(size = 12),
                   panel.background = element_rect(fill = "white")) + 
    
    # Annotation for posterior probabilities
    ggplot2::geom_text(data = p1, aes(x = x , y = y, label = label), size = p.size, hjust = 0) +
    ggplot2::geom_text(data = p2, aes(x = x , y = y, label = label), size = p.size) +
    
    ggplot2::guides(fill = guide_legend(nrow = 1)) +
    
    {if(!combine) ggplot2::theme(legend.position = "none")} +
    {if(!combine) ggtitle(paste0("Rank: ", no)) }
  
  return(out.plot)
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
          return(paste0("Density - ", x))
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
          return(paste0("Trace - ", x))
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
                       dplyr::mutate(y = ifelse(censored == 1,
                                     runif(n = 1, min = rc, max = rj.obj$dat$param$bounds["mu", 2]),
                                     runif(n = 1, min = rj.obj$dat$param$bounds["mu", 1], max = lc))) %>% 
                       dplyr::ungroup() %>% 
                       dplyr::pull(y) %>% 
                       mean(., na.rm = TRUE)
                   
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
                        utils::combn(x = speciesFrom, m = N) })
  
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
    unlist() %>% utils::combn(x = ., m = 2, simplify = FALSE)
  return(comb.list)
}

acceptance_rate <- function(AR.obj, rj.obj, mp){
  
  AR.obj$accept[1:(5 + rj.obj[[1]]$dat$covariates$n)] <- 
    purrr::map(.x = AR.obj$accept[1:(5 + rj.obj[[1]]$dat$covariates$n)],
               .f = ~round( .x/ mp$n.iter, 3))
  
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
      if(is.null(dim(.x))) .x[length(.x)] else .x[nrow(.x), ]
    })
    res$mlist <- dat$mlist
    return(res)
  } else {
    dat <- dat[!names(dat)%in%c("accept", "mlist", "config", "dat")]
    purrr::map(.x = dat, .f = ~ tail(get(f)(.x, end), (end - start) + 1))}
}

outOfbounds <- function(rj.obj, v, p){
  (any(v < rj.obj$dat$param$bounds[p, 1] | v > rj.obj$dat$param$bound[p, 2]))
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