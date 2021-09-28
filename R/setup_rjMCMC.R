#' Set up the rjMCMC sampler
#'
#' Set up a list object to hold the data from the MCMC sampler, and generate starting values for all model parameters. This function is called internally by \code{\link{run_rjMCMC}}.
#'
#' In addition to split/merge moves, three other types of MCMC samplers are implemented in \code{espresso} to facilitate convergence and avoid getting stuck in local maxima: two data-driven samplers (type I and type II), in which proposals are informed by the “cues” present in the original data, and one independent sampler, in which proposals are drawn at random from a Uniform distribution bounded by \code{range.dB} (see \code{\link{read_data}} or \code{\link{simulate_data}}), with no dependence on current values.
#'
#' @param rj.input Input dataset. Must be an object of class \code{rjconfig}.
#' @param n.burn Number of MCMC iterations to treat as burn-in.
#' @param n.iter Number of posterior samples.
#' @param n.chains Number of MCMC chains. Defaults to 3.
#' @param do.update Logical. Whether to update an existing sampler or set up a new one.
#' @param start.values Starting values to use when updating an existing sampler and \code{do.update = TRUE}.
#' @inheritParams run_rjMCMC
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}}
#' @keywords brs dose-response rjmcmc 

setup_rjMCMC <- function(rj.input, 
                         n.burn,
                         n.iter,
                         n.chains = 3,
                         p.split,
                         p.merge,
                         moves,
                         move.ratio,
                         m,
                         do.update = FALSE,
                         start.values = NULL) {
  
  #' -----------------------------------------------
  # Perform function checks
  #' -----------------------------------------------
  
  if(!do.update & !"rjconfig" %in% class(rj.input)) 
    stop("Input data must be initialised using <rj_initialise>.")
  if(is.null(start.values) & do.update) stop("Missing starting values.")
  if(!is.null(move.ratio) & !"integer" %in% class(move.ratio) & !"list" %in% class(move.ratio)) stop("<move.ratio> must be either an integer vector or a list.")
  if(is.list(move.ratio)) move.ratio <- unlist(move.ratio)
  if(do.update) tot.iter <- n.iter else tot.iter <- n.iter + n.burn # Total chain length
  
  #' -----------------------------------------------
  # Define move types
  #' -----------------------------------------------
  
  # Type 0: Split/merge
  # Type 1: Data-driven Type I
  # Type 2: Data-driven Type II
  # Type 3: Random
  
  modelmoves <- numeric(length = tot.iter)
  
  if(length(moves) == 1){
    modelmoves <- rep(moves, tot.iter)
  } else if (length(moves) > 1){
    if(m < tot.iter){
      move.iters <- seq(from = m, to = tot.iter, by = m)
    } else {
      move.iters <- rep(0, tot.iter)
    }
    move.ratio <- unlist(move.ratio)
    modelmoves[move.iters] <- rep(unlist(sapply(X = which(move.ratio > 0), 
                                                FUN = function(a) rep(a, each = move.ratio[a]))),
                                  length.out = length(move.iters))
  }
  
  modelmoves.tab <- table(modelmoves[(n.burn+1):tot.iter]) %>% 
    tibble::enframe() %>% 
    dplyr::rename(move = name, n = value) %>% 
    dplyr::mutate(move = dplyr::case_when(move == "0" ~ "split-merge",
                                          move == "1" ~ "data-driven (Type I)",
                                          move == "2" ~ "data-driven (Type II)",
                                          move == "3" ~ "random")) 
  
  rj <- list(mcmc = list(n.chains = n.chains,
                         n.burn = n.burn,
                         n.iter = n.iter,
                         tot.iter = tot.iter,
                         iter.rge = paste0(n.burn + 1, ":", tot.iter),
                         move = list(prob = setNames(c(p.split, p.merge), c("split", "merge")),
                                     m = modelmoves,
                                     tab = modelmoves.tab)))
  
  #' -----------------------------------------------
  # Create data matrices and vectors
  #' -----------------------------------------------
  
  rj <- append(rj, list(t.ij = matrix(data = 0, nrow = tot.iter, 
                                      ncol = rj.input$trials$n)))
  
  if(!rj.input$config$biphasic | rj.input$config$function.select){
    
    rj <- append(rj, list(mu.i = matrix(data = 0, nrow = tot.iter, ncol = rj.input$whales$n),
                          mu = matrix(data = 0, nrow = tot.iter, ncol = rj.input$species$n),
                          phi = numeric(tot.iter),
                          sigma = numeric(tot.iter)))
  }
  
  if(rj.input$config$biphasic){
    
    rj <- append(rj, list(mu.ij = array(data = 0, dim = c(tot.iter, rj.input$trials$n, 2)),
                          psi.i = matrix(data = 0, nrow = tot.iter, ncol = rj.input$whales$n),
                          k.ij = matrix(data = 0, nrow = tot.iter, ncol = rj.input$trials$n),
                          pi.ij = matrix(data = 0, nrow = tot.iter, ncol = rj.input$trials$n),
                          nu = array(data = 0, dim = c(tot.iter, rj.input$species$n, 2),
                                     dimnames = list(NULL, 
                                                     paste0(".", seq_len(rj.input$species$n)),
                                                     # paste0(" (", rj.input$species$names, ")"),
                                                     c("nu.lower", "nu.upper"))),
                          tau = matrix(data = 0, nrow = tot.iter, ncol = 2),
                          alpha = matrix(data = 0, nrow = tot.iter, ncol = rj.input$species$n),
                          psi = numeric(tot.iter),
                          omega = numeric(tot.iter)))

  }
 
  # Covariates, if present
  
  if (rj.input$covariates$n > 0) {
    for (j in rj.input$covariates$names) {
      rj[[j]] <- matrix(
        data = 0, nrow = tot.iter,
        ncol = ifelse(rj.input$covariates$fL[[j]]$nL == 0, 1, rj.input$covariates$fL[[j]]$nL)
      )
    }
    
    rj$include.covariates <- matrix(
      data = ifelse(!rj.input$config$covariate.select, 1, 0),
      nrow = tot.iter, ncol = rj.input$covariates$n
    )
    
  } else {
    rj$include.covariates <- matrix(data = 0, nrow = tot.iter, ncol = 1)
    rj$beta <- matrix(data = 0, nrow = tot.iter, ncol = 1)
  }
  
  rj$model <- rj$phase <- numeric(tot.iter)
  rj$mlist <- rj.input$config$mlist

  if(!rj.input$config$function.select){
    if(rj.input$config$biphasic) rj$phase <- rep(2, tot.iter) else rj$phase <- rep(1, tot.iter)
  }
  
  #' -----------------------------------------------
  # Assign column names
  #' -----------------------------------------------

  colnames(rj$t.ij) <- paste0("t", rep(seq_len(rj.input$whales$n), rj.input$trials$nper), ".", 
                              unlist(sapply(rj.input$trials$nper, seq_len)))
  
  if(!rj.input$config$biphasic | rj.input$config$function.select){
  colnames(rj$mu.i) <- paste0("mu.", 1:rj.input$whales$n)
  colnames(rj$mu) <- paste0("mu.", 1:rj.input$species$n)
  }
  
  if(rj.input$config$biphasic){
    # Note: cannot assign column names to mu.ij and nu as these are stored in arrays
    colnames(rj$psi.i) <- paste0("psi.i.", seq_len(rj.input$whales$n))
    colnames(rj$k.ij) <- paste0("k.ij.", seq_len(rj.input$trials$n))
    colnames(rj$pi.ij) <- paste0("pi.ij.", seq_len(rj.input$trials$n))
    colnames(rj$tau) <- c("tau.lower", "tau.upper")
    colnames(rj$alpha) <- paste0("alpha.", seq_len(rj.input$species$n))
    # colnames(rj$alpha) <- paste0("alpha (", rj.input$species$names, ")")
  }
  
  if(is.null(rj.input$species$abbrev)) {
    rj$abbrev <- list(rj.input$species$abbrev) 
    } else { rj$abbrev <- rj.input$species$abbrev }
  
  if (rj.input$covariates$n > 0) {
    for (j in rj.input$covariates$names) {
      if(is.null(rj.input$covariates$fL[[j]]$Lnames)){
        colnames(rj[[j]]) <- j
      } else{
        colnames(rj[[j]]) <- paste0(j, "_", rj.input$covariates$fL[[j]]$Lnames)
      }
    }
    colnames(rj$include.covariates) <- rj.input$covariates$names
  } else {
    colnames(rj$include.covariates) <- "beta"
  }
  
  # Acceptance rates
  rj$accept <- list()
  
  if(!rj.input$config$biphasic | rj.input$config$function.select){
    rj$accept <- append(rj$accept, 
                        list(t.ij = numeric(rj.input$trials$n),
                             mu.i = numeric(rj.input$whales$n), 
                             mu = 0, 
                             phi = 0, 
                             sigma = 0))
  }
  
  if(rj.input$config$biphasic){
    rj$accept <- append(rj$accept, list(mu.ij.1 = numeric(rj.input$trials$n), 
                                        mu.ij.2 = numeric(rj.input$trials$n), 
                                        psi.i = numeric(rj.input$whale$n), 
                                        k.ij = numeric(rj.input$trials$n), 
                                        nu = c(0, 0),
                                        tau = c(0, 0),
                                        alpha = 0, 
                                        omega = 0))}
                   
  if(rj.input$config$covariate.select){
    rj$accept <- append(rj$accept, list(move.covariates = numeric(1)))}
  
  if(rj.input$covariates$n > 0){
    cl <- as.list(numeric(rj.input$covariates$n))
    names(cl) <- rj.input$covariates$names
    rj$accept <- append(rj$accept, cl)
  }
  
  if(rj.input$config$model.select){
    ml <- as.list(rep(0, length(rj.input$config$move$moves)))
    names(ml) <- paste0("move.", rj.input$config$move$moves) 
    if(rj.input$config$function.select){
      ml <- append(ml, ml)
      names(ml) <- paste0(c("mono.", "bi."), names(ml))
    } else {
      if(rj.input$config$biphasic){
        names(ml) <- paste0("bi.", names(ml))
      } else {
        names(ml) <- paste0("mono.", names(ml)) 
      }
    }
    rj$accept <- append(rj$accept, ml)
  }
  
  if(rj.input$config$function.select){
    rj$accept <- append(rj$accept, list(to.monophasic = 0, to.biphasic = 0))
  }

  #' -----------------------------------------------
  # Generate starting values
  #' -----------------------------------------------
  
  if(do.update){
    
    for(it in names(start.values)){
      if(it == "mlist") rj[[it]] <- start.values[[it]]
      if(is.null(dim(rj[[it]]))) rj[[it]][1] <- start.values[[it]][1] 
      if(length(dim(rj[[it]])) == 2) rj[[it]][1, ] <- start.values[[it]]
      if(length(dim(rj[[it]])) == 3) rj[[it]][1, ,] <- start.values[[it]]
    }
    
  } else {
    
    # Model ID
    if (rj.input$config$model.select) {
      rj$current.model <- rj$model[1] <-
        sample(x = rj.input$config$clust[[1]]$model, size = 1, prob = rj.input$config$clust[[1]]$p_scale)
    } else {
      # Because the model list only has length 1
      rj$current.model <- rj$model[1:tot.iter] <- names(rj.input$config$mlist) 
    }
    
    # Functional form
    if(rj.input$config$function.select) rj$phase[1] <- sample(x = 1:2, size = 1)
    
    if(rj.input$covariates$n > 0) {
      
      if(rj.input$config$covariate.select) rj$include.covariates[1, ] <- 
          sample(x = 0:1, size = rj.input$covariates$n, replace = TRUE) 
      
      for (j in rj.input$covariates$names){
        v <- rnorm(n = ifelse(rj.input$covariates$fL[[j]]$nL == 0, 1, 
                              rj.input$covariates$fL[[j]]$nparam), mean = 0, sd = 3)
        if(rj.input$covariates$fL[[j]]$nL == 0){ rj[[j]][1, ] <- v  } else { rj[[j]][1, ] <- c(0, v)}}
      
    }
    
    # Model parameters 
    # -- MONOPHASIC ----
    
    if(!rj.input$config$biphasic | rj.input$config$function.select){
      
    rj$sigma[1] <- rj.input$config$var[1] + 
      rtnorm(n = 1, location = 0, scale = 3, 
             L = rj.input$config$priors["sigma", 1] - rj.input$config$var[1],
             U = rj.input$config$priors["sigma", 2] - rj.input$config$var[1])
    
    rj$phi[1] <- rj.input$config$var[2] + 
      rtnorm(n = 1, location = 0, scale = 3, 
             L = rj.input$config$priors["phi", 1] - rj.input$config$var[2],
             U = rj.input$config$priors["phi", 2] - rj.input$config$var[2])
    
    input.data <- tibble::tibble(species = rj.input$species$trials,
                   y = rj.input$obs$y_ij,
                   censored = rj.input$obs$censored,
                   rc = rj.input$obs$Rc,
                   lc = rj.input$obs$Lc)
    
    mu.start <- purrr::map_dbl(
      .x = seq_len(nb_groups(rj$mlist[[rj$current.model]])),
      .f = ~ {
        
        sp.data <- input.data %>% 
          dplyr::filter(species == .x)
        
        if(all(is.na(sp.data$y))){
          sp.data %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(y = ifelse(censored == 1,
                                     runif(n = 1, min = rc, max = rj.input$config$priors["mu", 2]),
                                     runif(n = 1, min = rj.input$config$priors["mu", 1], max = lc))) %>% 
            dplyr::ungroup() %>% 
            dplyr::pull(y) %>% 
            mean(., na.rm = TRUE)
        } else {
          mean(sp.data$y, na.rm = TRUE)
        }})

    mu.deviates <- rtnorm(n = length(mu.start), 
                   location = 0,
                   scale = 5, 
                   L = rj.input$param$dose.range[1] - mu.start,
                   U = rj.input$param$dose.range[2] - mu.start)
    
    rj$mu[1, ] <- (mu.start + mu.deviates)[rj$mlist[[rj$current.model]]]
    
    # rj$mu[1,] <- mu.start + mu.deviates[sapply(X = mu.start, 
    #                         FUN = function(x) which(unique(mu.start) == x), simplify = TRUE)]
    
    rj$mu.i[1, ] <- rtnorm(n = rj.input$whales$n,
                           location = rj$mu[1, rj.input$species$id], 
                           scale = rj$phi[1],
                           L = rj.input$param$dose.range[1], 
                           U = rj.input$param$dose.range[2])
    
    mu.ij.config <- rj$mu.i[1, rj.input$whales$id]
    
    # The scale needs to be large here to prevent numerical issues (returning Inf)
    rj$t.ij[1, ] <- rtnorm(n = rj.input$trials$n, 
                           location = mu.ij.config, 
                           scale = 30,
                           # rj$sigma[1],
                           L = rj.input$param$dose.range[1], 
                           U = rj.input$param$dose.range[2])

    } 
    
    if(rj.input$config$biphasic | rj.input$config$function.select){
      
    # Model parameters 
    # -- BIPHASIC ----
    
      # inits.bi <- purrr::map(.x = seq_len(nb_groups(rj$mlist[[rj$current.model]])),
      #       .f = ~sort(rj.input$obs$y_ij[rj.input$species$trials %in% which(rj$mlist[[rj$current.model]] == .x)]))
      
      inits.bi <- purrr::map(.x = seq_len(nb_groups(rj$mlist[[rj$current.model]])),
                             .f = ~{
                               censored <- rj.input$obs$censored[rj.input$species$trials %in% 
                                            which(rj$mlist[[rj$current.model]] == .x)]
                               y.obs <- rj.input$obs$y_ij[rj.input$species$trials %in% 
                                           which(rj$mlist[[rj$current.model]] == .x)]
                               Lc.obs <- rj.input$obs$Lc[rj.input$species$trials %in% 
                                           which(rj$mlist[[rj$current.model]] == .x)]
                               Rc.obs <- rj.input$obs$Rc[rj.input$species$trials %in% 
                                           which(rj$mlist[[rj$current.model]] == .x)]
                               y.obs[is.na(y.obs) & censored == 1] <- 
                                 runif(n = sum(is.na(y.obs) & censored == 1),
                                       min = Rc.obs[is.na(y.obs) & censored == 1],
                                       max = rj.input$param$dose.range[2])
                               y.obs[is.na(y.obs) & censored == -1] <- 
                                 runif(n = sum(is.na(y.obs) & censored == -1),
                                       min = rj.input$param$dose.range[1],
                                       max = Lc.obs[is.na(y.obs) & censored == -1])
                               
                               sort(y.obs)})

      rj$alpha[1, ] <- sapply(X = inits.bi, FUN = function(x){
        rtnorm(n = 1, location = median(x, na.rm = TRUE), scale = 1, 
                         L = rj.input$param$dose.range[1],
                         U = rj.input$param$dose.range[2])})[rj$mlist[[rj$current.model]]]
      
      rj$nu[1, , 1] <- 
        sapply(X = inits.bi, FUN = function(x){
          rtnorm(n = 1, location = quantile(x = x, probs = 0.25, na.rm = TRUE, names = FALSE), scale = 2, 
                 L = rj.input$param$dose.range[1],
                 U = median(x, na.rm = TRUE))})[rj$mlist[[rj$current.model]]]
      
      rj$nu[1, , 2] <- 
        sapply(X = inits.bi, FUN = function(x){
          rtnorm(n = 1, location = quantile(x = x, probs = 0.75, na.rm = TRUE, names = FALSE), scale = 2, 
                 L = median(x, na.rm = TRUE),
                 U = rj.input$param$dose.range[2])})[rj$mlist[[rj$current.model]]]

      inits.tau <- rbind(sapply(X = inits.bi, FUN = function(x){sd(x[x<median(x, na.rm = TRUE)], na.rm = TRUE)}), sapply(X = inits.bi, FUN = function(x){sd(x[x>median(x, na.rm = TRUE)], na.rm = TRUE)}))
      
      inits.tau[is.na(inits.tau)] <- runif(n = sum(is.na(inits.tau)), 
                                           min = rj.input$config$priors["tau", 1],
                                           max = rj.input$config$priors["tau", 2])
      
      rj$tau[1, ] <- rowMeans(inits.tau) + 
        rtnorm(n = 1, location = 0, scale = 2, 
               L = rj.input$config$priors["tau", 1] - rowMeans(inits.tau),
               U = rj.input$config$priors["tau", 2] - rowMeans(inits.tau))

      rj$psi[1] <- rnorm(n = 1, mean = rj.input$config$priors["psi", 1], 
                         sd = rj.input$config$priors["psi", 2])
      
      rj$omega[1] <- runif(n = 1, min = rj.input$config$priors["omega", 1], 
                           max = rj.input$config$priors["omega", 2])
      rj$psi.i[1, ] <- rnorm(n = rj.input$whales$n, mean = rj$psi[1], sd = rj$omega[1])
      
      included.cov <- rj$include.covariates[1, ]
      
      if (sum(included.cov) == 0) {
        covnames <- "nil"
        cov.effects <- 0
      } else {
        included.cov <- included.cov[included.cov == 1]
        covvalues <- do.call(c, purrr::map(
          .x = names(included.cov),
          .f = ~ {qnorm(pnorm(q = rj[[.x]][1, ], 
                              mean = rj.input$config$priors[.x, 1], 
                              sd = rj.input$config$priors[.x, 2]))}))
        cov.effects <- colSums(t(do.call(cbind, rj.input$covariates$dummy[names(included.cov)])) * covvalues)
      }
      
      psi_ij <- rj$psi.i[1, rj.input$whales$id] + cov.effects

    rj$pi.ij[1, ] <- pnorm(q = psi_ij)  
    
    rj$mu.ij[1, , 1] <- unname(rtnorm(n = rj.input$trials$n, 
           location = rj$nu[1, , 1][rj.input$species$trials],
           scale = rj$tau[1, 1], 
           L = rj.input$param$dose.range[1], 
           U = rj$alpha[1, rj.input$species$trials]))

    rj$mu.ij[1, , 2] <- unname(rtnorm(n = rj.input$trials$n, 
          location = rj$nu[1, , 2][rj.input$species$trials],
          scale = rj$tau[1, 2], 
          L = rj$alpha[1, rj.input$species$trials],
          U = rj.input$param$dose.range[2]))

    rj$k.ij[1, ] <- (1 - rbinom(n = rj.input$trials$n, size = 1, prob = rj$pi.ij[1, ])) + 1
    rj$t.ij[1, ] <- ((2 - rj$k.ij[1, ]) * rj$mu.ij[1, , 1]) +
      ((1 - (2 - rj$k.ij[1, ])) * rj$mu.ij[1, , 2])
    }
    
    # Censoring
    if (!all(rj.input$obs$censored == 0)) {

      if(rj$phase[1] == 1){
      
      # Right-censored
      rj$t.ij[1, rj.input$obs$censored == 1] <-
        rtnorm(n = sum(rj.input$obs$censored == 1),
              location = mu.ij.config[rj.input$obs$censored == 1],
              scale = 30,
              L = rj.input$obs$Rc[rj.input$obs$censored == 1], 
              U = rj.input$param$dose.range[2])
      
      # Left-censored
      rj$t.ij[1, rj.input$obs$censored == -1] <-
        rtnorm(n = sum(rj.input$obs$censored == -1),
               location = mu.ij.config[rj.input$obs$censored == -1],
               scale = 30,
               L = rj.input$param$dose.range[1],
               U = rj.input$obs$Lc[rj.input$obs$censored == -1])
      
      } else {
       
        rj$k.ij[1, rj.input$obs$censored == 1] <- 2
        rj$k.ij[1, rj.input$obs$censored == -1] <- 1
        
        for (a in 1:rj.input$trials$n) {
          
          if (rj.input$obs$censored[a] == 1){
            if(rj$mu.ij[1, a, 2] < rj.input$obs$Rc[a]){
              rj$mu.ij[1, a, 2] <- rtnorm(
                n = 1,
                location = rj$nu[1, rj.input$species$trials[a], 2],
                scale = 30,
                L = rj.input$obs$Rc[a],
                U = rj.input$param$dose.range[2])
            }} 
          
          if (rj.input$obs$censored[a] == -1){
            if(rj$mu.ij[1, a, 1] > rj.input$obs$Lc[a]){ 
              rj$mu.ij[1, a, 1] <- rtnorm(
                n = 1,
                location = rj$nu[1, rj.input$species$trials[a], 1],
                scale = 30,
                L = rj.input$param$dose.range[1],
                U = rj.input$obs$Lc[a])
              
            }}

          rj$t.ij[1, ] <- (2 - rj$k.ij[1, ]) * rj$mu.ij[1, , 1] + 
            ((1 - (2 - rj$k.ij[1, ])) * rj$mu.ij[1, , 2])
          
        }
      }} # End if censoring

    if(any(is.na(rj$t.ij[1, ]))) warning("NA(s) returned as initial values for t.ij")
    if(any(is.infinite(rj$t.ij[1, ]))) warning("Inf returned as initial values for t.ij")
    
  } # End if do.update
  
  # Iteration
  mod.par <- tibble::tibble(param = c("alpha", "k.ij", "mu", "mu.i", "mu.ij",     
                           "nu", "omega", "phi", "pi.ij", "psi", "psi.i",
                           "sigma", "t.ij", "tau",   
                           # if(rj.input$covariates$n > 0) 
                             "covariates", # Need this term even when no covariates
                           rj.input$covariates$names)) %>% 
    dplyr::mutate(monophasic = c(0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 
                                 # if (rj.input$covariates$n > 0)
                                   1, 
                                 rep(1, rj.input$covariates$n)),
                    biphasic = c(1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
                                 # if (rj.input$covariates$n > 0) 
                                   1,  
                                 rep(1, rj.input$covariates$n)))

  if(!rj.input$config$function.select){
    if(rj.input$config$biphasic) 
      mod.par <- mod.par %>% dplyr::filter(biphasic == 1) else 
        mod.par <- mod.par %>% dplyr::filter(monophasic == 1)
  } else {
    mod.par <- mod.par %>% dplyr::filter(biphasic + monophasic > 0)
  }
  
  rj$iter <- rep(1, nrow(mod.par))
  names(rj$iter) <- mod.par$param
  
  # Finally, store the data
  rj$dat <- rj.input
  rj$config <- rj$dat$config
  rj$dat$config <- NULL
  rj$update <- do.update
  
  return(rj)
  
}