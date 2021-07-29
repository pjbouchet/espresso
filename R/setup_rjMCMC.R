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
  if(!"integer" %in% class(move.ratio) & !"list" %in% class(move.ratio)) stop("<move.ratio> must be either an integer vector or a list.")

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
  
  if(m < tot.iter){
    move.iters <- seq(from = m, to = tot.iter, by = m)
  } else {
    move.iters <- rep(0, tot.iter)
  }
  modelmoves[move.iters] <- rep(rep(c(1, 2, 3), move.ratio, length.out = length(move.iters)))
  modelmoves.tab <- table(modelmoves[(n.burn+1):tot.iter]) %>% 
    tibble::enframe() %>% 
    dplyr::rename(move = name, n = value) %>% 
    dplyr::mutate(move = dplyr::case_when(move == "0" ~ "split-merge",
                                          move == "1" ~ "data-driven (Type I)",
                                          move == "2" ~ "data-driven (Type II)",
                                          move == "3" ~ "random"))
  
  if(do.update){
    rj.input <- append(rj.input[[1]][["dat"]], rj.input[[1]]["config"])
  }
  
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
  
  rj <- append(rj, 
               list(
                 y.ij = matrix(data = 0, nrow = tot.iter, ncol = rj.input$trials$n),
                 t.ij = matrix(data = 0, nrow = tot.iter, ncol = rj.input$trials$n),
                 mu.i = matrix(data = 0, nrow = tot.iter, ncol = rj.input$whales$n),
                 mu = matrix(data = 0, nrow = tot.iter, ncol = rj.input$species$n),
                 phi = numeric(tot.iter),
                 sigma = numeric(tot.iter)))
  
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
  
  #' -----------------------------------------------
  # Assign column names
  #' -----------------------------------------------
  
  colnames(rj$y.ij) <- paste0("y", rep(seq_len(rj.input$whales$n), rj.input$trials$nper), ".", 
                              unlist(sapply(rj.input$trials$nper, seq_len)))
  colnames(rj$t.ij) <- paste0("t", rep(seq_len(rj.input$whales$n), rj.input$trials$nper), ".", 
                              unlist(sapply(rj.input$trials$nper, seq_len)))
  colnames(rj$mu.i) <- paste0("mu.", 1:rj.input$whales$n)
  colnames(rj$mu) <- paste0("mu.", 1:rj.input$species$n)
  
  rj$model <- numeric(tot.iter)
  rj$mlist <- rj.input$config$mlist
  
  if(is.null(rj.input$species$abbrev)) 
    rj$abbrev <- list(rj.input$species$abbrev) else rj$abbrev <- rj.input$species$abbrev
  
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
  

  rj$accept <- list(t.ij = numeric(rj.input$trials$n),
                    mu.i = numeric(rj.input$whales$n), 
                    mu = 0, phi = 0, sigma = 0)
  if(rj.input$covariates$n > 0) rj$accept <- append(rj$accept, numeric(rj.input$covariates$n))
  rj$accept <- append(rj$accept, numeric(4))
  
  if(rj.input$config$covariate.select) rj$accept <- append(rj$accept, numeric(1))
  names(rj$accept) <- unlist(list("t.ij", "mu.i", "mu", "phi", "sigma", 
                                  rj.input$covariates$names, 
                                  paste0("move.", c(0, 1, 2, 3)),
                                  ifelse(rj.input$config$covariate.select, "move.covariates", list(NULL))))
  
  
  #' -----------------------------------------------
  # Generate starting values
  #' -----------------------------------------------
  
  if(do.update){
    
    for(it in names(start.values)){
      if(it == "mlist") rj[[it]] <- start.values[[it]]
      if(is.null(dim(rj[[it]]))) rj[[it]][1] <- start.values[[it]][1] else 
        rj[[it]][1, ] <- start.values[[it]]
    }
    
  } else {
    
    # Model ID
    
    if (rj.input$config$model.select) {
      rj$current.model <- rj$model[1] <-
        sample(x = rj.input$config$clust[[1]]$model, size = 1, prob = rj.input$config$clust[[1]]$p)
    } else {
      # Because the model list only has length 1
      rj$current.model <- rj$model[1:tot.iter] <- names(rj.input$config$mlist) 
    }
    
    # Model parameters 
    
    rj$phi[1] <- rj.input$config$var[2] + 
      rtnorm(n = 1, location = 0, scale = 3, 
             L = rj.input$param$bounds["phi", 1] - rj.input$config$var[2],
             U = rj.input$param$bounds["phi", 2] - rj.input$config$var[2])
    
    rj$sigma[1] <- rj.input$config$var[1] + 
      rtnorm(n = 1, location = 0, scale = 3, 
             L = rj.input$param$bounds["sigma", 1] - rj.input$config$var[1],
             U = rj.input$param$bounds["sigma", 2] - rj.input$config$var[1])
    
    input.data <- tibble::tibble(species = rj.input$species$trials,
                   y = rj.input$obs$y_ij,
                   censored = rj.input$obs$censored,
                   rc = rj.input$obs$Rc,
                   lc = rj.input$obs$Lc)
    
    mu.start <- purrr::map_dbl(
      .x = rj.input$config$boot[[rj$current.model]],
      .f = ~ {
        
        sp.data <- input.data %>% dplyr::filter(species == .x)
        
        if(all(is.na(sp.data$y))){
          sp.data %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(y = ifelse(censored == 1,
                                     runif(n = 1, min = rc, max = rj.input$param$bounds["mu", 2]),
                                     runif(n = 1, min = rj.input$param$bounds["mu", 1], max = lc))) %>% 
            dplyr::ungroup() %>% 
            dplyr::pull(y) %>% 
            mean(., na.rm = TRUE)
        } else {
          mean(sp.data$y, na.rm = TRUE)
        }})
    
    rj$mu[1, ] <- mu.start + sapply(X = seq_len(n_groups(vec = mu.start)), FUN = function(x){
      rtnorm(n = 1, location = 0, scale = 5, 
             L = rj.input$param$bounds["mu", 1] - unique(mu.start)[x],
             U = rj.input$param$bounds["mu", 2] - unique(mu.start)[x])})[sapply(X = rj.input$config$boot[[rj$current.model]], FUN = function(x){which(unique(rj.input$config$boot[[rj$current.model]]) == x)})]
    
    rj$mu.i[1, ] <- rtnorm(n = rj.input$whales$n,
                           location = rj$mu[1, rj.input$species$id], 
                           scale = rj$phi[1],
                           L = rj.input$param$bounds["mu.i", 1], 
                           U = rj.input$param$bounds["mu.i", 2])
    
    mu.ij.config <- rj$mu.i[1, rj.input$whales$id]
    
    if(rj.input$covariates$n > 0) {
      
      if(rj.input$config$covariate.select) rj$include.covariates[1, ] <- 
          sample(x = 0:1, size = rj.input$covariates$n, replace = TRUE) 
      
      for (j in rj.input$covariates$names){
        v <- rnorm(n = ifelse(rj.input$covariates$fL[[j]]$nL == 0, 1, rj.input$covariates$fL[[j]]$nparam), mean = 0, sd = 3)
        if(rj.input$covariates$fL[[j]]$nL == 0){ rj[[j]][1, ] <- v  } else { rj[[j]][1, ] <- c(0, v)}}
      
    }
    
    # The scale needs to be large here to prevent numerical issues (returning Inf)
    rj$t.ij[1, ] <- rtnorm(n = rj.input$trials$n, 
                           location = mu.ij.config, 
                           scale = rj$sigma[1],
                           L = rj.input$param$bounds["t.ij", 1], 
                           U = rj.input$param$bounds["t.ij", 2])
    
    # Censoring
    if (!all(rj.input$obs$censored == 0)) {
      
      # Right-censored
      rj$t.ij[1, rj.input$obs$censored == 1] <-
        rtnorm(n = sum(rj.input$obs$censored == 1),
              location = mu.ij.config[rj.input$obs$censored == 1],
              scale = rj$sigma[1],
              L = rj.input$obs$Rc[rj.input$obs$censored == 1], 
              U = rj.input$param$bounds["t.ij", 2])
      
      # Left-censored
      rj$t.ij[1, rj.input$obs$censored == -1] <-
        rtnorm(n = sum(rj.input$obs$censored == -1),
               location = mu.ij.config[rj.input$obs$censored == -1],
               scale = rj$sigma[1],
               L = rj.input$param$bounds["t.ij", 1],
               U = rj.input$obs$Lc[rj.input$obs$censored == -1])
      
      
    }
    
    if(any(is.infinite(rj$t.ij[1, ]))) warning("Inf returned as initial values for t.ij")
    
    rj$y.ij[1, ] <- rj.input$obs$y_ij
    
    if (sum(!rj.input$obs$censored == 0) > 0) {
      
      # Censored
      rj$y.ij[1, !rj.input$obs$censored == 0] <-
        rnorm(n = sum(!rj.input$obs$censored == 0),
              mean = rj$t.ij[1, !rj.input$obs$censored == 0],
              sd = rj.input$obs$sd)
    }
    
    if(any(is.na(rj$y.ij[1, ]))) warning("NAs returned as initial values for y.ij")
    
  } # End if do.update
  
  # Finally, store the data
  rj$dat <- rj.input
  rj$config <- rj$dat$config
  rj$dat$config <- NULL
  rj$update <- do.update
  
  return(rj)
  
}