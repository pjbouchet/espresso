#' Generate a simulated dataset
#'
#' Simulate cetacean responses to Navy sonar exposure, as in behavioural response studies (BRSs).
#'
#' @export
#' @param biphasic Logical. If TRUE, will simulate data from a biphasic dose-response functional form.
#' @param n.species Number of species.
#' @param n.whales Number of individuals. Can be specified as a vector of length \code{n.species}, allowing for uneven samples sizes across species, or as an integer multiple of \code{n.species} for equal sample sizes across species.
#' @param max.trials Maximum number of exposures per individual. Each animal is exposed \code{x} times, with \code{x} taken as a random draw from a Uniform distribution between 1 and \code{max.trials}.
#' @param covariates Contextual covariates and their associated coefficients. Must be supplied as a named list, in which the baseline levels of factor covariates are given a coefficient of zero.
#' @param mu Mean response threshold(s) for each species.
#' @param phi Between-whale variance in response thresholds.
#' @param nu Mean response thresholds from the context-dependent and dose-dependent mixtures in a biphasic model. Must be a list of length two, where each element is a vector of length \code{n.species}.
#' @param tau Scale parameters associated with \code{nu}. Must be a vector of length two.
#' @param psi Probability of exhibiting a context-dependent response, expressed on a probit scale.
#' @param omega Variance in \code{psi}.
#' @param alpha Upper bound on the context-dependent response threshold, which corresponds to the lower bound of the dose-dependent threshold. Must be a vector of length \code{n.species}.
#' @param sigma Within-whale between-exposure variance in response thresholds.
#' @param Lc Left-censoring interval. Values of the minimum realised dose for each exposure are generated as random draws from a Uniform distribution within the bounds defined by `\code{Lc}.
#' @param Rc Right-censoring interval. Values of the maximum realised dose for each exposure are generated as random draws from a Uniform distribution within the bounds defined by `\code{Rc}.
#' @param seed Random seed (for reproducible results). 
#' 
#' @inheritParams read_data
#' 
#' @return A list object of class \code{brsdata}, with the following elements:
#' 
#' \tabular{ll}{
#'   \code{dat} \tab Output dataset, after processing. \cr
#'   \code{species} \tab Species data, including species names, groups, sample sizes etc. \cr
#'   \code{whales} \tab Individual data, including animal IDs.  \cr
#'   \code{trials} \tab Exposure data, including exposure IDs. \cr
#'   \code{covariates} \tab Covariate data, including dummy coding, sonar groupings, factor levels etc. \cr
#'   \code{obs} \tab Observations, including left- and right-censoring cutoffs.  \cr  
#'   \code{param} \tab General parameters. \cr
#'  }
#'  
#' @author Phil J. Bouchet
#' @seealso \code{\link{read_data}} \code{\link{summary.brsdata}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Simulate data for two species 
#' # (no censoring, monophasic functional form)
#' mydat <- simulate_data(biphasic = FALSE,
#'                        n.species = 2,
#'                        n.whales = 16,
#'                        max.trials = 3,
#'                        covariates = list(exposed = c(0, 5), range = 0.5),
#'                        mu = c(125, 142),
#'                        phi = 20,
#'                        sigma = 20,
#'                        Lc = c(60, 65),
#'                        Rc = c(214, 215),
#'                        seed = 58697)
#'                        
#' # Simulate data for three species 
#' # (right- and left-censoring, biphasic functional form)
#' mydat <- simulate_data(biphasic = TRUE,
#'                        n.species = 3,
#'                        n.whales = c(10, 5, 8),
#'                        max.trials = 3,
#'                        alpha = c(125, 140, 139),
#'                        nu = list(c(98, 110, 105), c(149, 165, 152)),
#'                        tau = c(20, 20),
#'                        psi = 0.5,
#'                        omega = 1,
#'                        Lc = c(60, 70),
#'                        Rc = c(165, 195),
#'                        seed = 58697)
#'                        
#' }
#' @keywords brs gvs rjmcmc dose-response                      

simulate_data <- function(biphasic = FALSE,
                          n.species = 2, 
                          n.whales = 20, 
                          max.trials = 3,
                          covariates = NULL,
                          mu = NULL, 
                          phi = NULL,
                          sigma = NULL,
                          nu = NULL,
                          tau = NULL,
                          psi = 0.5,
                          omega = 1,
                          alpha = NULL,
                          range.dB = c(60, 215),
                          range.var = c(0, 45),
                          range.biphasic = list(alpha = c(110, 160), omega = c(0, 5), tau = c(0, 20)),
                          Lc = c(65, 75),
                          Rc = c(190, 195),
                          obs.sd = 2.5,
                          verbose = TRUE,
                          seed = NULL){
  
  
  # ---------------------------------------------
  # This function is for simulated data
  # ---------------------------------------------
  
  simulation <- TRUE
  
  # Set the random seed
  if(!is.null(seed)) set.seed(seed) 
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if (!is.null(covariates)) {
    covariate.names <- names(covariates)
    n.covariates <- length(covariates)
  } else {
    covariate.names <- NULL
    which.cov <- n.covariates <- 0
  }
  
  if (max.trials < 1) stop("One or more trials required.")
  if (length(n.whales) == 1 & !sum(n.whales) %% n.species == 0) stop("Arguments must match.")

  if(!biphasic) {  
    if (length(n.whales) > 1 & !length(n.whales) == length(mu)) stop("Arguments must match.")
    if (is.null(mu) | is.null(phi) | is.null(sigma)) stop("Some parameters are missing.") 
    if (any(mu < range.dB[1]) | any(mu > range.dB[2]) |
        any(phi < range.var[1]) | any(phi > range.var[2]) |
        any(sigma < range.var[1]) | any(sigma > range.var[2])) {
      stop("Parameter out of bounds.")
    }
  }
  if (!length(Rc) == 2) stop("Rc must be a vector of length 2.")
  if(any(Rc > range.dB[2])) stop("Rc out of range.")
  if(any(Rc < range.dB[1])) stop("Rc out of range.")
  if (n.covariates > 4) stop("Only four covariates are allowed: <exposed>, <sonar>, <behaviour>, and <range>.")
  
  if(!is.null(covariates) & !"list" %in% class(covariates)) stop("Erroneous format for input covariates.")
  
  if (n.covariates > 0) {
    covariate.L <- list(exposed = 2, sonar = 2, behaviour = 3, range = 1)
    for(cc in covariate.names){
      if (!length(covariates[[cc]]) == covariate.L[[cc]])
        stop(paste0("Wrong parameter input for <", cc, "> covariate."))
      if (cc %in% c("sonar", "behaviour") & covariates[[cc]][1] > 0)
        stop(paste0("Wrong parameter input for <", cc, "> covariate."))
    }
  }
  
  if(biphasic){
    if (n.species > 1 & !n.species == length(nu[[1]])) stop("Arguments must match.")
    if (n.species > 1 & !n.species == length(nu[[2]])) stop("Arguments must match.")
    if(any(nu[[1]] > nu[[2]])) stop(paste0("Wrong parameter input for <nu>."))
    if(any(tau < 0)) stop("Non-positive integer in variance parameter <tau>")
    if(any(alpha < nu[[1]]) | any(alpha > nu[[2]])) stop("Alpha must be greater than nu[1] and lesser than nu[2].")
    if(!length(alpha) == n.species) stop(paste0("Wrong parameter input for <alpha>."))
    if(sum(purrr::map2_lgl(.x = alpha, .y = nu[[1]], .f = ~ .x < .y)) > 0) stop("<alpha> and <nu> do not match.")
    if(sum(purrr::map2_lgl(.x = alpha, .y = nu[[2]], .f = ~ .x > .y)) > 0) stop("<alpha> and <nu> do not match.")
    if(any(is.null(nu), is.null(tau), is.null(psi), is.null(omega), is.null(alpha))) stop("Missing parameters.")
    if(length(psi) > 1) stop(paste0("<psi> must be of length 1."))
    if(length(omega) > 1) stop(paste0("<omega> must be of length 1."))
  }
  
  #' ---------------------------------------------
  # Parameterisation
  #' ---------------------------------------------
  
  # Number of whales
  if (length(n.whales) > 1) {
    n.per.species <- n.whales
    n.whales <- sum(n.whales)
  } else {
    n.per.species <- rep(n.whales / n.species, n.species)
  }
  
  # Species names
  sp.names <- as.character(seq_len(n.species))
  
  # Species IDs
  species.id <- sapply(X = 1:n.species, FUN = function(x) rep(x, each = n.per.species[x])) %>%
    unlist() %>%
    as.vector()
  
  # Covariates
  if(n.covariates > 0){
    covariate.types <- list(exposed = "f", sonar = "f", behaviour = "f", range = "d")
    covariate.types <- covariate.types[covariate.names] %>% unlist()
  }
  
  # Lower bounds for model parameters
  lower.bound <- list(
    y.ij = range.dB[1],
    t.ij = range.dB[1],
    mu = range.dB[1], 
    mu.i = range.dB[1], 
    phi = range.var[1], 
    sigma = range.var[1])
  
  # Upper bounds for model parameters
  upper.bound <- list(
    y.ij = range.dB[2],
    t.ij = range.dB[2],
    mu = range.dB[2],
    mu.i = range.dB[2],
    phi = range.var[2],
    sigma = range.var[2])
  
  if(biphasic){
    lower.bound$nu <- range.dB[1]
    lower.bound$alpha <- range.biphasic$alpha[1]
    lower.bound$tau <- range.biphasic$tau[1]
    lower.bound$omega <- range.biphasic$omega[1]
    
    upper.bound$nu <- range.dB[2]
    upper.bound$alpha <- range.biphasic$alpha[2]
    upper.bound$tau <- range.biphasic$tau[2]
    upper.bound$omega <- range.biphasic$omega[2]
  }
  
  for(j in seq_len(n.covariates)) lower.bound[covariate.names[j]] <- -Inf
  for(j in seq_len(n.covariates)) upper.bound[covariate.names[j]] <- Inf
  
  param.bounds <- cbind(purrr::map_df(.x = lower.bound, .f = ~.x) %>% t(.),
                        purrr::map_df(.x = upper.bound, .f = ~.x) %>% t(.))
  
  # Measurement precision (inverse of variance)
  measurement.precision <- 1 / (obs.sd^2)
  
  # Trials and individual IDs
  n.trials.per.whale <- sample(x = 1:max.trials, size = n.whales, replace = TRUE)
  n.trials <- sum(n.trials.per.whale)
  whale.id <- purrr::map(.x = seq_len(n.whales), 
                         .f = ~rep(.x, each = n.trials.per.whale[.x])) %>% do.call(c, .)
  
  #' ---------------------------------------------
  # Simulate data
  #' ---------------------------------------------
  
  if(n.covariates > 0){
    
    covariates.list <- list(
      exposed = purrr::map(
        .x = n.trials.per.whale,
        .f = ~ c(0, rep(1, each = .x - 1))
      ) %>% do.call(c, .),
      sonar = sample(x = c("LFAS", "MFAS"), size = n.trials, replace = TRUE),
      behaviour = sapply(seq_len(n.trials), sample,
                         x = c("Feed", "Migrate", "Rest"),
                         size = 1),
      range = runif(n = n.trials, min = 0, max = 20)
    )
    
    covariates.df <- as.data.frame(do.call(cbind, covariates.list[covariate.names]))
    colnames(covariates.df) <- covariate.names
    
    for(ii in 1:length(covariate.names)){
      if(covariate.types[ii] == "f") 
        covariates.df[, covariate.names[ii]] <- as.factor(covariates.df[, covariate.names[ii]])
      if(covariate.types[ii] %in% c("d", "i")) 
        covariates.df[, covariate.names[ii]] <- as.numeric(covariates.df[, covariate.names[ii]])
    }
    
    # Dummy coding
    dummy.cov <- purrr::map(.x = seq_len(n.covariates),
                            .f = ~{
                              if(covariate.types[.x] == "f"){
                                fastDummies::dummy_cols(.data = covariates.df, 
                                                        select_columns = covariate.names[.x],
                                                        remove_first_dummy = FALSE, 
                                                        remove_selected_columns = TRUE) %>% 
                                  dplyr::select(-tidyselect::any_of(covariate.names))
                                
                              } else { covariates.df[, .x, drop = FALSE] }}) %>% 
      purrr::set_names(., covariate.names)
    
    # Factor levels
    fL <- sapply(covariate.names, 
                 FUN = function(x) factor_levels(covname = x, dat = covariates.df), 
                 simplify = FALSE, USE.NAMES = TRUE)
    
    I.covariates <- purrr::map(.x = covariate.names, 
                               .f = ~ dummy.cov[[.x]][, fL[[.x]]$index, drop = FALSE])
    
    I.covariates <- purrr::map(.x = I.covariates,
                               .f = ~apply(.x, 2, list) %>% purrr::flatten(.)) %>% 
      purrr::flatten(.)
    
  }
  
  # Monophasic functional form
  
  if(!biphasic){
    
    mu.i <- rtnorm(n = n.whales, location = rep(mu, n.per.species), 
                   scale = phi, L = lower.bound$mu.i, U = upper.bound$mu.i)
    
    # Contextual covariates: 
    # History of exposure (factor, 2 levels)
    # Sonar signal frequency (factor, 2 levels)
    # Behavioural mode (factor, 3 levels)
    # Range (continuous)
  
    mu.ij <- mu.i[whale.id]
    
    if(n.covariates > 0) {
      
      if(any(covariate.names == "exposed")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "exposed"])), FUN = function(fl){dummy.cov[["exposed"]][, fl] * covariates[["exposed"]][fl]}))
      if(any(covariate.names == "sonar")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "sonar"])), FUN = function(fl){dummy.cov[["sonar"]][, fl] * covariates[["sonar"]][fl]}))
      if(any(covariate.names == "behaviour")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "behaviour"])), FUN = function(fl){dummy.cov[["behaviour"]][, fl] * covariates[["behaviour"]][fl]}))
      if(any(covariate.names == "range")) mu.ij <-  mu.ij + covariates.df[, "range"] * covariates[["range"]] 
      
    }
    
    t.ij <- rtnorm(n = n.trials, location = mu.ij, scale = sigma, L = lower.bound$mu, U = upper.bound$mu)
    y_ij <- rnorm(n = n.trials, mean = t.ij, sd = obs.sd)
    
  } else { # Biphasic functional form
    
    psi_i <- rnorm(n = n.whales, mean = psi, sd = omega) 
    psi_ij <- psi_i[whale.id]
    
    if(n.covariates > 0){
      
      # Translate to probit scale
      covariates.biphasic <- 
        purrr::map(.x = covariate.names,
                   .f = ~{sapply(X = seq_along(covariates[[.x]]),
                          FUN = function(ind){if(ind %in% fL[[.x]]$index){
                           qnorm(pnorm(covariates[[.x]][ind], 
                                       mean = priors[[.x]][1], sd = priors[[.x]][2]), mean = 0, sd = 1) 
                                   } else {
                                     0 
                                   }} )}) %>% purrr::set_names(., covariate.names)
      
      if(any(covariate.names == "exposed")) psi_ij <-  psi_ij + covariates.df[, "exposed"] * covariates.biphasic[["exposed"]] 
      
      if(any(covariate.names == "sonar")) psi_ij <-  psi_ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "sonar"])), FUN = function(fl){dummy.cov[["sonar"]][, fl] * covariates.biphasic[["sonar"]][fl]}))
      
      if(any(covariate.names == "behaviour")) psi_ij <-  psi_ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "behaviour"])), FUN = function(fl){dummy.cov[["behaviour"]][, fl] * covariates.biphasic[["behaviour"]][fl]}))
      
      if(any(covariate.names == "range")) psi_ij <-  psi_ij + covariates.df[, "range"] * covariates.biphasic[["range"]] 
      
      pi_ij <- pnorm(q = psi_ij)
      
    } else {
      
      # No covariates
      pi_ij <- pnorm(q = psi_i[whale.id])  
      
    }

    # Choose which part of the dose-response curve each whale responds to, with prob pi_ij
    k_ij <- (1 - rbinom(n = n.trials, size = 1, prob = pi_ij)) + 1
    
    # Simulate for each trial the value of the lower and upper mixture components
    mu_ij <- matrix(data = 0, nrow = n.trials, ncol = 2)
    mu_ij[, 1] <- rtnorm(n = n.trials, location =  nu[[1]][species.id[whale.id]], scale = tau[1], L = param.bounds["nu", 1], U = alpha)
    mu_ij[, 2] <- rtnorm(n = n.trials, location =  nu[[2]][species.id[whale.id]], scale = tau[2], L = alpha, U = param.bounds["nu", 2])
    
    t_ij <- numeric(n.trials)
    for (u in 1:n.trials) t_ij[u] <- mu_ij[u, k_ij[u]]
    y_ij <- rnorm(n = n.trials, mean = t_ij, sd = obs.sd)
    
  }
  
  #' ---------------------------------------------
  # Left- and right-censoring
  #' ---------------------------------------------
  
  min.dose <- runif(n = n.trials, min = Lc[1], Lc[2])
  max.dose <- runif(n = n.trials, min = Rc[1], Rc[2])

  is.lcensored <- ifelse(y_ij < min.dose, -1, 0)
  is.rcensored <- ifelse(y_ij > max.dose, 1, 0)
  is.censored <- is.lcensored + is.rcensored
  
  y_ij[!is.censored == 0] <- NA
  
  min.dose[is.censored == 0] <- upper.bound$mu
  max.dose[is.censored == 0] <- lower.bound$mu

  #' ---------------------------------------------
  # Species summary
  #' ---------------------------------------------
  
  species.trials <- sapply(X = seq_along(species.id), 
                           FUN = function(x) rep(species.id[x], n.trials.per.whale[x])) %>% do.call(c, .)
  
  brsdat <- data.frame(cbind(species.trials, y_ij))
  names(brsdat) <- c("species", "spl")
  
  suppressWarnings(species.summary <- brsdat %>% 
                     dplyr::group_by(species) %>% 
                     dplyr::summarise(N_trials = dplyr::n(), 
                                      censored = sum(is.na(spl)), 
                                      mean = mean(spl, na.rm = TRUE),
                                      min = min(spl, na.rm = TRUE),
                                      max = max(spl, na.rm = TRUE), .groups = "keep") %>% 
                     dplyr::ungroup())
  
  species.summary$N_ind <- n.per.species[species.summary$species] 
  species.summary <- species.summary %>% dplyr::select(species, N_ind, N_trials, mean, min, max)
  
  
  #' ---------------------------------------------
  # Store results in list
  #' ---------------------------------------------
  
  res <- list(
    
    biphasic = biphasic,
    
    # Species
    species = list(names = sp.names,
                   groups = NULL,
                   abbrev = NULL,
                   n = n.species,
                   nper = n.per.species,
                   id = species.id,
                   summary = species.summary,
                   trials = species.trials),
    # Individuals
    whales = list(n = n.whales,
                  id = whale.id),
    # Exposures
    trials = list(n = n.trials,
                  nper = n.trials.per.whale),
    # Contextual factors
    covariates = list(n = n.covariates,
                      names = covariate.names),
    # Observations
    obs = list(y_ij = y_ij,
               censored = is.censored,
               Lc = min.dose,
               Lc.range = Lc,
               Rc = max.dose,
               Rc.range = Rc,
               prec = measurement.precision,
               sd = obs.sd),
    
    param = list(sim = simulation,
                 bounds = param.bounds,
                 dose.range = range.dB))
  
  if(!biphasic){
    res$param <- append(res$param, list(mu = mu, phi = phi, sigma = sigma))
  } else {
    res$param <- append(res$param, 
                  list(nu = nu, tau = tau, psi = psi, 
                  omega = omega, alpha = alpha))
  }
  
  if(n.covariates > 0) {
    
    res$covariates <- append(res$covariates, 
                             list(coefs = covariates,
                                  values = covariates.list,
                                  index = I.covariates, 
                                  df = covariates.df,
                                  dummy = dummy.cov))
    
    names(res$covariates$coefs) <- covariate.names
    res$covariates$fL <- fL
    if(biphasic) res$covariates$biphasic <- covariates.biphasic
    # sapply(X = covariate.names, 
    #                             FUN = function(x) factor_levels(covname = x, dat = res$covariates$df), 
    #                             simplify = FALSE, USE.NAMES = TRUE)
    
  }
  
  if(verbose){
  if(any(n.per.species < 5)) 
    warning("N < 5 for some species. Species with low sample sizes may need to be grouped a priori.")}
  
  class(res) <- c("brsdata", class(res))
  return(res) 
  
}