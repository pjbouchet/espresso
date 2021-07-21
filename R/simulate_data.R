#' Generate a simulated dataset
#'
#' Simulate cetacean responses to Navy sonar exposure, as in behavioural response studies (BRSs).
#'
#' @export
#' @param n.species Number of species.
#' @param n.whales Number of individuals. Can be specified as a vector of length \code{n.species}, allowing for uneven samples sizes across species, or as an integer multiple of \code{n.species} for equal sample sizes across species.
#' @param max.trials Maximum number of exposures per individual. Each animal is exposed \code{x} times, with \code{x} taken as a random draw from a Uniform distribution between 1 and \code{max.trials}.
#' @param covariates Contextual covariates and their associated coefficients. Must be supplied as a named list, in which the baseline levels of factor covariates are given a coefficient of zero.
#' @param mu Mean response threshold(s) for each species.
#' @param phi Between-whale variance in response thresholds.
#' @param sigma Within-whale between-exposure variance in response thresholds.
#' @param Rc Right-censoring interval. Values of the maximum realised dose for each exposure are generated as random draws from a Uniform distribution within these bounds, such that \code{max.dose <- runif(n = n.trials, min = Rc[1], max = Rc[2])}.
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
#'   \code{obs} \tab Observations, including right-censoring cutoffs.  \cr  
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
#' mydat <- simulate_data(n.species = 2, 
#'                        n.whales = 16, 
#'                        max.trials = 3, 
#'                        covariates = list(exposed = c(0, 5), range = 0.5),
#'                        mu = c(125, 158), 
#'                        phi = 20, 
#'                        sigma = 20, 
#'                        Rc = c(210, 211), 
#'                        seed = 58697)
#' }
#' @keywords brs gvs rjmcmc dose-response                      

simulate_data <- function(n.species = 2, 
                          n.whales = 20, 
                          max.trials = 3,
                          covariates = NULL,
                          mu = NULL, 
                          phi = NULL,
                          sigma = NULL,
                          range.dB = c(60, 215),
                          range.var = c(0, 45),
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
    # which.cov <- allcov[which(allcov %in% which.cov)]
    n.covariates <- length(covariates)
  } else {
    covariate.names <- NULL
    which.cov <- n.covariates <- 0
  }
  
  if (is.null(mu) | is.null(phi) | is.null(sigma)) stop("Some parameters are missing.")
  if (max.trials < 1) stop("One or more trials required.")
  if (length(n.whales) == 1 & !sum(n.whales) %% n.species == 0) stop("Arguments must match.")
  if (length(n.whales) > 1 & !length(n.whales) == length(mu)) stop("Arguments must match.")
  if (any(mu < range.dB[1]) | any(mu > range.dB[2]) |
      any(phi < range.var[1]) | any(phi > range.var[2]) |
      any(sigma < range.var[1]) | any(sigma > range.var[2])) {
    stop("Parameter out of bounds.")
  }
  if (!length(Rc) == 2) stop("Rc must be a vector of length 2.")
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
    covariate.types <- covariate.types[covariate.names] %>% unlist()}
  
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
  
  mu.i <- rtnorm(n = n.whales, location = rep(mu, n.per.species), 
                 scale = phi, L = lower.bound$mu.i, U = upper.bound$mu.i)
  
  # Contextual covariates: 
  # History of exposure (factor, 2 levels)
  # Sonar signal frequency (factor, 2 levels)
  # Behavioural mode (factor, 3 levels)
  # Range (continuous)
  
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
  
  mu.ij <- mu.i[whale.id]
  
  if(n.covariates > 0) {
    
    if(any(covariate.names == "exposed")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "exposed"])), FUN = function(fl){dummy.cov[["exposed"]][, fl] * covariates[["exposed"]][fl]}))
    
    if(any(covariate.names == "sonar")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "sonar"])), FUN = function(fl){dummy.cov[["sonar"]][, fl] * covariates[["sonar"]][fl]}))
    
    if(any(covariate.names == "behaviour")) mu.ij <-  mu.ij + Reduce(f = "+", lapply(X = seq_len(nlevels(covariates.df[, "behaviour"])), FUN = function(fl){dummy.cov[["behaviour"]][, fl] * covariates[["behaviour"]][fl]}))
    
    if(any(covariate.names == "range")) mu.ij <-  mu.ij + covariates.df[, "range"] * covariates[["range"]] 
    
  }
  
  t.ij <- rtnorm(n = n.trials, location = mu.ij, scale = sigma, L = lower.bound$mu, U = upper.bound$mu)
  y_ij <- rnorm(n = n.trials, mean = t.ij, sd = obs.sd)
  
  #' ---------------------------------------------
  # Right-censoring
  #' ---------------------------------------------
  
  max.dose <- runif(n = n.trials, min = Rc[1], Rc[2])
  is.censored <- ifelse(y_ij > max.dose, 1, 0)
  y_ij[is.censored == 1] <- NA
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
               Rc = max.dose,
               Rc.range = Rc,
               prec = measurement.precision,
               sd = obs.sd),
    # Parameters
    param = list(sim = simulation,
                 bounds = param.bounds,
                 dose.range = range.dB,
                 mu = mu,
                 phi = phi,
                 sigma = sigma))
  
  if(n.covariates > 0) {
    
    res$covariates <- append(res$covariates, 
                             list(coefs = covariates,
                                  values = covariates.list,
                                  index = I.covariates, 
                                  df = covariates.df,
                                  dummy = dummy.cov))
    
    names(res$covariates$coefs) <- covariate.names
    
    res$covariates$fL <- sapply(X = covariate.names, 
                                FUN = function(x) factor_levels(covname = x, dat = res$covariates$df), 
                                simplify = FALSE, USE.NAMES = TRUE)
    
  }
  
  if(verbose){
  if(any(n.per.species < 5)) 
    warning("N < 5 for some species. Species with low sample sizes may need to be grouped a priori.")}
  
  class(res) <- c("brsdata", class(res))
  return(res) 
  
}