#' Diagnostic plots
#'
#' Generate trace, density, and autocorrelation plots from an \code{rjtrace} object.
#'
#' @export
#' @param rj.obj rjMCMC trace object of class \code{rjtrace}.
#' @param covariates.incl Logical. If \code{TRUE}, the trace is filtered to only retain posterior estimates obtained when the contextual covariates were included in the model. Only relevant when \code{covariate.select = TRUE} in \code{\link{configure_rjMCMC}}.
#' @inheritParams plot.gvs
#' 
#' @details Adapted from Casey Youngflesh's function \code{\link[MCMCvis]{MCMCtrace}}.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}} \code{\link{trace_rjMCMC}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Import the example data, excluding species with sample sizes < 5
#' # and considering the sonar covariate
#' mydat <- read_data(file = NULL, min.N = 5, covariates = "sonar") 
#' summary(mydat)
#' 
#' # Configure the sampler
#' mydat.config <- configure_rjMCMC(dat = mydat,
#'                                  model.select = TRUE,
#'                                  covariate.select = FALSE,
#'                                  proposal.mh = list(t.ij = 10, mu.i = 10, 
#'                                                     mu = 7, phi = 10, sigma = 10),
#'                                  proposal.rj = list(dd = 20, cov = 7),
#'                                  prior.covariates = c(0, 30),
#'                                  n.rep = 100)
#' summary(mydat.config)
#' 
#' # Run the reversible jump MCMC
#' rj <- run_rjMCMC(dat = mydat.config,
#'                  n.chains = 2,
#'                  n.burn = 100,
#'                  n.iter = 100,
#'                  do.update = FALSE)
#' 
#' # Burn and thin
#' rj.trace <- trace_rjMCMC(rj.dat = rj)
#' 
#' # Get density and trace plots
#' plot(rj.trace)
#' }
#' @keywords brs dose-response rjmcmc 

plot.rjtrace <- function(rj.obj, 
                         param.name = NULL, 
                         covariates.incl = FALSE,
                         autocorr = FALSE, 
                         individual = TRUE){
  
  # Redo this in simpler way.
  # 1 determine if covariates
  # run plot on all params except cov
  # for cov: run below with covariates.incl
  # then lastly run autocorr
  
  if(rj.obj$dat$covariates$n == 0) covariates.incl <- FALSE
  
  #' ---------------------------------------------
  # Extract the trace and rename columns
  #' ---------------------------------------------
  mcmc.trace <- rj.obj$trace
  
  for(nc in 1:length(mcmc.trace)){
    colnames(mcmc.trace[[nc]])[which(startsWith(colnames(mcmc.trace[[1]]), prefix = "mu"))] <- 
      paste0("mu (", rj.obj$dat$species$names, ")")}
  
  if(is.null(param.name)){
    
    if(covariates.incl){
      
      do.filter <- TRUE
      
    } else {
    
      do.filter <- FALSE
      MCMC_trace(mcmc.trace, 
               iter = rj.obj$mcmc$n.iter, 
               pdf = FALSE, 
               ind = individual,
               params = "all")
      
    if(autocorr) bayesplot::mcmc_acf(x = mcmc.trace, pars = bpars)
      
    }
    
    
  } else {
    
    pn <- purrr::map(.x = param.name, 
                     .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[1]])))) %>%
      unlist()
    
    mcmc.trace <- mcmc.trace[, pn, drop = FALSE]
    
    if(any(param.name %in% rj.obj$dat$covariates$names) & covariates.incl){
      
      do.filter <- TRUE
      
    } else {
      
      do.filter <- FALSE
      
      MCMC_trace(mcmc.trace, 
                 iter = rj.obj$mcmc$n.iter, 
                 pdf = FALSE, 
                 ind = individual,
                 params = "all")
      
    }
    
  }
  
  
  if(do.filter){
    
    
    cov.cols <- purrr::map(.x = param.name[param.name %in% rj.obj$dat$covariates$names],
                           .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[1]]))))
    
    cov.trace <- purrr::map(.x = cov.cols, .f = ~mcmc.trace[, .x, drop = FALSE])
    main.trace <- mcmc.trace[, seq_len(ncol(mcmc.trace[[1]]))[-unlist(cov.cols)], drop = FALSE]
    
    cov.trace <- purrr::map_depth(.x = cov.trace, 
                                  .depth = 2,
                                  .f = ~
                                    tibble::as_tibble(.x) %>% 
                                    dplyr::filter(dplyr::across(contains("incl."), ~. == 1)) %>% 
                                    dplyr::select_at(dplyr::vars(-contains("incl."))))
    
    cov.n <- purrr::map_depth(.x = cov.trace, 
                              .depth = 2,
                              .f = ~nrow(.x)) %>% 
      tibble::enframe() %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(minimum = min(unlist(value))) %>% 
      dplyr::ungroup() %>% 
      dplyr::pull(minimum)
    
    cov.trace.final <- cov.trace
    
    for(co in seq_along(cov.trace.final)){
      for(ch in seq_len(rj.obj$mcmc$n.chains)){
        cov.trace.final[[co]][[ch]] <- cov.trace.final[[co]][[ch]][(nrow(cov.trace.final[[co]][[ch]]) - cov.n[co] + 1):nrow(cov.trace.final[[co]][[ch]]), ] %>% coda::as.mcmc(.)}}
    
    cov.trace.final <- purrr::map(.x = cov.trace.final, .f = ~coda::as.mcmc.list(.x))
    
    MCMC_trace(main.trace, 
               iter = rj.obj$mcmc$n.iter, 
               pdf = FALSE, 
               ind = individual,
               params = "all")
    
    purrr::walk2(.x = cov.trace.final,
                 .y = cov.n,
                 .f = ~MCMC_trace(.x,  iter = .y, pdf = FALSE, ind = individual, params = "all"))
  }
  
  if(autocorr){
    
    if(is.null(param.name)) bpars <- character() else bpars <- purrr::map(.x = param.name, .f = ~colnames(mcmc.trace[[1]])[grepl(pattern = .x, x = colnames(mcmc.trace[[1]]))]) %>% unlist()
    bayesplot::mcmc_acf(x = mcmc.trace, pars = bpars)
    }
  
}