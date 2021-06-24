#' Diagnostic plots
#'
#' Generate trace, density, and autocorrelation plots from an \code{rjtrace} object.
#'
#' @param rj.obj rjMCMC trace object of class \code{rjtrace}.
#' @param param.name Parameter name(s). Defaults to \code{all}, which returns plots for all parameters in the model.
#' @param covariates.incl Logical. If \code{TRUE}, the trace is filtered to only retain posterior estimates obtained when the contextual covariates were included in the model. Only relevant when \code{covariate.select = TRUE} in \code{\link{configure_rjMCMC}}.
#' @param autocorr Logical. Whether to output chain autocorrelation plots.
#' @param individual Logical. If \code{TRUE}, separate density lines will be plotted for each chain. If \code{FALSE}, one density line will be plotted for all chains.
#' 
#' @details Adapted from Casey Youngflesh's \code{\link[MCMCvis]{MCMCtrace}}.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}} \code{\link{trace_rjMCMC}}
#' @examples
#' library(espresso)
#' 
#' # Import the example data
#' 
#' mydat <- read_data(file = NULL) 
#' 
#' # Import a real dataset with the sonar and range covariates, 
#' # excluding sperm whales and any other species with a sample size
#' # smaller than two
#' 
#' mydat <- read_data(file = "path/to/my/data.csv", 
#'                   exclude.species = "Sperm whale",
#'                   min.N = 2) 
#' 
#' @keywords brs rjmcmc 

plot.rjtrace <- function(rj.obj, 
                         param.name = "all", 
                         covariates.incl = TRUE,
                         autocorr = FALSE, 
                         individual = TRUE){
  
  #' ---------------------------------------------
  # Extract the trace
  #' ---------------------------------------------
  mcmc.trace <- rj.obj$trace
  
  #' ---------------------------------------------
  # Rename columns
  #' ---------------------------------------------
  for(nc in 1:length(mcmc.trace)){
    colnames(mcmc.trace[[nc]])[which(startsWith(colnames(mcmc.trace[[1]]), prefix = "mu"))] <- 
      paste0("mu (", rj.obj$dat$species$names, ")")}
  
  #' ---------------------------------------------
  # Generate plots
  #' ---------------------------------------------
  
  if(!is.null(param.name)){
    
    if(param.name %in% rj.obj$dat$covariates$names){
      
      if(covariates.incl){
        
        cov.cols <- which(grepl(pattern = param.name, x = colnames(mcmc.trace[[1]])))
        
        cov.trace <- mcmc.trace[, cov.cols]
        cov.trace <- purrr::map(.x = cov.trace, 
                                .f = ~tibble::as_tibble(.x) %>% 
                                  dplyr::filter(.[[2]] == 1) %>% 
                                  dplyr::select_at(vars(-contains("incl."))))
        
        cov.n <- purrr::map(.x = cov.trace, .f = ~nrow(.x)) %>% 
          unlist(.) %>% min(.)
        
        for(j in seq_along(cov.trace)){
          cov.trace[[j]] <- cov.trace[[j]][(nrow(cov.trace[[j]]) - cov.n + 1):nrow(cov.trace[[j]]), ] %>% 
            coda::as.mcmc(.)}
        
        cov.trace <- as.mcmc.list(cov.trace)
        
        MCMC_trace(cov.trace, 
                   iter = cov.n, 
                   pdf = FALSE, 
                   ind = individual,
                   params = "all")
        
        if(autocorr) bayesplot::mcmc_acf(x = cov.trace, pars = colnames(cov.trace[[1]]))
        
      } else {
        
        param.name <- colnames(mcmc.trace[[1]])[startsWith(colnames(mcmc.trace[[1]]), prefix = param.name)]
        
        MCMC_trace(mcmc.trace, 
                   iter = rj.obj$mcmc$n.iter, 
                   pdf = FALSE, 
                   ind = individual,
                   params = ifelse(is.null(param.name), "all", param.name))
        
        if(autocorr) bayesplot::mcmc_acf(x = mcmc.trace, pars = param.name)
        
      }
      
    } else {
      
      MCMC_trace(mcmc.trace, 
                 iter = rj.obj$mcmc$n.iter, 
                 pdf = FALSE, 
                 ind = individual,
                 params = ifelse(is.null(param.name), "all", param.name))
      
      if(autocorr) bayesplot::mcmc_acf(x = cov.trace, pars = param.name)
      
    }
    
  } else {
    
    if(covariates.incl){
      
      if(rj.obj$dat$covariates$n > 0){
        
        cov.cols <- purrr::map(.x = rj.obj$dat$covariates$names,
                               .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[1]])))) 
        
        cov.trace <- purrr::map(.x = cov.cols, .f = ~mcmc.trace[, .x])
        cov.trace <- purrr::map_depth(.x = cov.trace, .depth = 2, 
                                      .f = ~tibble::as_tibble(.x) %>% 
                                        dplyr::filter(.[[2]] == 1) %>% 
                                        dplyr::select_at(vars(-contains("incl."))))
        
        cov.n <- purrr::map_depth(.x = cov.trace, .depth = 2, .f = ~nrow(.x)) %>% 
          purrr::map_dbl(.x = ., .f = ~min(unlist(.x)))
        
        for(u in seq_len(rj.obj$dat$covariates$n)){
          for(j in seq_along(cov.trace[[u]])){
            cov.trace[[u]][[j]] <- cov.trace[[u]][[j]][(nrow(cov.trace[[u]][[j]]) - cov.n[u] + 1):nrow(cov.trace[[u]][[j]]), ] %>% 
              coda::as.mcmc(.)}}
        
        cov.trace <- purrr::map(.x = cov.trace, .f = ~as.mcmc.list(.x))
        
        mcmc.trace <- mcmc.trace[, -unlist(cov.cols)]
        
        MCMC_trace(mcmc.trace, 
                   iter = rj.obj$mcmc$n.iter, 
                   pdf = FALSE, 
                   ind = individual,
                   params = "all")
        
        
        purrr::walk2(.x = cov.trace,
                     .y = cov.n,
                     .f = ~MCMC_trace(.x, 
                                      iter = .y, 
                                      pdf = FALSE, 
                                      ind = individual,
                                      params = "all"))
        
        if(autocorr) bayesplot::mcmc_acf(x = mcmc.trace)
        if(autocorr) {
          for(j in seq_along(cov.trace)) temp <- bayesplot::mcmc_acf(x = cov.trace[[j]]); print(temp)
        }
        
        
      }
      
    } else {
      
      MCMC_trace(mcmc.trace, 
                 iter = rj.obj$mcmc$n.iter, 
                 pdf = FALSE, 
                 ind = individual,
                 params = ifelse(is.null(param.name), "all", param.name))
      
      if(autocorr) bayesplot::mcmc_acf(x = mcmc.trace, pars = ifelse(is.null(param.name), character(), param.name))
      
    }
    
    
  }
}