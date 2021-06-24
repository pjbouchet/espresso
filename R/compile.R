#' Dose-response curves
#'
#' Compute dose-response functions from a fitted rjMCMC model.
#'
#' @param rj.object Input rjMCMC object of class \code{rjtrace}, as returned by \code{\link{trace_rjMCMC}}.
#' @param by.model Logical. If \code{TRUE}, the functions subsets posterior parameter estimates to produce separate dose-response curves for each model.
#' @param covariate Covariate name. This argument can be used to generate dose-response curves for specific contextual covariates, conditioned on the species (group) given in \code{species}.
#' @param covariate.values A vector of values for which dose-response curves are required. Only valid for continuous covariates.
#' @param species Species name. 
#' 
#' @return A list object of class \code{dose_response}.
#'  
#'  
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.rjdata}}
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
#' @export
#' @keywords brs rjmcmc 

compile_rjMCMC <- function(rj.object, 
                           by.model = FALSE, 
                           covariate = NULL, 
                           covariate.values = NULL,
                           species = NULL){
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if(!"rjtrace" %in% class(rj.object)) stop("Input must be an object of class <rj_trace>.")
  # if(!"rjdata" %in% class(rj.object$dat)) stop("Input must be an object of class <rj_data>.")
  if(!is.null(covariate)){ if(!covariate %in% rj.object$dat$covariates$names) stop("Unrecognised covariate.")}
  if(length(covariate) > 1) stop("<covariate> must be of length 1.")
  if(!all(species %in% rj.object$dat$species$names)) stop("Unrecognised species.")
  if(is.null(species) & !is.null(covariate)) stop("<species> cannot be set to NULL when a covariate is specified.")
  if(length(species) > 1 & !is.null(covariate)) stop("Max of 1 species (group) allowed")
  
  if(rj.object$dat$covariates$n == 0) covariate <- NULL
  
  if(!is.null(covariate)){
    
    if(is.null(covariate.values) & 
       rj.object$dat$config$fL[[covariate]]$nL == 0) stop("Must provide covariate values")
    
    if(is.null(covariate.values) & 
       rj.object$dat$config$fL[[covariate]]$nL > 0) stop("<covariate.values> ignored for factor covariates.")
    
    if(any(covariate.values < range(rj.object$dat$covariates$values[[covariate]])[1]) | 
       any(covariate.values > range(rj.object$dat$covariates$values[[covariate]])[2])) stop("Covariate values out of range")
  }
  
  #' ---------------------------------------------
  # Define dose range and quantiles
  #' ---------------------------------------------

  dose.range <- seq(rj.object$dat$param$bounds["mu", 1], rj.object$dat$param$bounds["mu", 2], length = 100)
  quants <- seq(95, 1, by = -5)
  
  #' ---------------------------------------------
  # Extract trace
  #' ---------------------------------------------
  
  mcmc.trace <- do.call(rbind, rj.object$trace) %>% tibble::as_tibble()
  
  if(!is.null(covariate.values)){
    
    for(tt in covariate.values){
      mcmc.trace <- mcmc.trace %>% 
        dplyr::mutate(!!paste0(covariate, "_", tt) := .data[[covariate]] *  tt)
    }
    mcmc.trace <- mcmc.trace %>% dplyr::select(-which(names(mcmc.trace) == covariate))
  }
  
  #' ---------------------------------------------
  # Split trace by model ID
  #' ---------------------------------------------
  
  if(by.model) mcmc.list <- purrr::map(.x = unique(mcmc.trace$model_ID), 
                                       .f = ~dplyr::filter(mcmc.trace, model_ID == .x)) else 
                                         mcmc.list <- list(mcmc.trace)
  
  #' ---------------------------------------------
  # Compute dose-response curves
  #' ---------------------------------------------
  
  pb <- progress::progress_bar$new(format = "Computing [:bar] :percent eta: :eta ",
                                   total = length(mcmc.list), 
                                   width = getOption("width") + 5)
  
  res <- purrr::map(.x = seq_along(mcmc.list), 
                    .f = function(mcl){
                      
                      if(by.model) pb$tick() # Launch progress bar
                      
                      # Identify columns corresponding to species means
                      if (!is.null(species)) {
                        mu.index <- which(colnames(mcmc.list[[mcl]]) %in%
                                            paste0("mu.", which(rj.object$dat$species$names %in% species)))
                        sp.names <- species
                      } else {
                        mu.index <- which(startsWith(x = colnames(mcmc.list[[mcl]]), prefix = "mu"))
                        sp.names <- unlist(rj.object$dat$species$names)
                      }
                      
                      dr.raw <- lapply(X = mu.index, 
                                       FUN = function(mu.x) {
                                         
                                         if(!is.null(covariate)){
                                           
                                           cov.index <- which(startsWith(x = colnames(mcmc.list[[mcl]]),
                                                                         prefix = covariate))
                                           
                                           if(rj.object$config$fL[[covariate]]$nL > 0) 
                                             cov.values <- list(rep(0, nrow(mcmc.list[[mcl]]))) else 
                                             cov.values <- list()
                                           
                                           cov.values <- append(cov.values,
                                                                lapply(X = cov.index, 
                                                                       FUN = function(cc) dplyr::pull(mcmc.list[[mcl]][, cc])))
                                           
                                           dr.mean <- purrr::map(.x = cov.values,
                                                                 .f = ~dplyr::pull(mcmc.list[[mcl]][, mu.x]) + .x) %>% 
                                             purrr::set_names(x = ., nm = rj.object$config$fL[[covariate]]$Lnames)
                                           
                                           
                                         } else {
                                           
                                           dr.mean <- list(dplyr::pull(mcmc.list[[mcl]][, mu.x]))
                                         }
                                         
                                         purrr::map(.x = dr.mean,
                                                    .f = ~truncnorm::ptruncnorm(
                                                      q = rep(dose.range, each = nrow(mcmc.list[[mcl]])),
                                                      a = rj.object$dat$param$bounds["mu", 1],
                                                      b = rj.object$dat$param$bounds["mu", 2],
                                                      mean = .x,
                                                      sd = sqrt((dplyr::pull(mcmc.list[[mcl]][, "phi"])^2) + 
                                                                  (dplyr::pull(mcmc.list[[mcl]][, "sigma"])^2)))
                                                    
                                         )}) # End dr.raw
                      
                      dr.raw <- purrr::set_names(x = dr.raw, nm = sp.names)
                      
                      doseresp.values <- purrr::map_depth(
                        .x = dr.raw,
                        .depth = 2,
                        .f = ~ {
                          lapply(X = 1:nrow(mcmc.list[[mcl]]), FUN = function(a) {
                            nth_element(vector = .x, 
                                        starting.position = a, 
                                        n = nrow(mcmc.list[[mcl]]))
                            
                          }) %>% do.call(rbind, .)}) %>%
                        purrr::set_names(x = ., nm = sp.names) 
                      
                      # Median
                      p.median <- purrr::map_depth(.x = doseresp.values, 
                                                   .depth = 2,
                                                   .f = ~ apply(X = .x, MARGIN = 2, FUN = median))
                      
                      # Lower quantiles
                      q.low <- purrr::map_depth(
                        .x = doseresp.values, 
                        .depth = 2,
                        .f = ~ {
                          tmp <- .x
                          purrr::map(.x = quants, .f = ~ apply(
                            X = tmp, MARGIN = 2,
                            FUN = quantile, (50 - .x / 2) / 100
                          ))}) 
                      
                      # Upper quantiles
                      q.up <- purrr::map_depth(
                        .x = doseresp.values,
                        .depth = 2,
                        .f = ~ {
                          tmp <- .x
                          purrr::map(.x = quants, .f = ~ apply(
                            X = tmp, MARGIN = 2,
                            FUN = quantile, (50 + .x / 2) / 100
                          ))}) 
                      
                      
                      # Return all results in a list
                      list(median = p.median, 
                           lower = q.low, 
                           upper = q.up)
                      
                    })
  
  output <- purrr::map(.x = res, 
                       .f = ~{
                         
                         out <- list()
                         out$median <- tibble::enframe(.x$median) %>% 
                           tidyr::unnest(cols = c(value)) %>% 
                           dplyr::rename(species = name) %>% 
                           dplyr::mutate(param = "median")
                         
                         out$lower <- tibble::enframe(.x$lower) %>% 
                           tidyr::unnest(cols = c(value)) %>% 
                           dplyr::rename(species = name) %>% 
                           dplyr::mutate(param = "lower")
                         
                         out$upper <- tibble::enframe(.x$upper) %>% 
                           tidyr::unnest(cols = c(value)) %>% 
                           dplyr::rename(species = name) %>% 
                           dplyr::mutate(param = "upper")
                         
                         out$posterior <- dplyr::bind_rows(out$median, out$lower, out$upper)
                         out$lower <- out$upper <- out$median <- NULL
                         out
                       })
  
  if(by.model) output <- 
    purrr::set_names(x = output, 
                     nm =  sapply(X = unique(mcmc.trace$model_ID), 
                                  FUN = function(m) rj.object$mlist[rj.object$mlist$ID == m, ]$model)) else 
                                    output <- purrr::set_names(x = output, nm = "output")
  
  output$dose.range <- dose.range
  output$quants <- quants
  
  output <- list(dat = output, 
                 sim = rj.object$dat$param$sim,
                 sp.groups = rj.object$dat$species$groups,
                 abbrev = rj.object$abbrev,
                 by.model = by.model, 
                 mlist = rj.object$mlist, 
                 p.med = rj.object$p.med,
                 covariate = covariate,
                 covariate.values = covariate.values,
                 species = species)
  
  class(output) <- c("dose_response", class(output))
  return(output)
  
}