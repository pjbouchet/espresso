#' Extract trace
#'
#' Burn and thin the MCMC chain(s).
#' 
#' @export
#' @param rj.dat MCMC object of class \code{rjmcmc}, as returned by \code{\link{run_rjMCMC}}.
#' @param burn Number of iterations to discard as burn-in.
#' @param thin Thinning interval. Retains every \code{thin} iteration(s).
#' 
#' @details Adapted from the \code{fitR} package \url{https://github.com/sbfnk/fitR}.
#' 
#' @return A list object of class \code{rjtrace} containing the trace in \code{mcmc.list} format as well as additional information needed for the derivation of dose-response curves.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}}
#' @examples
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
#' @keywords brs dose-response rjmcmc 

trace_rjMCMC <- function(rj.dat, 
                         burn = NULL, 
                         thin = 1) {

  #' ---------------------------------------------
  # Function checks and initial parameterisation
  #' ---------------------------------------------
  
  if(!"rjmcmc" %in% class(rj.dat)) stop("Input must be of class <rjmcmc>.")
  
  # MCMC parameters
  mcmc.params <- rj.dat[[1]]$mcmc
  is.update <- rj.dat[[1]]$update
  covariate.names <- rj.dat[[1]]$dat$covariates$names
  covariate.levels <- lapply(X = covariate.names,
                     FUN = function(x) paste0(x, "_", rj.dat[[1]]$dat$covariates$fL[[x]]$Lnames)) %>% 
    unlist()
  
  # If this is not the first run of the sampler, update the relevant details
  if(is.update){
    if(is.null(burn)) burn <- 0 else mcmc.params$n.burn <- burn
  } else {
    if(is.null(burn)) burn <- mcmc.params$n.burn else mcmc.params$n.burn <- burn
    mcmc.params$n.iter <- mcmc.params$tot.iter - burn
    mcmc.params$iter.rge <- paste0(mcmc.params$n.burn + thin, ":", mcmc.params$tot.iter)
  }
  
  if(burn > rj.dat[[1]]$mcmc$n.iter) stop("Burn-in cannot exceed number of iterations.")
  mcmc.params$thin <- thin
  
  # Define model parameters
  params <- c("mu", "phi", "sigma", covariate.names, "model_ID", "model_size")
  if(rj.dat[[1]]$dat$covariates$n > 0) params <- c(params, "include.covariates")
  
  run_times <- purrr::map(.x = 1:length(rj.dat), .f = ~rj.dat[[.x]]$run_time) %>% 
    purrr::set_names(x = ., nm = paste0("Chain ", 1:length(rj.dat))) %>% 
    tibble::enframe() %>% 
    dplyr::rename(Chain = name, Time = value) %>% 
    dplyr::mutate(Time = do.call(c, Time))
  
  #' ---------------------------------------------
  # Compile model list and calculate model size
  #' ---------------------------------------------
  
  model.list <- purrr::map(.x = rj.dat, .f = "mlist") %>%
    purrr::map(.x = ., .f = ~tibble::enframe(.x)) %>% 
    do.call(rbind, .) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(model = name, group = value) %>% 
    dplyr::arrange(model) %>% 
    tibble::rowid_to_column("ID")
  
  for(j in 1:mcmc.params$n.chains){
    
    # Model IDs
    rj.dat[[j]]$model_ID <- sapply(X = rj.dat[[j]]$model, 
                                   FUN = function(a) model.list[model.list$model == a, ]$ID,
                                   USE.NAMES = FALSE)
    
    # Number of groupings
    rj.dat[[j]]$model_size <- sapply(X = rj.dat[[j]]$model, 
                                     FUN = function(b) n_groups(rj.dat[[j]]$mlist[[b]]),
                                     USE.NAMES = FALSE)
  }
  
  
  #' ---------------------------------------------
  # Extract parameters and burn / thin
  #' ---------------------------------------------
  mcmc.trace <- purrr::map(.x = seq_len(length(rj.dat)),
                           .f = ~do.call(cbind, rj.dat[[.x]][params]))
  
  mu.indices <- grep(pattern = "mu.", x = colnames(mcmc.trace[[1]]))
  
  # Discard burn-in and thin chain if desired
  for(nc in seq_len(length(mcmc.trace))){
    if (is.data.frame(mcmc.trace[[nc]]) || is.matrix(mcmc.trace[[nc]])) {
      if (burn > 0) mcmc.trace[[nc]] <- mcmc.trace[[nc]][-(1:burn), ] # Remove burn-in
      if(is.vector(mcmc.trace[[nc]])) mcmc.trace[[nc]] <- matrix(data = mcmc.trace[[nc]], ncol = 1)
      # Thin if needed
      mcmc.trace[[nc]] <- mcmc.trace[[nc]][seq(from = 1, to = nrow(mcmc.trace[[nc]]), thin), ]
      colnames(mcmc.trace[[nc]])[mu.indices] <- paste0("mu.", rj.dat[[1]]$dat$species$names)
    }
  }
  
  #' ---------------------------------------------
  # Covariates
  #' ---------------------------------------------
  
  # if(rj.dat[[1]]$dat$covariates$n > 0){
  #   
  #   for(nc in seq_len(length(mcmc.trace))){
  #     
  #     colnames(mcmc.trace[[nc]])[rj.dat[[1]]$dat$species$n + 2 + 1:sum(purrr::map_dbl(.x = rj.dat[[1]]$config$fL, .f = ~max(.x$index)))] <- Reduce("c", purrr::map(.x = rj.dat[[1]]$dat$covariates$dummy, .f = ~colnames(.x)))
  #     
  #     colnames(mcmc.trace[[nc]])[(ncol(mcmc.trace[[nc]]) - rj.dat[[1]]$dat$covariates$n + 1):ncol(mcmc.trace[[nc]])] <- paste0("incl.", rj.dat[[1]]$dat$covariates$names)
  #     
  #   }}
  
  
  #' ---------------------------------------------
  # Posterior medians
  #' ---------------------------------------------
  
  posterior.medians <- do.call(rbind, mcmc.trace) %>% 
    tibble::as_tibble(., .name_repair = "minimal") %>% 
    dplyr::summarise(dplyr::across(where(is.numeric), ~median(.x))) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    dplyr::rename(param = rowname, pmed = V1) 
  
  if(rj.dat[[1]]$dat$covariates$n > 0){
    
    for(nc in seq_len(length(mcmc.trace))){
      
      # Remove baseline levels of factor covariates
      baseline.lvls <- 
        colnames(mcmc.trace[[nc]])[intersect(x = purrr::map(.x = paste0(rj.dat[[1]]$dat$covariates$names, "_"),
                                                            .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[nc]])))) %>% do.call(c, .),
                                             y = which(as.integer(colSums(mcmc.trace[[nc]]) == 0)==1))]
      mcmc.trace[[nc]] <- mcmc.trace[[nc]][, !colnames(mcmc.trace[[nc]]) %in% baseline.lvls]
    } # End for loop
    
  } else {
    
    for(nc in seq_len(length(mcmc.trace))){
      mcmc.trace[[nc]] <- mcmc.trace[[nc]][, !colnames(mcmc.trace[[nc]])=="beta"]} # End for loop
    
  } # End if
  
  
  #' ---------------------------------------------
  # Create trace
  #' ---------------------------------------------
  
  mcmc.trace <- purrr::map(.x = mcmc.trace, .f = ~coda::mcmc(.x)) %>% 
    coda::mcmc.list()
  
  mcmc.params <- append(mcmc.params, list(thin = thin, run_time = run_times))
  
  # Acceptance rates
  AR <- purrr::map(.x = seq_len(mcmc.params$n.chains), .f = ~rj.dat[[.x]]["accept"])
  mcmc.params$move <- mcmc.params$move$m[burn:mcmc.params$tot.iter]
  mcmc.params$move <- table(mcmc.params$move)
  names(mcmc.params$move) <- paste0("move.", names(mcmc.params$move))
  
  AR <- purrr::map(.x = AR, .f = ~acceptance_rate(AR.obj = .x, rj.obj = rj.dat, mp = mcmc.params)) %>% 
    tibble::enframe(.) %>% 
    tidyr::unnest_wider(value) %>% 
    tidyr::unnest_wider(accept) %>% 
    dplyr::rename(chain = name) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(t.ij = mean(t.ij), mu.i = mean(mu.i)) %>% 
    dplyr::ungroup()
  
  # Return output as a named list
  res <- list(dat = rj.dat[[1]]$dat,
              config = rj.dat[[1]]$config,
              trace = mcmc.trace, 
              mcmc = mcmc.params, 
              accept = AR, 
              mlist = model.list,
              p.med = posterior.medians,
              abbrev = rj.dat[[1]]$abbrev)
  
  class(res) <- c("rjtrace", class(res))
  return(res) 
}