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
#'                                  function.select = FALSE,
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
#' }
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
    if(burn > rj.dat[[1]]$mcmc$n.iter) stop("Burn-in cannot exceed number of iterations.")
  } else {
    if(is.null(burn)) burn <- mcmc.params$n.burn else mcmc.params$n.burn <- burn
    mcmc.params$n.iter <- (mcmc.params$tot.iter - burn) / thin
    mcmc.params$iter.rge <- paste0(mcmc.params$n.burn + thin, ":", mcmc.params$tot.iter)
  }

  mcmc.params$thin <- thin
  
  # Define model parameters
  if(rj.dat[[1]]$config$function.select){
    params <- c("mu", "phi", "sigma", "alpha", "nu", "tau", "omega", "psi")
  } else {
  if(rj.dat[[1]]$config$biphasic){
    params <- c("alpha", "nu", "tau", "omega", "psi")
  } else {
    params <- c("mu", "phi", "sigma")}
  }
  
  if(rj.dat[[1]]$dat$covariates$n > 0) params <- c(params, covariate.names)
  params <- c(params, "include.covariates")
  params <- c(params,  "model_ID", "model_size")
  params <- c(params,  "phase")
  if(rj.dat[[1]]$config$model.select) params <- c(params,  "accept.model")
  if(rj.dat[[1]]$config$covariate.select) params <- c(params,  "accept.covariates")
  if(rj.dat[[1]]$config$function.select) params <- c(params,  "accept.phase")
  
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
    dplyr::rowwise() %>% 
    dplyr::mutate(N = dplyr::n_distinct(unlist(group))) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(N) %>% 
    tibble::rowid_to_column("ID")
  
  model.l <- as.list(model.list$ID) %>% purrr::set_names(x = ., nm = model.list$model)
  model.g <- as.list(model.list$N) %>% purrr::set_names(x = ., nm = model.list$model)
  
  for(j in 1:mcmc.params$n.chains){
    rj.dat[[j]]$model_ID <- unname(unlist(model.l[rj.dat[[j]]$model])) # Model IDs
    rj.dat[[j]]$model_size <- unname(unlist(model.g[rj.dat[[j]]$model])) # Number of groupings
    if(rj.dat[[j]]$dat$covariates$n > 0){
      colnames(rj.dat[[j]]$include.covariates) <- paste0("incl.", colnames(rj.dat[[j]]$include.covariates))
    }
    if(length(rj.dat[[j]]$model_ID) == 1) 
      rj.dat[[j]]$model_ID <- rep(rj.dat[[j]]$model_ID, mcmc.params$n.iter)
    if(length(rj.dat[[j]]$model_size) == 1) 
      rj.dat[[j]]$model_size <- rep(rj.dat[[j]]$model_size, mcmc.params$n.iter)
  }
  
  #' ---------------------------------------------
  # Extract parameters and burn / thin
  #' ---------------------------------------------

  mcmc.trace <- purrr::map(.x = seq_len(length(rj.dat)),
                           .f = ~{
                             tmp <- rj.dat[[.x]][params]
                             for(ii in seq_along(tmp)){
                               if(!is.null(dim(tmp[[ii]]))) {   
                                 if(length(dim(tmp[[ii]])) > 2){
                                 tmp[[ii]] <- matrix(tmp[[ii]], 
                                                     nrow = dim(tmp[[ii]])[1], 
                                                     ncol = prod(dim(tmp[[ii]])[2], dim(tmp[[ii]])[3]),
                                                     dimnames = list(NULL,
                                as.vector(t(outer(dimnames(tmp[[ii]])[[3]], 
                                                dimnames(tmp[[ii]])[[2]], FUN = "paste0")))))
                                 }
                               } else {
                                 tmp[[ii]] <- matrix(tmp[[ii]], ncol = 1, 
                                                     dimnames = list(NULL, names(tmp)[ii]))
                               } 
                             }
                             do.call(cbind, tmp)
                           })
  
  if(rj.dat[[1]]$config$biphasic){
    nu.indices <- grep(pattern = "nu.", x = colnames(mcmc.trace[[1]]))
    tau.indices <- grep(pattern = "tau.", x = colnames(mcmc.trace[[1]]))
  } else {
    mu.indices <- grep(pattern = "mu.", x = colnames(mcmc.trace[[1]]))
  }
  
  
  # mcmc.trace <- purrr::map(.x = seq_len(length(mcmc.trace)),
  #            .f = ~ {
  #              if(rj.dat[[1]]$config$model.select) 
  #                mcmc.trace[[.x]][, "accept.model"] <- cumsum(mcmc.trace[[.x]][, "accept.model"])
  #              
  #              if(rj.dat[[1]]$config$covariate.select) 
  #                mcmc.trace[[.x]][, "accept.covariates"] <- cumsum(mcmc.trace[[.x]][, "accept.covariates"])
  #              
  #              if(rj.dat[[1]]$config$function.select) 
  #                mcmc.trace[[.x]][, "accept.phase"] <- cumsum(mcmc.trace[[.x]][, "accept.phase"])
  #              
  #              mcmc.trace[[.x]]
  #            })
  
  # Discard burn-in and thin chain(s) if desired
  for(nc in seq_len(length(mcmc.trace))){
    if (is.data.frame(mcmc.trace[[nc]]) || is.matrix(mcmc.trace[[nc]])) {
      if (burn > 0) mcmc.trace[[nc]] <- mcmc.trace[[nc]][-(1:burn), ] # Remove burn-in
      if(is.vector(mcmc.trace[[nc]])) mcmc.trace[[nc]] <- matrix(data = mcmc.trace[[nc]], ncol = 1)
      # Thin if needed
      mcmc.trace[[nc]] <- mcmc.trace[[nc]][seq(from = 1, to = nrow(mcmc.trace[[nc]]), thin), ]
    }
  }
  
  #' ---------------------------------------------
  # Model ranks
  #' ---------------------------------------------

  model.ranks <- purrr::map(.x = rj.dat, .f = "model")
  for(nc in seq_along(model.ranks)){
      if (burn > 0) model.ranks[[nc]] <- model.ranks[[nc]][-(1:burn)] # Remove burn-in
      model.ranks[[nc]] <- model.ranks[[nc]][seq(from = 1, to = length(model.ranks[[nc]]), thin)]
    }
  
  model.ranks <- do.call(c, model.ranks) %>% 
    table(.) %>% 
    tibble::enframe(.) %>% 
    dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(value))) %>% 
    dplyr::select(-value) %>% 
    dplyr::arrange(rank) %>% 
    dplyr::rename(model = name)
  
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
  
  if(rj.dat[[1]]$config$model.select){
  posterior.medians.by.model <- do.call(rbind, mcmc.trace) %>% 
    tibble::as_tibble(., .name_repair = "minimal") %>% 
    split(x = ., f = factor(.$model_ID)) %>% 
    purrr::map(.x = ., .f = ~{
      .x %>% 
        dplyr::summarise(dplyr::across(where(is.numeric), ~median(.x))) %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column() %>% 
        dplyr::rename(param = rowname, pmed = V1) 
    })
  names(posterior.medians.by.model) <- sapply(X = as.numeric(names(posterior.medians.by.model)), FUN = function(b) model.list[model.list$ID == b, ]$model)
  } else {
    posterior.medians.by.model <- list(posterior.medians) %>% 
      purrr::set_names(x = ., nm = model.list$model)
  }
  
  # Remove baseline levels of factor covariates
  if(rj.dat[[1]]$dat$covariates$n > 0){
    for(nc in seq_len(length(mcmc.trace))){
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
  
  #' ---------------------------------------------
  # Acceptance rates
  #' ---------------------------------------------
  
  AR <- 1 - coda::rejectionRate(mcmc.trace)
      
  #' ---------------------------------------------
  # Effective sample sizes
  #' ---------------------------------------------
  
  ess.params <- colnames(mcmc.trace[[1]])
  if(!rj.dat[[1]]$config$model.select) ess.params <- ess.params[!ess.params %in% c("model_ID", "model_size")]
  if(!rj.dat[[1]]$config$function.select) ess.params <- ess.params[!ess.params == "phase"]
  sink(tempfile())
  ESS <- mcmcse::ess(mcmc.trace[, ess.params]) %>% 
    tibble::enframe(.) %>% 
    dplyr::rename(parameter = name, ESS = value)
  sink()
  
  if(!rj.dat[[1]]$config$model.select) 
    ESS <- ESS %>% dplyr::filter(!parameter %in% c("model_ID", "model_size"))
  
  if(!rj.dat[[1]]$config$covariate.select) 
    ESS <- ESS %>% dplyr::filter(!stringr::str_detect(parameter, "incl."))
  
  if(!rj.dat[[1]]$config$function.select) 
    ESS <- ESS %>% dplyr::filter(!parameter == "phase")
  
  # Return output as a named list
  res <- list(dat = rj.dat[[1]]$dat,
              config = rj.dat[[1]]$config,
              trace = mcmc.trace, 
              mcmc = mcmc.params, 
              accept = AR, 
              ess = ESS,
              ranks = model.ranks,
              mlist = model.list,
              p.med = posterior.medians,
              p.med.bymodel = posterior.medians.by.model,
              abbrev = rj.dat[[1]]$abbrev)
  
  class(res) <- c("rjtrace", class(res))
  return(res) 
}