#' Diagnostic plots
#'
#' Generate trace, density, and autocorrelation plots from an \code{rjtrace} object.
#'
#' @export
#' @param rj.obj rjMCMC trace object of class \code{rjtrace}.
#' @param covariates.incl Logical. If \code{TRUE}, the trace is filtered to only retain posterior estimates obtained when the contextual covariates were included in the model. Only relevant when \code{covariate.select = TRUE} in \code{\link{configure_rjMCMC}}.
#' @param phase Integer. If used, will only generate plots for the parameters of the monophasic (1) or biphasic (2) model.
#' @inheritParams plot.gvs
#' @import ggplot2
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
#' 
#' # Get density and trace plots
#' plot(rj.trace)
#' }
#' @keywords brs dose-response rjmcmc 

plot.rjtrace <- function(rj.obj, 
                         param.name = NULL,
                         phase = NULL,
                         type = "both", # or trace or density
                         adjust = 2,
                         gvals = NULL,
                         priors = NULL,
                         covariates.incl = FALSE,
                         autocorr = FALSE,
                         individual = TRUE){
  
  if(rj.obj$dat$covariates$n == 0) covariates.incl <- FALSE
  if(!is.null(param.name)) { if(any(param.name %in% c("phi", "sigma", "mu"))) phase <- 1 }
  if(!is.null(param.name)) { if(any(param.name %in% c("alpha", "tau", "psi", "omega", "nu"))) phase <- 2 }

  
  #' ---------------------------------------------
  # Extract the trace and rename columns
  #' ---------------------------------------------
  mcmc.trace <- rj.obj$trace
  all.indices <- NULL
  if(rj.obj$config$function.select | rj.obj$config$biphasic){
    nu.indices <- grep(pattern = "nu.", x = colnames(mcmc.trace[[1]]))
    tau.indices <- grep(pattern = "tau.", x = colnames(mcmc.trace[[1]]))
    alpha.indices <- grep(pattern = "alpha.", x = colnames(mcmc.trace[[1]]))
    all.indices <- c(all.indices, nu.indices, alpha.indices)
  }
  if(rj.obj$config$function.select | !rj.obj$config$biphasic) {
    mu.indices <- grep(pattern = "mu.", x = colnames(mcmc.trace[[1]]))
    all.indices <-  c(all.indices, mu.indices)
  }
  
  all.indices <- sort(all.indices)
  
  cov.names <- purrr::map2(.x = names(rj.obj$dat$covariates$fL), 
                           .y = purrr::map(.x = rj.obj$dat$covariates$fL, "Lnames"),
                           .f = ~{if(!is.null(.y)) paste0(.x, "_", .y[2:length(.y)])}) %>% 
    purrr::compact() %>% 
    do.call(c, .)
  
  cov.indices <- which(colnames(mcmc.trace[[1]]) %in% cov.names)
                         
  for(nc in seq_len(length(mcmc.trace))){
    
    colnames(mcmc.trace[[nc]])[cov.indices] <-
      gsub(pattern = "_", replacement = " (", x = colnames(mcmc.trace[[nc]])[cov.indices])
    colnames(mcmc.trace[[nc]])[cov.indices] <- paste0(colnames(mcmc.trace[[nc]])[cov.indices], ")")
    
    for(nn in c("mu", "alpha", "nu.lower", "nu.upper")){
      colnames(mcmc.trace[[nc]])[grepl(pattern = nn, x = colnames(mcmc.trace[[nc]]))] <- 
        paste0(nn, " (", rj.obj$dat$species$names, ")")
    }
  }
  
  if(is.null(param.name)) mpars <- colnames(mcmc.trace[[1]]) else 
    mpars <- colnames(mcmc.trace[[1]])[unlist(purrr::map(.x = param.name, .f = ~which(startsWith(colnames(mcmc.trace[[1]]), prefix = .x))))]
  
  if(!rj.obj$config$model.select) mpars <- mpars[!mpars %in% c("model_size", "model_ID")]
  if(!rj.obj$config$covariate.select) mpars <- mpars[!mpars %in% paste0("incl.", rj.obj$dat$covariates$names)]
  if(!rj.obj$config$function.select) mpars <- mpars[!mpars %in% "phase"]
  
  if(!is.null(phase)){
    
    if(phase == 1){
      
      mpars <- mpars[!grepl(pattern = "alpha", x = mpars)]
      mpars <- mpars[!grepl(pattern = "nu", x = mpars)]
      mpars <- mpars[!grepl(pattern = "tau", x = mpars)]
      mpars <- mpars[!grepl(pattern = "omega", x = mpars)]
      mpars <- mpars[!grepl(pattern = "psi", x = mpars)]
      
    } else if (phase == 2){
      
      mpars <- mpars[!grepl(pattern = "mu", x = mpars)]
      mpars <- mpars[!grepl(pattern = "phi", x = mpars)]
      mpars <- mpars[!grepl(pattern = "sigma", x = mpars)]
    }
    
    phase.n <- purrr::map_dbl(.x = 1:length(mcmc.trace), .f = ~mcmc.trace[[.x]][mcmc.trace[[.x]][,"phase"] == phase, ] %>% nrow()) %>%
      min(.)
    
    mcmc.trace <- purrr::map(.x = 1:length(mcmc.trace), 
                             .f = ~mcmc.trace[[.x]][mcmc.trace[[.x]][, "phase"] == phase, mpars, drop = FALSE])
    
    for(ch in seq_len(rj.obj$mcmc$n.chains)){
      mcmc.trace[[ch]] <- mcmc.trace[[ch]][(nrow(mcmc.trace[[ch]]) - phase.n + 1):nrow(mcmc.trace[[ch]]), , drop = FALSE]}
    
    mcmc.trace <- purrr::map(.x = mcmc.trace, .f = ~coda::as.mcmc(.x))
    mcmc.trace <- coda::as.mcmc.list(mcmc.trace)
    
  }
  
  if(is.null(param.name)) bpars <- mpars else bpars <- purrr::map(.x = param.name, .f = ~colnames(mcmc.trace[[1]])[grepl(pattern = .x, x = colnames(mcmc.trace[[1]]))]) %>% unlist()
  
  mpars <- mpars[!grepl(pattern = "accept", x = mpars)]
  bpars <- bpars[!grepl(pattern = "accept", x = bpars)]
  
  # gvals <- rj.post$dat$param[c("mu", "sigma", "phi")] |> unlist()
  
  if(is.null(param.name)){
    
    if(covariates.incl){
      
      do.filter <- TRUE
      
    } else {
    
      do.filter <- FALSE
      MCMC_trace(mcmc.trace, 
                 type = type,
                 priors = priors,
                 gvals = gvals,
                 adjust = adjust,
                 iter = rj.obj$mcmc$n.iter, 
                 pdf = FALSE, 
                 ind = individual,
                 params = mpars)
      
    }
    
  } else {
    
    pn <- purrr::map(.x = param.name, 
                     .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[1]])))) %>%
      unlist()
    
    mcmc.trace <- mcmc.trace[, pn, drop = FALSE]
    mpars <- colnames(mcmc.trace[[1]])
    
    if(any(param.name %in% rj.obj$dat$covariates$names) & covariates.incl){
      
      do.filter <- TRUE
      
    } else {
      
      do.filter <- FALSE
      
      MCMC_trace(mcmc.trace, 
                 type = type,
                 priors = priors,
                 gvals = gvals,
                 adjust = adjust,
                 iter = rj.obj$mcmc$n.iter, 
                 pdf = FALSE, 
                 ind = individual,
                 params = mpars)
      
    }
    
  }
  
  if(do.filter){
    
    cov.cols <- purrr::map(.x = rj.obj$dat$covariates$names,
                           .f = ~which(grepl(pattern = .x, x = colnames(mcmc.trace[[1]]))))
    cov.names <- colnames(mcmc.trace[[1]])[cov.cols[[1]]]
    
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
               type = type,
               priors = priors,
               adjust = adjust,
               gvals = gvals,
               iter = rj.obj$mcmc$n.iter, 
               pdf = FALSE, 
               ind = individual,
               params = mpars[!mpars %in% cov.names])
    
    purrr::walk2(.x = cov.trace.final,
                 .y = cov.n,
                 .f = ~MCMC_trace(.x, type = type, gvals = gvals, adjust = adjust, priors = priors, iter = .y, pdf = FALSE, ind = individual))
  }
  
  if(autocorr){
    
    ac.plot <- bayesplot::mcmc_acf(x = mcmc.trace, pars = bpars[!bpars == "phase"])
    ac.plot$data$Parameter <- as.character(ac.plot$data$Parameter)
        
    # Create dummy faceted plot - to be able to calculate the required number of pages
    tmp.plot <- tryCatch(expr = {ggplot2::ggplot(data = na.omit(gdat$data), ggplot2::aes(x = Lag, y = AC)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(Chain)), alpha = 0.5) +
        ggforce::facet_wrap_paginate(ggplot2::vars(Parameter), ncol = 4, nrow = 4, page = 1)},
        error = function(cond) NULL)
 
    if(is.null(tmp.plot)) n.pages <- 1 else n.pages <- ggforce::n_pages(tmp.plot)
    
    for(pg in 1:n.pages){
    ac.gg <- ggplot2::ggplot(data = na.omit(ac.plot$data), ggplot2::aes(x = Lag, y = AC)) + 
      ggplot2::geom_point(ggplot2::aes(colour = factor(Chain)), alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(colour = factor(Chain))) +
      {if(is.null(tmp.plot)) ggforce::facet_wrap_paginate(ggplot2::vars(Parameter)) } +
      {if(!is.null(tmp.plot)) ggforce::facet_wrap_paginate(ggplot2::vars(Parameter), 
                                                           ncol = 4, nrow = 4, page = pg) } + 
      ggplot2::scale_colour_manual(values = gg_color_hue(length(mcmc.trace))) +
      ggplot2::theme(axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            plot.margin = margin(t = 0.15, r = 1, b = 0.15, l = 0.15, "cm"),
            legend.position = "top",
            legend.text = element_text(size = 11)) + 
      ggplot2::labs(colour = "Chain", y = "Autocorrelation")
    print(ac.gg)
    }
  }
  
}