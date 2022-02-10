#' Dose-response curves
#'
#' Compute dose-response functions from a fitted rjMCMC model.
#'
#' @export
#' @param rj.object Input rjMCMC object of class \code{rjtrace}, as returned by \code{\link{trace_rjMCMC}}.
#' @param phase Dose-response functional form: monophasic (1) or biphasic (2).
#' @param by.model Logical. If \code{TRUE}, the function subsets posterior parameter estimates to produce separate dose-response curves for each candidate model.
#' @param model.rank Rank of the model to generate curves for when \code{by.model = TRUE}.
#' @param covariate Covariate name. This argument can be used to generate dose-response curves for specific contextual covariates, conditioned on the species (group) given by \code{species}.
#' @param covariate.values A vector of values for which dose-response curves are required. Only valid for continuous covariates.
#' @param species Species name. 
#' @param credible.intervals Credible intervals. Must be a integer vector in \code{(0, 100]}. Defaults to 5-95% in 5% increments.
#' @param npts Number of quadrature points to use to integrate out the random effects when computing dose-response curves for biphasic models. Defaults to 20. 
#' 
#' 
#' @return A list object of class \code{dose_response}.
#'  
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}} \code{\link{plot.dose_response}}
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
# 
#' # Burn and thin
#' rj.trace <- trace_rjMCMC(rj.dat = rj)
#' 
#' # Get dose-response functions
#' doseR <- compile_rjMCMC(rj.trace)
#' }
#' @keywords brs dose-response rjmcmc 

compile_rjMCMC <- function(rj.object, 
                           phase = 1,
                           by.model = FALSE,
                           model.rank = 1,
                           covariate = NULL,
                           covariate.values = NULL,
                           species = NULL,
                           credible.intervals = c(95, 50, 5),
                           # seq(95, 5, by = -5)
                           npts = 20){
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if (!"rjtrace" %in% class(rj.object)) stop("Input must be an object of class <rj_trace>.")
  if (!is.null(covariate)) {
    if (!covariate %in% rj.object$dat$covariates$names) stop("Unrecognised covariate.")
  }
  if (length(covariate) > 1) stop("<covariate> must be of length 1.")
  if (!all(species %in% rj.object$dat$species$names)) stop("Unrecognised species.")
  if (is.null(species) & !is.null(covariate)) stop("<species> cannot be set to NULL when a covariate is specified.")
  if (length(species) > 1 & !is.null(covariate)) stop("Max of 1 species (group) allowed.")
  if (any(credible.intervals > 100) | any(credible.intervals <= 0)) stop("Credible intervals must lie in the interval (0, 100].")
  if (!50 %in% credible.intervals) { # Must contain the median
    credible.intervals <- sort(c(credible.intervals, 50), decreasing = TRUE)
  }
  
  if (!phase %in% c(1, 2)) stop("Unrecognised input for argument <phase>.")
  if (!rj.object$config$function.select) {
    if (rj.object$config$biphasic & phase == 1) stop("Cannot extract monophasic curves. Please set <phase> to 2.")
    if (!rj.object$config$biphasic & phase == 2) stop("Cannot extract biphasic curves. Please set <phase> to 1.")
  }
  
  if (rj.object$dat$covariates$n == 0) covariate <- NULL
  
  if (!is.null(covariate)) {
    if (is.null(covariate.values) &
        rj.object$dat$covariates$fL[[covariate]]$nL == 0) {
      stop("Must provide covariate values")
    }
    
    if (!is.null(covariate.values) &
        rj.object$dat$covariates$fL[[covariate]]$nL > 0) {
      stop("<covariate.values> ignored for factor covariates.")
    }
    
    if (!is.null(covariate.values)) {
      if (any(covariate.values < range(rj.object$dat$covariates$df[, covariate])[1]) |
          any(covariate.values > range(rj.object$dat$covariates$df[, covariate])[2])) {
        stop("Covariate values out of range")
      }
    }
  }
  
  #' ---------------------------------------------
  # Define dose range and credible intervals
  #' ---------------------------------------------
  
  if(phase == 1) pn <- "mu" else pn <- "nu"
  dose.range <- seq(rj.object$config$priors[pn, 1], rj.object$config$priors[pn, 2], length = 100)
  credible.intervals  <- sort(credible.intervals, decreasing = TRUE)
  
  #' ---------------------------------------------
  # Extract trace
  #' ---------------------------------------------
  
  mcmc.trace <- do.call(rbind, rj.object$trace) %>%
    tibble::as_tibble()
  
  if(!is.null(covariate.values)){
    
    for(tt in covariate.values){
      mcmc.trace <- mcmc.trace %>% 
        dplyr::mutate(!!paste0(covariate, "_", tt) := .data[[covariate]] *  tt)
    }
    mcmc.trace <- mcmc.trace %>% dplyr::select(-which(names(mcmc.trace) == covariate))
  }
  
  #' ---------------------------------------------
  # Split trace by model ID / phase
  #' ---------------------------------------------
  
  mcmc.trace <- dplyr::filter(mcmc.trace, phase == phase)
  m.filtered <- rj.object$ranks[rj.object$ranks$rank == model.rank, ]$model
  m.ID <- rj.object$mlist[rj.object$mlist$model == m.filtered, ]$ID
  m.groups <- rj.object$mlist[rj.object$mlist$model == m.filtered, ]$group[[1]]
  if(by.model) mcmc.trace <- dplyr::filter(mcmc.trace, model_ID == m.ID)
  
  #' ---------------------------------------------
  # Compute dose-response curves
  #' ---------------------------------------------
  
  psi.index <- which(startsWith(x = colnames(mcmc.trace), prefix = "psi"))
  
  # Identify columns corresponding to species means
  if (!is.null(species)) {
    col.index <- 
      which(colnames(mcmc.trace) %in% paste0(pn, ".", ifelse(phase == 1, "", "lower."), 
                                             which(rj.object$dat$species$names == species)))
    sp.names <- species
  } else {
    col.index <- which(startsWith(x = colnames(mcmc.trace), prefix = pn))
    sp.names <- unlist(rj.object$dat$species$names)
  }
  
  if(by.model){
    if(is.null(covariate)){
    Mg <- sp.names[sapply(X = 1:length(unique(m.groups)), FUN = function(x) which(m.groups == x)[1])]
    } else {
    Mg <- species
    }
  } else {
    Mg <- sp.names[1:length(sp.names)]
  }

    # MONOPHASIC ----------------
    
    if(phase == 1){
      
      cat("Calculating ...\n")
      
      dr.raw <- lapply(X = col.index, FUN = function(col.x) {
        
        if(!is.null(covariate)){
          
          cov.index <- which(startsWith(x = colnames(mcmc.trace), prefix = covariate))
          
          if(rj.object$dat$covariates$fL[[covariate]]$nL > 0) 
            cov.values <- list(rep(0, nrow(mcmc.trace))) else cov.values <- list()
            
            cov.values <- append(cov.values, 
                       lapply(X = cov.index, FUN = function(cc) dplyr::pull(mcmc.trace[, cc])))
            
            dr.mean <- purrr::map(.x = cov.values,
                                  .f = ~dplyr::pull(mcmc.trace[, col.x]) + .x) %>% 
              purrr::set_names(x = ., nm = rj.object$dat$covariates$fL[[covariate]]$Lnames)
            
            
        } else {
          
          dr.mean <- list(dplyr::pull(mcmc.trace[, col.x]))
        }
        
        purrr::map(.x = dr.mean,
                   .f = ~truncnorm::ptruncnorm(
                     q = rep(dose.range, each = nrow(mcmc.trace)),
                     a = rj.object$config$priors["mu", 1],
                     b = rj.object$config$priors["mu", 2],
                     mean = .x,
                     sd = sqrt((dplyr::pull(mcmc.trace[, "phi"])^2) + (dplyr::pull(mcmc.trace[, "sigma"])^2)))
                   
        )}) # End dr.raw
      
      dr.raw <- purrr::set_names(x = dr.raw, nm = sp.names)
      
      doseresp.values <- purrr::map_depth(
        .x = dr.raw,
        .depth = 2,
        .f = ~ {
          lapply(X = 1:nrow(mcmc.trace), FUN = function(a) {
            nth_element(vector = .x, starting.position = a, n = nrow(mcmc.trace))
          }) %>% do.call(rbind, .)}) %>%
        purrr::set_names(x = ., nm = sp.names) 
      
    } else if (phase == 2){
      
      # BIPHASIC ----------------

      doseresp.values <- lapply(X = Mg, FUN = function(sp.x) {
        
        dr.raw <- dplyr::select(mcmc.trace, 
                                c("omega", "psi", "tau.lower", "tau.upper", 
                                  paste0("nu.lower.", sp.x),
                                  paste0("nu.upper.", sp.x),
                                  paste0("alpha.", sp.x))) %>% 
          as.matrix()
        
        if(!is.null(covariate)){
          
          cov.index <- which(startsWith(x = colnames(mcmc.trace), prefix = covariate))
          
          if(rj.object$dat$covariates$fL[[covariate]]$nL > 0) 
            cov.values <- list(rep(0, nrow(mcmc.trace))) else cov.values <- list()
            
            cov.values <-
              append(
                cov.values,
                lapply(
                  X = cov.index,
                  FUN = function(cc) {
                    qnorm(pnorm(
                      q = dplyr::pull(mcmc.trace[, cc]),
                      mean = rj.object$config$priors[covariate, 1],
                      sd = rj.object$config$priors[covariate, 2]
                    ))
                  }
                )
              )
            
            dr.mean <- 
              purrr::map(.x = cov.values, .f = ~dplyr::pull(mcmc.trace[, psi.index]) + .x) %>% 
              purrr::set_names(x = ., nm = rj.object$dat$covariates$fL[[covariate]]$Lnames)
            
            dr.raw <- purrr::map(.x = dr.mean, .f = ~{
              dr.raw[, "psi"] <- .x
              dr.raw
            })
            
        } else {
          
          dr.raw <- list(dr.raw)
          
        }

        dr.out <- purrr::map(.x = dr.raw, .f = ~{
        
        # Set up matrices
        p.response <- matrix(data = 0, nrow = nrow(mcmc.trace), ncol = length(dose.range))
        
        for (i in 1:nrow(.x)) {
          
          cat("\r", ifelse(by.model, "Species group ", "Species "),
              sp.x, " (", sprintf("%03d", round(100*i/nrow(.x), 0)), "%)", sep = "")
          cat("\n")

          # Integrate out the random effect
          pi.individ <- pnorm(q = qnorm(p = seq(0, 1, length = (npts + 2))[-c(1, (npts + 2))],
                                        mean = .x[i, "psi"], sd = .x[i, "omega"]))
          
          # Dose-response curve corresponding to lower mixture component
          p.response.individ.lower <- 
            truncnorm::ptruncnorm(q = dose.range, 
                                  a = rj.object$dat$param$dose.range[1], 
                                  b = .x[i, which(grepl("alpha", colnames(.x)))],
                                  mean = .x[i, which(grepl("nu.lower", colnames(.x)))], 
                                  sd = .x[i, which(grepl("tau.lower", colnames(.x)))])
          
          # Dose-response curve corresponding to higher mixture component
          p.response.individ.upper <- 
            truncnorm::ptruncnorm(q = dose.range, 
                                  a = .x[i, which(grepl("alpha", colnames(.x)))], 
                                  b = rj.object$dat$param$dose.range[2], 
                                  mean = .x[i, which(grepl("nu.upper", colnames(.x)))], 
                                  sd = .x[i, which(grepl("tau.upper", colnames(.x)))])
          
          p.response.individ.lower.mat <-
            matrix(data = p.response.individ.lower,
                   nrow = npts, ncol = length(dose.range), byrow = TRUE)
          p.response.individ.lower.mat <-
            sweep(p.response.individ.lower.mat, MARGIN = 1, pi.individ, `*`)
          
          p.response.individ.upper.mat <-
            matrix(data = p.response.individ.upper,
                   nrow = npts, ncol = length(dose.range), byrow = TRUE)
          p.response.individ.upper.mat <-
            sweep(p.response.individ.upper.mat, MARGIN = 1, STATS = 1 - pi.individ, FUN = `*`)
          
          p.response.individ <- p.response.individ.lower.mat + p.response.individ.upper.mat
          
          p.response[i, ] <- Rfast::colmeans(p.response.individ)
          
        } # End for i loop
        p.response
        })
      }) %>% 
        purrr::set_names(x = ., nm = Mg)

    }
    
    # Median
    p.median <-
      purrr::map_depth(
        .depth = 2,
        .x = doseresp.values, 
        .f = ~ Rfast::colMedians(x = .x, na.rm = TRUE)
      )
    
    if(by.model & is.null(covariate)) p.median <- p.median[m.groups] %>% purrr::set_names(., sp.names)
    
    # Lower quantiles
    q.low <- purrr::map_depth(
      .x = doseresp.values, 
      .depth = 2,
      .f = ~ {
        tmp <- .x
        purrr::map(.x = credible.intervals, .f = ~ apply(
          X = tmp, MARGIN = 2,
          FUN = quantile, (50 - .x / 2) / 100, na.rm = TRUE
        ))}) 
    
    if(by.model & is.null(covariate)) q.low <- q.low[m.groups] %>% purrr::set_names(., sp.names)
    
    # Upper quantiles
    q.up <- purrr::map_depth(
      .x = doseresp.values,
      .depth = 2,
      .f = ~ {
        tmp <- .x
        purrr::map(.x = credible.intervals, .f = ~ apply(
          X = tmp, MARGIN = 2,
          FUN = quantile, (50 + .x / 2) / 100, na.rm = TRUE
        ))}) 
    
    if(by.model & is.null(covariate)) q.up <- q.up[m.groups] %>% purrr::set_names(., sp.names)
    
    # Return all results in a list
    res <- list(median = p.median, 
                lower = q.low, 
                upper = q.up)

  
    output <- purrr::map(
      .x = list(res),
      .f = ~ {
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
      }
    )
  
  if(by.model) output <- purrr::set_names(x = output, nm = m.filtered)
  
  output$dose.range <- dose.range
  output$cred.int <- credible.intervals
  
  output <- list(dat = output, 
                 phase = phase,
                 sim = rj.object$dat$param$sim,
                 sp.groups = rj.object$dat$species$groups,
                 abbrev = rj.object$abbrev,
                 by.model = by.model, 
                 mlist = rj.object$mlist, 
                 ranks = rj.object$ranks,
                 model.rank = model.rank,
                 p.med = rj.object$p.med,
                 p.med.bymodel = rj.object$p.med.bymodel,
                 covariate = covariate,
                 covariate.values = covariate.values,
                 species = species,
                 names = rj.object$dat$species$names)
  
  if(!is.null(covariate)) output$fL <- rj.object$dat$covariates$fL[[covariate]] else output$fL <- rj.object$dat$covariates$fL
  
  cat("Done!")
  class(output) <- c("dose_response", class(output))
  return(output)
  
}