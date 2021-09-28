#' Update sampler
#'
#' Update the Markov chain(s) associated with an existing rjMCMC sampler.
#'
#' @param rjdat Input rjMCMC sampler. Must be an object of class \code{rjmcmc}, as returned by \code{\link{run_rjMCMC}}.
#' @param n.iter Number of posterior samples.
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
#' rj.update <- update_rjMCMC(rjdat = rj, n.iter = 500)
#' }
#' @keywords brs dose-response rjmcmc 

update_rjMCMC <- function(rjdat, n.iter = 1000){
  
  #' -----------------------------------------------
  #' Perform function checks
  #' -----------------------------------------------
  if(!"rjmcmc" %in% class(rjdat)) stop("Input must be of class <rjmcmc>.")
  
  #' -----------------------------------------------
  #' Extract values of the parameters that require updating
  #' -----------------------------------------------
  init.pars <- purrr::map(.x = rjdat[1], .f = ~.x[["mcmc"]][c("n.burn", "n.iter", "tot.iter")]) %>% 
    purrr::flatten()
  init.ar <- purrr::map(.x = rjdat[1], .f = ~.x[["accept"]]) %>% purrr::flatten()
  init.times <- purrr::map(.x = rjdat[1], .f = ~.x[["run_time"]])
  
  #' -----------------------------------------------
  #' Run the rjMCMC sampler, using the last values of the existing object 
  #' as starting values for the updated object
  #' -----------------------------------------------
  update.res <- run_rjMCMC(dat = rjdat,
                           n.chains = length(rjdat),
                           n.burn = 0,
                           n.iter = n.iter,
                           do.update = TRUE)
  
  #' -----------------------------------------------
  #' Make necessary updates
  #' -----------------------------------------------
  for(nc in seq_along(update.res)){
    update.res[[nc]]$mcmc$n.burn <- init.pars$tot.iter
    update.res[[nc]]$mcmc$n.iter <- n.iter
    update.res[[nc]]$mcmc$tot.iter <- init.pars$tot.iter + n.iter
    update.res[[nc]]$mcmc$iter.rge <- paste0(update.res[[nc]]$mcmc$n.burn + 1, ":", 
                                             update.res[[nc]]$mcmc$tot.iter)
    update.res[[nc]]$run_time <- hms::as_hms(init.times[[1]] + update.res[[nc]]$run_time)
    
  }
  
  return(update.res)
}