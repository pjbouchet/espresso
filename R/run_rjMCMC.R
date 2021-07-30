#' Reversible jump MCMC
#'
#' Run the rjMCMC algorithm on dose-response data.
#'
#' @export
#' @param dat A configured rjMCMC object of class \code{rjconfig}, as returned by \code{\link{configure_rjMCMC}}.
#' @param n.chains Number of MCMC chains.
#' @param n.burn Number of MCMC iterations to use as burn-in.
#' @param n.iter Number of posterior samples.
#' @param do.update Logical. If \code{TRUE}, updates an existing rjMCMC object.
#' @importFrom foreach `%dopar%`
#' @note  A progress bar is used to monitor code execution on Mac and Linux operating systems. This feature does not currently work on Windows.
#' @return A list object of class \code{rjmcmc}.
#' @author Phil J. Bouchet
#' @seealso \code{\link{configure_RJMCMC}} \code{\link{plot.rjtrace}} \code{\link{update_rjMCMC}}
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
#' }
#' @keywords brs dose-response rjmcmc 

run_rjMCMC <- function(dat, 
                       n.chains = 3,
                       n.burn = 1000, 
                       n.iter = 1000,
                       do.update = FALSE) {
  
  # If this is an update, grab the last values of the chain(s) to use as starting values
  if(do.update){
    last.iter <- purrr::map(.x = seq_along(dat), .f = ~glance(dat = dat, which.chain = .x, f = "update"))
  } else {
    last.iter <- vector(mode = "list", length = n.chains)
  }
  
  # Parallel computing
  if(n.chains > parallel::detectCores(all.tests = FALSE, logical = TRUE) -1) stop("Insufficient cores.")
  
  # Create a parallel cluster and register the backend
  cl <- suppressMessages(parallel::makeCluster(n.chains, outfile = ""))
  doParallel::registerDoParallel(cl)
  
  # Setup the rjMCMC
  rj.list <- purrr::map(.x = seq_len(n.chains), 
                        .f = ~setup_rjMCMC(rj.input = dat,
                                           n.chains = n.chains,
                                           n.burn = n.burn, 
                                           n.iter = n.iter,
                                           p.split = dat$config$move$prob[1],
                                           p.merge = dat$config$move$prob[2],
                                           m = dat$config$move$freq,
                                           move.ratio = dat$config$move$ratio,
                                           start.values = last.iter[[.x]]))
  
  # Needed to set up the progress bar on Windows
  # See https://stackoverflow.com/questions/7349570/wrapper-to-for-loops-with-progress-bar
  # sink("/dev/null")
  sink(tempfile())
  pb <- utils::txtProgressBar(min = 0, max = rj.list[[1]]$mcmc$tot.iter, style = 3)
  sink()
  
  # Launch the loop on multiple cores if needed
  rj.res <- foreach::foreach(nc = seq_len(n.chains)) %dopar% {
                               
                               rj <- rj.list[[nc]]
                               
                               for (i in 2:rj$mcmc$tot.iter) {
                                 
                                 # Print progress bar
                                 utils::setTxtProgressBar(pb, i)
                                 
                                 # Start timer
                                 if(i == 2) start.time <- Sys.time()
                                 if(i == rj$mcmc$tot.iter){
                                   end.time <- Sys.time()
                                   rj$run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, 
                                                                       time2 = start.time,
                                                                       units = "auto")), 1)}
                                 
                                 #' ---------------------------------------------------------------------
                                 # || Censoring ----
                                 #' ---------------------------------------------------------------------
                                 
                                 rj$y.ij[i, ] <- rj$y.ij[i - 1, ]
                                 
                                 if (!all(rj$dat$obs$censored == 0)) {
     
                                   rj$y.ij[i, !rj$dat$obs$censored == 0] <- 
                                     rnorm(n = sum(!rj$dat$obs$censored == 0),
                                           mean = rj$t.ij[i - 1, !rj$dat$obs$censored == 0],
                                           sd = rj$dat$obs$sd)
                                   
                                   }
                                 
                                 if(rj$config$model.select & rj$dat$species$n > 1){ 

                                   #' ---------------------------------------------------------------------
                                   # || Step 1: Split / merge ----
                                   #' ---------------------------------------------------------------------

                                   # Propose jump and new parameters
                                   proposed.jump <- propose_jump(rj.obj = rj, move.type = rj$mcmc$move$m[i])
                                   new.params <- proposal_rj(rj.obj = rj, jump = proposed.jump, iter = i)
                                   
                                   # Priors
                                   logprior.new <- uniform_prior(rj.obj = rj,
                                                                 param.name = "mu",
                                                                 param = new.params$mu)
                                   
                                   logprior.cur <- uniform_prior(rj.obj = rj,
                                                                 param.name = "mu",
                                                                 param = rj$mu[i - 1, ])
                                   
                                   # Likelihoods
                                   loglik.new <- likelihood(rj.obj = rj,
                                                            iter = i,
                                                            model = proposed.jump$model$id, 
                                                            values = list(mu = new.params$mu), 
                                                            RJ = TRUE)

                                   loglik.cur <- likelihood(rj.obj = rj,
                                                            iter = i,
                                                            model = rj$mlist[[rj$current.model]], 
                                                            RJ = TRUE)
                                   
                                   # Proposal densities
                                   logpropdens <- propdens_rj(rj.obj = rj,
                                                              param = new.params, 
                                                              jump = proposed.jump,
                                                              iter = i)

                                
                                   # Posterior ratio
                                   lognum <- loglik.new + 
                                     logprior.new + 
                                     logpropdens[1] + 
                                     proposed.jump$p[2] + 
                                     proposed.jump$jacobian
                                   
                                   logden <- loglik.cur + 
                                     logprior.cur + 
                                     logpropdens[2] + 
                                     proposed.jump$p[1]
                                   
                                   if (runif(1) < exp(lognum - logden)) {
                                     
                                     new.model <- vec_to_model(input.vector = proposed.jump$model$id,
                                                               sp.names = rj$dat$species$names)
                                     
                                     rj$model[i] <- rj$current.model <- new.model
                                     if(!new.model %in% names(rj$mlist)) rj$mlist[[new.model]] <- 
                                       proposed.jump$model$id
                                     
                                     if(!new.model %in% rj$config$clust[[1]]$model) 
                                       
                                       rj$config$clust[[1]] <- 
                                       rbind(rj$config$clust[[1]], 
                                             tibble::tibble(model = new.model, p = 0, p_scale = 0)) %>% 
                                       dplyr::mutate(p_scale = rescale_p(p))
                                     
                                     rj$mu[i, ] <- new.params$mu
                                     
                                     if(i > rj$mcmc$n.burn) rj$accept[paste0("move.", proposed.jump$type)] <- 
                                       rj$accept[[paste0("move.", proposed.jump$type)]] + 1
                                     
                                   } else {
                                     
                                     rj$model[i] <- rj$current.model
                                     rj$mu[i, ] <- rj$mu[i - 1, ]
                                   }
                                   
                                 } else {
                                   
                                   rj$mu[i, ] <- rj$mu[i - 1, ]
                                   
                                 }
                                 
                                 #' ---------------------------------------------------------------------
                                 # || Step 2: Add/remove covariates ----
                                 #' ---------------------------------------------------------------------

                                 if(rj$dat$covariates$n > 0){
                                   
                                   if(rj$config$covariate.select){
                                     
                                     # We consider all models to be equally likely when it comes to covariates.
                                     # Therefore, the priors on models cancel out.
                                     
                                     # Select a covariate at random and switch it on/off
                                     chosen.covariate <- sample(x = rj$dat$covariates$names, size = 1)
                                     proposed.covariates <- rj$include.covariates[i - 1, ]
                                     add.remove <- 1 - rj$include.covariates[i - 1, ][chosen.covariate]
                                     proposed.covariates[chosen.covariate] <- add.remove
                                     
                                     if(add.remove){
                                       
                                       # Coefficients for covariate terms
                                       beta.params <- sapply(X = chosen.covariate, 
                                                             FUN = function(w) as.matrix(rj[[w]])[i - 1, ],
                                                             simplify = FALSE, USE.NAMES = TRUE)
                                       
          # Propose new values when a covariate is added
          beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index] <- rnorm(n = rj$dat$covariates$fL[[chosen.covariate]]$nparam, mean = 0, sd = rj$config$prop$cov)
                                       
          # Prior
          logprior <- normal_prior(
            rj.obj = rj,
            param.name = chosen.covariate,
            param = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index])
                                       
          # Proposal density
          logpropdens <- sum(dnorm(x = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index], mean = 0, sd = rj$config$prop$cov, log = TRUE))
                                       
                                       # Ratio of prior / proposal density
                                       R <- logprior - logpropdens
                                       
                                     } else {
                                       
                                       beta.params <- list()
                                       beta.params[[chosen.covariate]] <- NULL
                                       
                                       # Prior
                                       logprior <- normal_prior(rj.obj = rj,
                                       param.name = chosen.covariate, 
                                       param = rj[[chosen.covariate]][i - 1, rj$dat$covariates$fL[[chosen.covariate]]$index])
                                       
                                       # Proposal density
                                       logpropdens <- sum(dnorm(x = rj[[chosen.covariate]][i - 1,                                                       rj$dat$covariates$fL[[chosen.covariate]]$index], 
                                       mean = 0, sd = rj$config$prop$cov, log = TRUE))
                                       
                                       
                                       # Ratio of prior / proposal density
                                       R <- logpropdens - logprior
                                       
                                     }
                                     
                                     # Probabilities of addition / removal
                                     pmove <- log(c(switch(as.character(add.remove), 
                                                           "1" = 1/sum(rj$include.covariates[i - 1, ] == 0), 
                                                           "0" = 1/sum(rj$include.covariates[i - 1, ] == 1)), 
                                                    1/sum(proposed.covariates == add.remove)))
                                     
                                     # The `values` argument will be the proposed values when
                                     # adding a covariate and the current values when removing
                                     # one (but these will be ignored, as the corresponding
                                     # indicator in include.cov will be 0).
                                     # include.cov corresponds to the include.covariates 
                                     # matrix in rj - it is filled with either 1 or 0 
                                     # depending on whether a covariate is included or 
                                     # omitted from the model.
                                     
                                     loglik.new <- 
                                       likelihood(rj.obj = rj,
                                                  iter = i,
                                                  model = rj$mlist[[rj$current.model]], 
                                                  include.cov = proposed.covariates,
                                           values = stats::setNames(list(beta.params[[chosen.covariate]]), 
                                                                           chosen.covariate),
                                                  RJ = TRUE)
                                     
                                     loglik.cur <- likelihood(rj.obj = rj,
                                                              iter = i,
                                                              model = rj$mlist[[rj$current.model]], 
                                                              include.cov = rj$include.covariates[i - 1, ], 
                                                              RJ = TRUE)
                                     
                                     # Posterior ratio - Model move probabilities are given in proposed.jump
                                     lognum <- loglik.new + R + pmove[2]
                                     logden <- loglik.cur + R + pmove[1]
                                     
                                     for (k in rj$dat$covariates$names) rj[[k]][i, ] <- rj[[k]][i - 1, ]
                                     
                                     if (runif(1) < exp(lognum - logden)) {
                                       
                                       rj$include.covariates[i, ] <- proposed.covariates
                                       if(add.remove) rj[[chosen.covariate]][i, ] <- beta.params[[chosen.covariate]]
                                       if(i > rj$mcmc$n.burn) rj$accept[["move.covariates"]] <- rj$accept[["move.covariates"]] + 1
                                       
                                     } else {
                                       
                                       rj$include.covariates[i, ] <- rj$include.covariates[i - 1, ]
                                       
                                     }
                                     
                                   } else {
                                     
                                     for (k in rj$dat$covariates$names) rj[[k]][i, ] <- rj[[k]][i - 1, ]
                                     
                                   }  # End if (covariate.select)
                                 } # End if (n.covariates > 0)
                                 
                                 
                                 #' ---------------------------------------------------------------------
                                 # || Step 3: Functional form ----
                                 #' ---------------------------------------------------------------------
                                 
                                 # We can generalize the standard 'monophasic' dose-response function
                                 # by allowing the threshold to come from one of two truncated normal 
                                 # distributions, one with lower exposure values than the other. This gives
                                 # a 'biphasic' function with a context-dependent part and a dose-dependent part.
                                 
                                 #' ---------------------------------------------------------------------
                                 # || Step 4: Update parameters ----
                                 #' ---------------------------------------------------------------------
                                 
                                 #'------------------------------
                                 # // t.ij ----
                                 #'------------------------------
                                 proposed.t.ij <- proposal_mh(rj.obj = rj, param.name = "t.ij", iter = i)

                                 loglik.proposed <- likelihood(rj.obj = rj,
                                                               iter = i,
                                                               model = rj$mlist[[rj$current.model]], 
                                                               values = list(t.ij = proposed.t.ij), 
                                                               lprod = FALSE)
                                 
                                 loglik.current <- likelihood(rj.obj = rj,
                                                              iter = i,
                                                              model = rj$mlist[[rj$current.model]], 
                                                              param.name = "t.ij",
                                                              lprod = FALSE)
                                 
                                 prop.forward <- propdens_mh(rj.obj = rj,
                                                             param.name = "t.ij", 
                                                             dest = proposed.t.ij, 
                                                             orig = rj$t.ij[i - 1, ])
                                 
                                 prop.backward <- propdens_mh(rj.obj = rj,
                                                              param.name = "t.ij", 
                                                              dest = rj$t.ij[i - 1, ],
                                                              orig = proposed.t.ij)
                                 
                                 accept.prob <- runif(1) < exp((loglik.proposed + prop.backward) - 
                                                                 (loglik.current + prop.forward))
                                 
                                 rj$t.ij[i, ] <- sapply(seq_len(rj$dat$trials$n), 
                                   function(x) ifelse(accept.prob[x], proposed.t.ij[x], rj$t.ij[i - 1, x]))
                                 
                                 if(i > rj$mcmc$n.burn) rj$accept[["t.ij"]] <- 
                                   rj$accept[["t.ij"]] + unname(ifelse(accept.prob, 1, 0))
                                 
                                 
                                 #'------------------------------
                                 # // sigma ----
                                 #'------------------------------
                                 proposed.sigma <- proposal_mh(rj.obj = rj, param.name = "sigma", iter = i)
                                 
                                 if(outOfbounds(rj.obj = rj, v = proposed.sigma, p = "sigma")){
                                   accept.prob <- 0
                                 } else {
                                   
                                   loglik.proposed <- likelihood(rj.obj = rj,
                                                                 iter = i,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = list(sigma = proposed.sigma))
                                   
                                   loglik.current <- likelihood(rj.obj = rj,
                                                                iter = i,
                                                                model = rj$mlist[[rj$current.model]],
                                                                param.name = "sigma")
                                   
                                   accept.prob <- exp(loglik.proposed - loglik.current)}
                                 
                                 if (runif(1) < accept.prob){
                                   rj$sigma[i] <- proposed.sigma
                                   if(i > rj$mcmc$n.burn) rj$accept[["sigma"]] <- rj$accept[["sigma"]] + 1
                                 } else {rj$sigma[i] <- rj$sigma[i - 1]}
                                 
                                 #'------------------------------
                                 # // covariates ----
                                 #'------------------------------
                                 # Only covariates currently in the model are updated
                                 # The loop below provides a generalised version of the code that
                                 # is applicable to any input covariates but impairs execution speed.
                                 # As of June 2021, only four covariates are considered in the analysis,
                                 # meaning the the loop is not hugely necessary here. It can be commented
                                 # out if more covariates are to be added in the future. 
                                 
                                 # if(rj$dat$covariates$n > 0){
                                 #   for (a in rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                 #     # Normal proposal
                                 #     proposed.value <- proposal_mh(rj.obj = rj, param.name = a) 
                                 #     if(rj$dat$covariates$fL[[a]]$nL >= 2) proposed.value <- c(0, proposed.value)
                                 #     prop.list <- list()
                                 #     prop.list[[a]] <- proposed.value
                                 #     loglik.proposed <- likelihood(rj.obj = rj,
                                 #                                   model = rj$mlist[[rj$current.model]],
                                 #                                   values = prop.list)
                                 #     loglik.current <- likelihood(rj.obj = rj,
                                 #                                  model = rj$mlist[[rj$current.model]],
                                 #                                  param.name = a)
                                 #      logprior.proposed <- normal_prior(rj.obj = rj,
                                 #                                   param.name = a,
                                 #                                   param = proposed.value)
                                 #      logprior.current <- normal_prior(rj.obj = rj,
                                 #                                       param.name = a, 
                                 #                                       param = rj[[a]][i, ])
                                 #       accept.prob <- exp((loglik.proposed + logprior.proposed) -
                                 #                      (loglik.current + logprior.current))
                                 # if (runif(1) < accept.prob){
                                 #   rj[[a]][i, ] <- proposed.value
                                 #   if(i > rj$mcmc$n.burn) rj$accept[[a]] <- rj$accept[[a]] + 1}
                                 #   }
                                 # }
                                 
                                 if(rj$dat$covariates$n > 0){
                                   
                                   #'------------------------------
                                   # // exposure history ----
                                   #'------------------------------
                                   
                                   if("exposed" %in% rj$dat$covariates$names){
                                     
                                     if("exposed" %in% rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       proposed.value <- 
                                         proposal_mh(rj.obj = rj, param.name = "exposed", iter = i)
                                       
                                       if(rj$dat$covariates$fL[["exposed"]]$nL >= 2) 
                                         proposed.value <- c(0, proposed.value)
                                       
                                       prop.list <- list()
                                       prop.list[["exposed"]] <- proposed.value
                                       
                                       loglik.proposed <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    values = prop.list)
                                       loglik.current <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    param.name = "exposed")
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "exposed",
                                                      param = proposed.value)
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "exposed",
                                                      param = rj[["exposed"]][i, ])
                                       
                                       accept.prob <- exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (runif(1) < accept.prob){
                                         rj[["exposed"]][i, ] <- proposed.value 
                                         if(i > rj$mcmc$n.burn) rj$accept[["exposed"]] <- 
                                             rj$accept[["exposed"]] + 1}
                                       
                                     }
                                   }
                                   
                                   #'------------------------------
                                   # // sonar signal ----
                                   #'------------------------------
                                   
                                   if("sonar" %in% rj$dat$covariates$names){
                                     
                                     if("sonar" %in% rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       proposed.value <- 
                                         proposal_mh(rj.obj = rj, param.name = "sonar", iter = i)
                                       
                                       if(rj$dat$covariates$fL[["sonar"]]$nL >= 2) 
                                         proposed.value <- c(0, proposed.value)
                                       
                                       prop.list <- list()
                                       prop.list[["sonar"]] <- proposed.value
                                       
                                       loglik.proposed <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    values = prop.list)
                                       
                                       loglik.current <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    param.name = "sonar")
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "sonar",
                                                      param = proposed.value)
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "sonar",
                                                      param = rj[["sonar"]][i, ])
                                       
                                       accept.prob <- exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (runif(1) < accept.prob){
                                         rj[["sonar"]][i, ] <- proposed.value 
                                         if(i > rj$mcmc$n.burn) rj$accept[["sonar"]] <- 
                                             rj$accept[["sonar"]] + 1}
                                       
                                     }
                                   }
                                   
                                   
                                   #'------------------------------
                                   # // behavioural mode ----
                                   #'------------------------------
                                   
                                   if("behaviour" %in% rj$dat$covariates$names){
                                     
                                     if("behaviour" %in% rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       proposed.value <- 
                                         proposal_mh(rj.obj = rj, param.name = "behaviour", iter = i)
                                       
                                       if(rj$dat$covariates$fL[["behaviour"]]$nL >= 2) 
                                         proposed.value <- c(0, proposed.value)
                                       
                                       prop.list <- list()
                                       prop.list[["behaviour"]] <- proposed.value
                                       
                                       loglik.proposed <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    values = prop.list)
                                       
                                       loglik.current <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    param.name = "behaviour")
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "behaviour",
                                                      param = proposed.value)
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "behaviour",
                                                      param = rj[["behaviour"]][i, ])
                                       
                                       accept.prob <- exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (runif(1) < accept.prob){
                                         rj[["behaviour"]][i, ] <- proposed.value 
                                         if(i > rj$mcmc$n.burn) rj$accept[["behaviour"]] <- 
                                             rj$accept[["behaviour"]] + 1}
                                       
                                     }
                                   }
                                   
                                   #'------------------------------
                                   # // Source-whale range ----
                                   #'------------------------------
                                   
                                   if("range" %in% rj$dat$covariates$names){
                                     
                                     if("range" %in% 
                                        rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       proposed.value <- 
                                         proposal_mh(rj.obj = rj, param.name = "range", iter = i)
                                       
                                       if(rj$dat$covariates$fL[["range"]]$nL >= 2) 
                                         proposed.value <- c(0, proposed.value)
                                       
                                       prop.list <- list()
                                       prop.list[["range"]] <- proposed.value
                                       
                                       loglik.proposed <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    values = prop.list)
                                       
                                       loglik.current <- 
                                         likelihood(rj.obj = rj,
                                                    iter = i,
                                                    model = rj$mlist[[rj$current.model]],
                                                    param.name = "range")
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "range",
                                                      param = proposed.value)
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "range",
                                                      param = rj[["range"]][i, ])
                                       
                                       accept.prob <- exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (runif(1) < accept.prob){
                                         rj[["range"]][i, ] <- proposed.value 
                                         if(i > rj$mcmc$n.burn) rj$accept[["range"]] <- 
                                             rj$accept[["range"]] + 1}
                                       
                                     }
                                   }
                                 }
                                 
                                 #'------------------------------
                                 # // mu.i ----
                                 #'------------------------------
                                 proposed.mu.i <- proposal_mh(rj.obj = rj, param.name = "mu.i", iter = i)
                                 
                                 loglik.proposed <- likelihood(rj.obj = rj,
                                                               iter = i,
                                                               model = rj$mlist[[rj$current.model]],
                                                               values = list(mu.i = proposed.mu.i),
                                                               lprod = FALSE)
                                 
                                 loglik.current <- likelihood(rj.obj = rj,
                                                              iter = i,
                                                              model = rj$mlist[[rj$current.model]],
                                                              param.name = "mu.i", 
                                                              lprod = FALSE)
                                 
                                 prop.forward <- propdens_mh(rj.obj = rj,
                                                             param.name = "mu.i", 
                                                             dest = proposed.mu.i,
                                                             orig = rj$mu.i[i - 1, ])
                                 
                                 prop.backward <- propdens_mh(rj.obj = rj,
                                                              param.name = "mu.i", 
                                                              dest = rj$mu.i[i - 1, ],
                                                              orig = proposed.mu.i)
                                 
                                 accept.prob <- runif(1) < exp((loglik.proposed + prop.backward) -
                                                                 (loglik.current + prop.forward))
                                 
                                 rj$mu.i[i, ] <- sapply(seq_len(rj$dat$whales$n), 
                                     function(x) ifelse(accept.prob[x], proposed.mu.i[x], rj$mu.i[i - 1, x]))
                                 
                                 if(i > rj$mcmc$n.burn) 
                                   rj$accept[["mu.i"]] <- 
                                   rj$accept[["mu.i"]] + unname(ifelse(accept.prob, 1, 0))
                                 
                                 #'------------------------------
                                 # // mu ----
                                 #'------------------------------
                                 proposed.mu <- proposal_mh(rj.obj = rj, param.name = "mu", iter = i)
                                 
                                 if(outOfbounds(rj.obj = rj, v = proposed.mu, p = "mu")){
                                   
                                   accept.prob <- 0
                                   
                                 } else {
                                   
                                   loglik.proposed <- likelihood(rj.obj = rj,
                                                                 iter = i,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = list(mu = proposed.mu))
                                   
                                   loglik.current <- likelihood(rj.obj = rj,
                                                                iter = i,
                                                                model = rj$mlist[[rj$current.model]],
                                                                param.name = "mu")
                                   
                                   accept.prob <- exp(loglik.proposed - loglik.current)}
                                 
                                 if (runif(1) < accept.prob){
                                   rj$mu[i, ] <- proposed.mu 
                                   if(i > rj$mcmc$n.burn) rj$accept[["mu"]] <- rj$accept[["mu"]] + 1}
                                 
                                 #'------------------------------
                                 # // phi ----
                                 #'------------------------------
                                 proposed.phi <- proposal_mh(rj.obj = rj, param.name = "phi", iter = i)
                                 
                                 if(outOfbounds(rj.obj = rj, v = proposed.phi, p = "phi")){
                                   
                                   accept.prob <- 0
                                   
                                 } else {
                                   
                                   loglik.proposed <- likelihood(rj.obj = rj,
                                                                 iter = i,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = list(phi = proposed.phi))
                                   
                                   loglik.current <- likelihood(rj.obj = rj,
                                                                iter = i,
                                                                model = rj$mlist[[rj$current.model]],
                                                                param.name = "phi")
                                   
                                   accept.prob <- exp(loglik.proposed - loglik.current)}
                                 
                                 if (runif(1) < accept.prob) {
                                   rj$phi[i] <- proposed.phi 
                                   if(i > rj$mcmc$n.burn) rj$accept[["phi"]] <- rj$accept[["phi"]] + 1
                                 } else {rj$phi[i] <- rj$phi[i - 1]}
                                 
                               } # End RJMCMC
                               
                               rj
                               
                                      
                             } # End foreach
  
  # When the function terminates
  on.exit({
    parallel::stopCluster(cl = cl) # Stop cluster
    })
  
  class(rj.res) <- c("rjmcmc", class(rj.res))
  return(rj.res)
  
  # }
  # ) # End with
  
}
