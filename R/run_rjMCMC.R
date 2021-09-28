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
  
  if(do.update){
    dat <- append(dat[[1]][["dat"]], dat[[1]]["config"])
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
                                           moves = dat$config$move$moves,
                                           m = dat$config$move$freq,
                                           move.ratio = dat$config$move$ratio,
                                           do.update = do.update,
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
                                 
                                 #' -----------------------------------------------------
                                 # Step 1: Split / merge ----
                                 #' -----------------------------------------------------

                                 if(rj$config$model.select & rj$dat$species$n > 1) { 

                                     # Propose jump and new parameters
                                     proposed.jump <- propose_jump(rj.obj = rj, 
                                                                   move.type = rj$mcmc$move$m[i],
                                                                   phase = rj$phase[i - 1])
                                     
                                     new.params <- proposal_rj(rj.obj = rj,
                                                               jump = proposed.jump,
                                                               phase = rj$phase[i - 1],
                                                               huelsenbeck = FALSE)
                                     
                                     # Priors
                                     logprior.new <- 
                                       uniform_prior(rj.obj = rj,
                                       param.name = if(rj$phase[i - 1] == 1) "mu" else c("nu", "alpha"),
                                            param = new.params)
                                     
                                     logprior.cur <- 
                                       uniform_prior(rj.obj = rj,
                                       param.name = if(rj$phase[i - 1] == 1) "mu" else c("nu", "alpha"),
                                       param = if(rj$phase[i - 1] == 1) 
                                         list(out = list(mu = rj$mu[i - 1, ])) else 
                                                       list(out = list(nu = t(rj$nu[i - 1, ,]), 
                                                                       alpha = rj$alpha[i - 1, ])))
                                     
                                     # Likelihoods
                                     loglik.new <-
                                       likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                                                  param.name = if(rj$phase[i - 1] == 1) "mu" else
                                         c("nu", "alpha"),
                                                  rj.obj = rj,
                                                  model = proposed.jump$model$id,
                                                  values = new.params$out,
                                                  included.cov = NULL,
                                                  RJ = TRUE,
                                                  lprod = TRUE)
                                     
                                     loglik.cur <- 
                                       likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                                                  param.name = if(rj$phase[i - 1] == 1) "mu" else
                                         c("nu", "alpha"),
                                                  rj.obj = rj,
                                                  model = rj$mlist[[rj$current.model]], 
                                                  values = NULL,
                                                  included.cov = NULL,
                                                  RJ = TRUE, 
                                                  lprod = TRUE)
                                     
                                     # Proposal densities
                                     logpropdens <- propdens_rj(rj.obj = rj,
                                                                param = new.params,
                                                                jump = proposed.jump,
                                                                phase = rj$phase[i - 1])
                                     
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
                                       
                                       # Add model to list
                                       if(!new.model %in% names(rj$mlist)) rj$mlist[[new.model]] <- 
                                         proposed.jump$model$id
                                       
                                       if(!new.model %in% rj$config$clust[[1]]$model){ 
                                         rj$config$clust[[1]] <- 
                                         rbind(rj$config$clust[[1]], 
                                               tibble::tibble(model = new.model, p = 0, p_scale = 0)) %>% 
                                         dplyr::mutate(p_scale = rescale_p(p))
                                       }
                                       
                                       if(rj$phase[i - 1] == 1){
                                       rj$mu[i, ] <- new.params$out$mu
                                       if(rj$config$function.select){
                                         rj$nu[i, ,] <- rj$nu[i - 1, ,]
                                         rj$alpha[i, ] <- rj$alpha[i - 1, ]
                                       }}
                                       
                                       if(rj$phase[i - 1] == 2){
                                       rj$nu[i, ,] <- t(new.params$out$nu)
                                       rj$alpha[i, ] <- new.params$out$alpha
                                       if(rj$config$function.select){
                                         rj$mu[i, ] <- rj$mu[i - 1, ]
                                       }}
                                        
                                       if(i > rj$mcmc$n.burn){
                                         rj$accept[paste0(ifelse(rj$phase[i] == 1, "mono.", "bi."),
                                                          "move.", proposed.jump$type)] <- 
                                         rj$accept[[paste0(ifelse(rj$phase[i] == 1, "mono.", "bi."),
                                                           "move.", proposed.jump$type)]] + 1 }
                                       
                                     } else { # If proposal not accepted
                                       
                                       rj$model[i] <- rj$current.model <- rj$model[i - 1]
                                       if(rj$config$function.select | rj$phase[i - 1] == 1){
                                       rj$mu[i, ] <- rj$mu[i - 1, ]}
                                       if(rj$config$function.select | rj$phase[i - 1] == 2){
                                       rj$nu[i, ,] <- rj$nu[i - 1, ,]
                                       rj$alpha[i, ] <- rj$alpha[i - 1, ]
                                       }
                                        
                                     }
 
                                 } else { # If no model selection
                                   
                                   if(rj$config$function.select | rj$phase[i - 1] == 1){
                                   rj$mu[i, ] <- rj$mu[i - 1, ]}
                                   if(rj$config$function.select | rj$phase[i - 1] == 2){
                                   rj$nu[i, ,] <- rj$nu[i - 1, ,]
                                   rj$alpha[i, ] <- rj$alpha[i - 1, ]}
                                   
                                 }
                                 
                                 if(rj$config$function.select | rj$phase[i - 1] == 1){
                                 rj$iter["mu"] <- rj$iter["mu"] + 1}
                                 if(rj$config$function.select | rj$phase[i - 1] == 2){
                                 rj$iter["nu"] <- rj$iter["nu"] + 1
                                 rj$iter["alpha"] <- rj$iter["alpha"] + 1}
               
                                 #' -----------------------------------------------------
                                 # Step 2: Add/remove covariates ----
                                 #' -----------------------------------------------------

                                 if(rj$dat$covariates$n > 0){
                                   
                                   if(rj$config$covariate.select){
                                     
                                     # We consider all models to be equally likely 
                                     # when it comes to covariates.
                                     # Therefore, the priors on models cancel out.
                                     
                                     # Select a covariate at random and switch it on/off
                                     chosen.covariate <- sample(x = rj$dat$covariates$names, size = 1)
                                     proposed.covariates <- rj$include.covariates[i - 1, ]
                                     add.remove <- 1 - rj$include.covariates[i - 1, ][chosen.covariate]
                                     proposed.covariates[chosen.covariate] <- add.remove
                                     
                                     if(add.remove){
                                       
                                       # Coefficients for covariate terms
                                       beta.params <- sapply(X = chosen.covariate,
                                                             FUN = function(w)
                                                               as.matrix(rj[[w]])[i - 1, ],
                                                             simplify = FALSE, USE.NAMES = TRUE)
                                       
          # Propose new values when a covariate is added
          beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index] <- rnorm(n = rj$dat$covariates$fL[[chosen.covariate]]$nparam, mean = 0, sd = rj$config$prop$cov)
              
                                                              
                                       # Prior
                                       logprior <- normal_prior(
                                         rj.obj = rj,
                                         param.name = chosen.covariate,
                                         param = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index])
                                       
                                       # Proposal density
                                       logpropdens <- sum(dnorm(x = beta.params[[chosen.covariate]][rj$dat$covariates$fL[[chosen.covariate]]$index], mean = 0, 
                                                      sd = rj$config$prop$cov, log = TRUE))
                                       
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
                                     # indicator in included.cov will be 0).
                                     # included.cov corresponds to the include.covariates 
                                     # matrix in rj - it is filled with either 1 or 0 
                                     # depending on whether a covariate is included or 
                                     # omitted from the model.
                                     
                                     loglik.new <- 
                                       likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                                                  param = "covariates",
                                                  rj.obj = rj,
                                                  model = rj$mlist[[rj$current.model]], 
                                                  values = stats::setNames(list(beta.params[[chosen.covariate]]), 
                                                                           chosen.covariate),
                                                  included.cov = proposed.covariates,
                                                  RJ = TRUE,
                                                  lprod = TRUE)
                                     
                                     loglik.cur <- 
                                       likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                                                  param = "covariates",
                                                  rj.obj = rj,
                                                  model = rj$mlist[[rj$current.model]],
                                                  values = NULL,
                                                  included.cov = rj$include.covariates[i - 1, ], 
                                                  RJ = TRUE,
                                                  lprod = TRUE)
                                     
                                     # Posterior ratio - Model move probabilities are given in proposed.jump
                                     lognum <- loglik.new + R + pmove[2]
                                     logden <- loglik.cur + R + pmove[1]

                                     
                                     for (k in rj$dat$covariates$names) rj[[k]][i, ] <- rj[[k]][i - 1, ]
                                     
                                     if (runif(1) < exp(lognum - logden)) {
                                       
                                       rj$include.covariates[i, ] <- proposed.covariates
                                       if(add.remove) rj[[chosen.covariate]][i, ] <- 
                                           beta.params[[chosen.covariate]]
                                       if(i > rj$mcmc$n.burn) 
                                         rj$accept[["move.covariates"]] <- rj$accept[["move.covariates"]] + 1
                                       
                                     } else {
                                       
                                       rj$include.covariates[i, ] <- rj$include.covariates[i - 1, ]
                                       
                                     }
                                     
                                   } else {
                                     
                                     for (k in rj$dat$covariates$names){
                                       rj[[k]][i, ] <- rj[[k]][i - 1, ]}
                                     
                                   }  # End if (covariate.select)
                                 } # End if (n.covariates > 0)

                                 rj$iter[c("covariates", rj$dat$covariates$names)] <- 
                                   rj$iter[c("covariates", rj$dat$covariates$names)] + 1
                                 
                                 #' ---------------------------------------------------------------------
                                 # Step 3: Functional form ----
                                 #' ---------------------------------------------------------------------
                                 
                                 if(rj$config$function.select){
                                   
                                   # We can generalize the standard 'monophasic' dose-response function
                                   # by allowing the threshold to come from one of two truncated normal 
                                   # distributions, one with lower exposure values than the other. 
                                   # This gives a 'biphasic' function with a context-dependent part 
                                   # and a dose-dependent part.
                                   
                                   ff.prop <- proposal_ff(rj.obj = rj, from.phase = rj$phase[i - 1])
                                   
                                   # Likelihoods
                                   loglik.new <-
                                     likelihood(biphasic = ifelse(ff.prop$to.phase == 1, FALSE, TRUE),
                                                param.name = if(ff.prop$to.phase == 1) c("mu", "mu.i", "phi", "sigma") else c("nu", "alpha", "tau", "omega", "psi", "mu.ij", "psi.i", "k.ij"),
                                                rj.obj = rj,
                                                model = rj$mlist[[rj$current.model]],
                                                values = ff.prop,
                                                included.cov = NULL,
                                                RJ = TRUE,
                                                lprod = TRUE)
                                   
                                   loglik.current <-
                                     likelihood(biphasic = ifelse(rj$phase[i - 1] == 1, FALSE, TRUE),
                                                param.name = NULL,
                                                rj.obj = rj,
                                                model = rj$mlist[[rj$current.model]],
                                                values = NULL,
                                                included.cov = NULL,
                                                RJ = TRUE,
                                                lprod = TRUE)
                                   
                                   # Priors
                                   logprior.new <- ff_prior(rj.obj = rj, 
                                                            param = ff.prop, 
                                                            phase = ff.prop$to.phase)
                                   
                                   logprior.cur <- ff_prior(rj.obj = rj, 
                                                            param = NULL,
                                                            phase = rj$phase[i - 1])
                                   
                                   # Proposal densities
                                   logpropdens <- propdens_ff(rj.obj = rj, param = ff.prop)
                                   
                                   # Posterior ratio (Jacobian is 1 and p is 1)
                                   lognum <- loglik.new + 
                                     logprior.new + 
                                     logpropdens[1]
                                   
                                   logden <- loglik.cur + 
                                     logprior.cur + 
                                     logpropdens[2]
                                   
                                   if (runif(1) < exp(lognum - logden)) {
                                     
                                     rj$phase[i] <- ff.prop$to.phase
                                     
                                     if(i > rj$mcmc$n.burn){
                                       
                                       if(ff.prop$to.phase == 1)
                                       rj$accept["to.monophasic"] <- rj$accept["to.monophasic"] + 1
                                       
                                       if(ff.prop$to.phase == 2)
                                         rj$accept["to.biphasic"] <- rj$accept["to.biphasic"] + 1
                                     }
                                     

                                     if(ff.prop$to.phase == 1){
                                       
                                       rj$sigma[i] <- ff.prop$sigma
                                       rj$mu.i[i, ] <- ff.prop$mu.i
                                       rj$mu[i, ] <- ff.prop$mu
                                       rj$phi[i] <- ff.prop$phi

                                       rj$tau[i, ] <- rj$tau[i - 1, ]
                                       rj$omega[i] <- rj$omega[i - 1]
                                       rj$psi[i] <- rj$psi[i - 1]
                                       rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
                                       rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
                                       rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
                                       rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]
                                       
                                     } else {
                                       
                                       rj$alpha[i, ] <- ff.prop$alpha
                                       rj$nu[i, , 1] <- ff.prop$nu[1, ]
                                       rj$nu[i, , 2] <- ff.prop$nu[2, ]
                                       rj$tau[i, ] <- ff.prop$tau
                                       rj$omega[i] <- ff.prop$omega
                                       rj$psi[i] <- ff.prop$psi
                                       rj$mu.ij[i, ,] <- ff.prop$mu.ij
                                       rj$psi.i[i, ] <- ff.prop$psi.i
                                       rj$k.ij[i, ] <- ff.prop$k.ij
                                       rj$pi.ij[i, ] <- ff.prop$pi.ij
                                       
                                       rj$sigma[i] <- rj$sigma[i - 1]
                                       rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
                                       rj$phi[i] <- rj$phi[i - 1]
                                     }
                                     
                                   } else {
                                     
                                     rj$phase[i] <- rj$phase[i - 1]
                                     
                                     rj$sigma[i] <- rj$sigma[i - 1]
                                     rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
                                     rj$phi[i] <- rj$phi[i - 1]
                                     
                                     rj$tau[i, ] <- rj$tau[i - 1, ]
                                     rj$omega[i] <- rj$omega[i - 1]
                                     rj$psi[i] <- rj$psi[i - 1]
                                     rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
                                     rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
                                     rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
                                     rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]
                                     
                                   }
                                   
                                 } else {
                                   
                                   rj$phase[i] <- rj$phase[i - 1]
                                   
                                   if(rj$config$function.select | rj$phase[i - 1] == 1){
                                   rj$sigma[i] <- rj$sigma[i - 1]
                                   rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
                                   rj$phi[i] <- rj$phi[i - 1]
                                   }
                                   
                                   if(rj$config$function.select | rj$phase[i - 1] == 2){
                                   rj$tau[i, ] <- rj$tau[i - 1, ]
                                   rj$omega[i] <- rj$omega[i - 1]
                                   rj$psi[i] <- rj$psi[i - 1]
                                   rj$mu.ij[i, ,] <- rj$mu.ij[i - 1, ,]
                                   rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
                                   rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
                                   rj$pi.ij[i, ] <- rj$pi.ij[i - 1, ]

                                   }
                                   
                                 }
                                 
                                 rj$t.ij[i, ] <- rj$t.ij[i - 1, ] # t.ij stay the same
                                 rj$iter["t.ij"] <- rj$iter["t.ij"] + 1
                                 
                                 if(rj$config$function.select | rj$phase[i - 1] == 1){
                                 rj$iter[c("sigma", "mu.i", "phi")] <- 
                                   rj$iter[c("sigma", "mu.i", "phi")] + 1
                                 }
                                 
                                 if(rj$config$function.select | rj$phase[i - 1] == 2){
                                 rj$iter[c("tau", "omega", "psi", "mu.ij", "k.ij", "psi.i", "pi.ij")] <- 
                                rj$iter[c("tau", "omega", "psi", "mu.ij", "k.ij", "psi.i", "pi.ij")] + 1
                                 }
                                 
                                 #' ---------------------------------------------------------------------
                                 # Step 4: Update parameters ----
                                 #' ---------------------------------------------------------------------
                                 
                                 #'------------------------------
                                 ### covariates ----
                                 #'------------------------------
                                 
                                 if(rj$dat$covariates$n > 0){
                                   
                                   #'------------------------------
                                   # // exposure history ----
                                   #'------------------------------
                                   
                                   if("exposed" %in% rj$dat$covariates$names){
                                     
                                     if("exposed" %in% 
                                        rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       prop.list <- proposal_mh(rj.obj = rj, param.name = "exposed")

                                         loglik.proposed <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "exposed",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = prop.list,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                         
                                         loglik.current <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "exposed",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = NULL,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                      
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "exposed",
                                                      param = prop.list[["exposed"]][2])
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "exposed",
                                                      param = rj[["exposed"]][rj$iter["exposed"], 2])
                                       
                                       accept.prob <- 
                                         runif(1) < exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (accept.prob){
                                         rj[["exposed"]][i, ] <- unlist(prop.list)
                                         if(i > rj$mcmc$n.burn) rj$accept[["exposed"]] <- 
                                             rj$accept[["exposed"]] + 1
                                       } 
                                     }
                                   }
                                   
                                   #'------------------------------
                                   # // sonar signal ----
                                   #'------------------------------
                                   
                                   if("sonar" %in% rj$dat$covariates$names){
                                     
                                     if("sonar" %in% 
                                        rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       prop.list <- proposal_mh(rj.obj = rj, param.name = "sonar")

                                         loglik.proposed <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "sonar",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = prop.list,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                         
                                         loglik.current <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "sonar",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = NULL,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "sonar",
                                                      param = prop.list[["sonar"]][2])
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "sonar",
                                                      param = rj[["sonar"]][rj$iter["sonar"], 2])
                                       
                                       accept.prob <- 
                                         runif(1) < exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (accept.prob){
                                         rj[["sonar"]][i, ] <- unlist(prop.list)
                                         if(i > rj$mcmc$n.burn) rj$accept[["sonar"]] <- 
                                             rj$accept[["sonar"]] + 1
                                       } 
                                     }
                                   }
                                   
                                   #'------------------------------
                                   # // behavioural mode ----
                                   #'------------------------------
                                   
                                   if("behaviour" %in% rj$dat$covariates$names){
                                     
                                     if("behaviour" %in% 
                                        rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       prop.list <- proposal_mh(rj.obj = rj, param.name = "behaviour")
                                       
                                         loglik.proposed <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "behaviour",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = prop.list,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                         
                                         loglik.current <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "behaviour",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = NULL,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "behaviour",
                                                      param = prop.list[["behaviour"]][2])
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "behaviour",
                                                      param = rj[["behaviour"]][rj$iter["behaviour"], 2])
                                       
                                       accept.prob <- 
                                         runif(1) < exp((loglik.proposed + logprior.proposed) - 
                                                            (loglik.current + logprior.current))
                                       
                                       if (accept.prob){
                                         rj[["behaviour"]][i, ] <- unlist(prop.list)
                                         if(i > rj$mcmc$n.burn) rj$accept[["behaviour"]] <- 
                                             rj$accept[["behaviour"]] + 1
                                       } 
                                     }
                                   }
                                   
                                   #'------------------------------
                                   # // Source-whale range ----
                                   #'------------------------------
                                   
                                   if("range" %in% rj$dat$covariates$names){
                                     
                                     if("range" %in% 
                                        rj$dat$covariates$names[rj$include.covariates[i, ] == 1]) {
                                       
                                       prop.list <- proposal_mh(rj.obj = rj, param.name = "range")
                                         
                                       loglik.proposed <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "range",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = prop.list,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                         
                                         loglik.current <- 
                                           likelihood(biphasic = rj$phase[i] == 2,
                                                      param.name = "range",
                                                      rj.obj = rj,
                                                      model = rj$mlist[[rj$current.model]],
                                                      values = NULL,
                                                      included.cov = NULL,
                                                      RJ = FALSE,
                                                      lprod = TRUE)
                                       
                                       logprior.proposed <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "range",
                                                      param = prop.list[["range"]])
                                       
                                       logprior.current <- 
                                         normal_prior(rj.obj = rj,
                                                      param.name = "range",
                                                      param = rj[["range"]][rj$iter["range"]])
                                       
                                       accept.prob <- runif(1) < 
                                         exp((loglik.proposed + logprior.proposed) - 
                                             (loglik.current + logprior.current))
                                       
                                       if (accept.prob){
                                         rj[["range"]][i, ] <- unlist(prop.list)
                                         if(i > rj$mcmc$n.burn) rj$accept[["range"]] <- 
                                             rj$accept[["range"]] + 1
                                       } 
                                     }
                                   }
                                   
                                 } # End covariates.n > 0
                                 
                                 if(rj$phase[i] == 1){
                                   
                                   #'------------------------------
                                   ### t.ij (monophasic) ----
                                   #'------------------------------
                                   proposed.t.ij <- proposal_mh(rj.obj = rj, param.name = "t.ij")
                                   
                                   loglik.proposed <- likelihood(param.name = "t.ij",
                                                                 rj.obj = rj,
                                                                 model = rj$mlist[[rj$current.model]], 
                                                                 values = proposed.t.ij, 
                                                                 included.cov = NULL,
                                                                 RJ = FALSE,
                                                                 lprod = FALSE)
                                   
                                   loglik.current <- likelihood(param.name = "t.ij",
                                                                rj.obj = rj,
                                                                model = rj$mlist[[rj$current.model]], 
                                                                values = NULL,
                                                                included.cov = NULL,
                                                                RJ = FALSE,
                                                                lprod = FALSE)
                                   
                                   prop.forward <- propdens_mh(rj.obj = rj,
                                                               param.name = "t.ij", 
                                                               dest = proposed.t.ij[[1]], 
                                                               orig = rj$t.ij[i - 1, ])
                                   
                                   prop.backward <- propdens_mh(rj.obj = rj,
                                                                param.name = "t.ij", 
                                                                dest = rj$t.ij[i - 1, ],
                                                                orig = proposed.t.ij[[1]])
                                   
                                   accept.prob <- runif(rj$dat$trials$n) < 
                                     exp((loglik.proposed + prop.backward) - 
                                           (loglik.current + prop.forward))
                                   
                                   # rj$t.ij[i, ] <- rj$t.ij[i - 1, ]
                                   rj$t.ij[i, accept.prob] <- proposed.t.ij[[1]][accept.prob]

                                   if(i > rj$mcmc$n.burn) rj$accept[["t.ij"]] <- 
                                     rj$accept[["t.ij"]] + unname(ifelse(accept.prob, 1, 0))
                                   
                                   # rj$iter["t.ij"] <- rj$iter["t.ij"] + 1
                                   
                                   #'------------------------------
                                   # // sigma ----
                                   #'------------------------------
                                   proposed.sigma <- proposal_mh(rj.obj = rj, param.name = "sigma")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.sigma[[1]], p = "sigma")){
                                     accept.prob <- 0
                                   } else {
                                     
                                     loglik.proposed <- likelihood(param.name = "sigma",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.sigma,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(param.name = "sigma",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob){
                                     rj$sigma[i] <- proposed.sigma[[1]]
                                     if(i > rj$mcmc$n.burn) rj$accept[["sigma"]] <- 
                                         rj$accept[["sigma"]] + 1
                                   }
                                   
                                   # rj$iter["sigma"] <- rj$iter["sigma"] + 1
                                   
                                   #'------------------------------
                                   # // mu.i ----
                                   #'------------------------------
                                   proposed.mu.i <- proposal_mh(rj.obj = rj, param.name = "mu.i")
                                   
                                   loglik.proposed <- likelihood(param.name = "mu.i",
                                                                 rj.obj = rj,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = proposed.mu.i,
                                                                 included.cov = NULL,
                                                                 RJ = FALSE,
                                                                 lprod = FALSE)
                                   
                                   loglik.current <- likelihood(param.name = "mu.i", 
                                                                rj.obj = rj,
                                                                model = rj$mlist[[rj$current.model]],
                                                                values = NULL,
                                                                included.cov = NULL,
                                                                RJ = FALSE,
                                                                lprod = FALSE)
                                   
                                   prop.forward <- propdens_mh(rj.obj = rj,
                                                               param.name = "mu.i", 
                                                               dest = proposed.mu.i[[1]],
                                                               orig = rj$mu.i[i - 1, ])
                                   
                                   prop.backward <- propdens_mh(rj.obj = rj,
                                                                param.name = "mu.i", 
                                                                dest = rj$mu.i[i - 1, ],
                                                                orig = proposed.mu.i[[1]])
                                   
                                   accept.prob <- 
                                     runif(rj$dat$whales$n) < 
                                     exp((loglik.proposed + prop.backward) -
                                           (loglik.current + prop.forward))
                                   
                                   # rj$mu.i[i, ] <- rj$mu.i[i - 1, ]
                                   rj$mu.i[i, accept.prob] <- proposed.mu.i[[1]][accept.prob]
                                   
                                   if(i > rj$mcmc$n.burn) 
                                     rj$accept[["mu.i"]] <- 
                                     rj$accept[["mu.i"]] + unname(ifelse(accept.prob, 1, 0))
                                   
                                   # rj$iter["mu.i"] <- rj$iter["mu.i"] + 1
                                   
                                   #'------------------------------
                                   # // mu ----
                                   #'------------------------------
                                   proposed.mu <- proposal_mh(rj.obj = rj, param.name = "mu")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.mu[[1]], p = "mu")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(param = "mu",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.mu,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(param.name = "mu",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob){
                                     rj$mu[i, ] <- proposed.mu[[1]] 
                                     if(i > rj$mcmc$n.burn) rj$accept[["mu"]] <- rj$accept[["mu"]] + 1}
                                   
                                   #'------------------------------
                                   # // phi ----
                                   #'------------------------------
                                   proposed.phi <- proposal_mh(rj.obj = rj, param.name = "phi")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.phi[[1]], p = "phi")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(param.name = "phi",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.phi,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(param.name = "phi",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$phi[i] <- proposed.phi[[1]] 
                                     if(i > rj$mcmc$n.burn) rj$accept[["phi"]] <- rj$accept[["phi"]] + 1
                                   }
                                   
                                   # rj$iter["phi"] <- rj$iter["phi"] + 1
                                   
                                 } else {

                                   if(rj$dat$covariates$n > 0){
                                     
                                     rj$pi.ij[i, ] <- 
                                       pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id] + 
                                               cov_effects(rj.obj = rj))
                                     
                                   } else {
                                     
                                     rj$pi.ij[i, ] <- pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id])
                                     
                                   }
                                   
                                   # rj$iter["pi.ij"] <- rj$iter["pi.ij"] + 1
                                   
                                   #'------------------------------
                                   # // alpha ----
                                   #'------------------------------
                                   proposed.alpha <- proposal_mh(rj.obj = rj, param.name = "alpha")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.alpha[[1]], p = "alpha")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "alpha",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.alpha,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "alpha",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$alpha[i, ] <- proposed.alpha[[1]]
                                     if(i > rj$mcmc$n.burn) rj$accept[["alpha"]] <-
                                         rj$accept[["alpha"]] + 1
                                   }
                                  
                                   
                                   #'------------------------------
                                   # // nu1 ----
                                   #'------------------------------
                                   
                                   proposed.nu <- proposal_mh(rj.obj = rj, param.name = "nu")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.nu[[1]][1, ], p = "nu1")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "nu1",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.nu,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "nu1",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$nu[i, , 1] <- proposed.nu[[1]][1, ]
                                     if(i > rj$mcmc$n.burn) rj$accept[["nu"]][1] <- 
                                         rj$accept[["nu"]][1] + 1
                                   }
                                   
                                   #'------------------------------
                                   # // nu2 ----
                                   #'------------------------------
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.nu[[1]][2, ], p = "nu2")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "nu2",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.nu,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "nu2",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$nu[i, , 2] <- proposed.nu[[1]][2, ]
                                     if(i > rj$mcmc$n.burn) rj$accept[["nu"]][2] <- 
                                         rj$accept[["nu"]][2] + 1
                                   }
                                   
                                   #'------------------------------
                                   # // tau1 ----
                                   #'------------------------------
                                   
                                   proposed.tau <- proposal_mh(rj.obj = rj, param.name = "tau")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.tau[[1]][1], p = "tau")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "tau1",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.tau,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "tau1",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$tau[i, 1] <- proposed.tau[[1]][1]
                                     # rj$tau[i, 2] <- rj$tau[i - 1, 2]
                                     if(i > rj$mcmc$n.burn) rj$accept[["tau"]][1] <-
                                         rj$accept[["tau"]][1] + 1
                                   }
                                   
                                   # rj$iter["tau"] <- rj$iter["tau"] + 1
                                   
                                   #'------------------------------
                                   # // tau2 ----
                                   #'------------------------------
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.tau[[1]][2], p = "tau")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "tau2",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.tau,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "tau2",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$tau[i, 2] <- proposed.tau[[1]][2]
                                     if(i > rj$mcmc$n.burn) rj$accept[["tau"]][2] <-
                                         rj$accept[["tau"]][2] + 1
                                   }
                                   
                                   #'------------------------------
                                   # // omega ----
                                   #'------------------------------
                                   
                                   proposed.omega <- proposal_mh(rj.obj = rj, param.name = "omega")
                                   
                                   if(outOfbounds(rj.obj = rj, v = proposed.omega[[1]], p = "omega")){
                                     
                                     accept.prob <- 0
                                     
                                   } else {
                                     
                                     loglik.proposed <- likelihood(biphasic = TRUE,
                                                                   param.name = "omega",
                                                                   rj.obj = rj,
                                                                   model = rj$mlist[[rj$current.model]],
                                                                   values = proposed.omega,
                                                                   included.cov = NULL,
                                                                   RJ = FALSE,
                                                                   lprod = TRUE)
                                     
                                     loglik.current <- likelihood(biphasic = TRUE,
                                                                  param.name = "omega",
                                                                  rj.obj = rj,
                                                                  model = rj$mlist[[rj$current.model]],
                                                                  values = NULL,
                                                                  included.cov = NULL,
                                                                  RJ = FALSE,
                                                                  lprod = TRUE)
                                     
                                     accept.prob <- exp(loglik.proposed - loglik.current)}
                                   
                                   if (runif(1) < accept.prob) {
                                     rj$omega[i] <- proposed.omega[[1]]
                                     if(i > rj$mcmc$n.burn) rj$accept[["omega"]] <- 
                                         rj$accept[["omega"]] + 1
                                   }
                                   
                                   # rj$iter["omega"] <- rj$iter["omega"] + 1
                                   
                                   #'------------------------------
                                   # // psi ----
                                   #'------------------------------
                                   
                                   # Gibbs sampler as we have a Normal prior and a Normal likelihood
                                   tmp.calc <- rj$dat$whales$n / rj$omega[rj$iter["omega"]] ^ 2
                                   sigma2.post <- 1 / (rj$config$psi.gibbs[1] + tmp.calc) 
                                   mu.post <- sigma2.post * (rj$config$psi.gibbs[2] + 
                                               mean(rj$psi.i[rj$iter["psi.i"], ]) * tmp.calc)
                                   rj$psi[i] <- rnorm(1, mean = mu.post, sd = sqrt(sigma2.post))
                                   # rj$iter["psi"] <- rj$iter["psi"] + 1
                                   
                                   #'------------------------------
                                   # // mu.ij ----
                                   #'------------------------------
                                   
                                   proposed.mu.ij <- proposal_mh(rj.obj = rj, param.name = "mu.ij")
                                   
                                   loglik.proposed <- likelihood(biphasic = TRUE,
                                                                 param.name = "mu.ij",
                                                                 rj.obj = rj,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = proposed.mu.ij,
                                                                 included.cov = NULL,
                                                                 RJ = FALSE,
                                                                 lprod = FALSE)
                                   
                                   loglik.current <- likelihood(biphasic = TRUE,
                                                                param.name = "mu.ij",
                                                                rj.obj = rj,
                                                                model = rj$mlist[[rj$current.model]],
                                                                values = NULL,
                                                                included.cov = NULL,
                                                                RJ = FALSE,
                                                                lprod = FALSE)
                                   
                                   logprop.forward <- propdens_mh(rj.obj = rj, param.name = "mu.ij",
                                                                  dest = proposed.mu.ij[[1]],
                                                                  orig = rj$mu.ij[rj$iter["mu.ij"], ,])
                                   
                                   logprop.backward <- propdens_mh(rj.obj = rj, param.name = "mu.ij",
                                                                   dest = rj$mu.ij[rj$iter["mu.ij"], ,],
                                                                   orig = proposed.mu.ij[[1]])
                                   
                                   accept.prob <-  runif(2 * rj$dat$trials$n) < 
                                     exp((loglik.proposed + logprop.backward) - 
                                           (loglik.current + logprop.forward))
                                   
                                   rj$mu.ij[i, accept.prob[, 1], 1] <- 
                                     proposed.mu.ij[[1]][accept.prob[, 1], 1]
                                   
                                   # rj$mu.ij[i, !accept.prob[, 1], 1] <- 
                                   #   rj$mu.ij[i - 1, !accept.prob[, 1], 1]
                                   
                                   rj$mu.ij[i, accept.prob[, 2], 2] <- 
                                     proposed.mu.ij[[1]][accept.prob[, 2], 2]
                                   
                                   # rj$mu.ij[i, !accept.prob[, 2], 2] <- 
                                   #   rj$mu.ij[i - 1, !accept.prob[, 2], 2]
                                   
                                   if(i > rj$mcmc$n.burn){
                                     rj$accept[["mu.ij.1"]] <- rj$accept[["mu.ij.1"]] + accept.prob[, 1]
                                     rj$accept[["mu.ij.2"]] <- rj$accept[["mu.ij.2"]] + accept.prob[, 2]}
                                   
                                   # rj$iter["mu.ij"] <- rj$iter["mu.ij"] + 1
                                   # rj$iter["t.ij"] <- rj$iter["t.ij"] + 1
                                   
                                   rj$t.ij[i, ] <- ((2 - rj$k.ij[rj$iter["k.ij"], ]) * 
                                                      rj$mu.ij[i, , 1]) +
                                     ((1 - (2 - rj$k.ij[rj$iter["k.ij"], ])) * 
                                        rj$mu.ij[i, , 2])
                                   
                                   #'------------------------------
                                   # // psi.i ----
                                   #'------------------------------
                                   
                                   proposed.psi.i <- proposal_mh(rj.obj = rj, param.name = "psi.i")
                                   
                                   loglik.proposed <- likelihood(biphasic = TRUE,
                                                                 param.name = "psi.i",
                                                                 rj.obj = rj,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = proposed.psi.i,
                                                                 included.cov = NULL,
                                                                 RJ = FALSE,
                                                                 lprod = FALSE)
                                   
                                   loglik.current <- likelihood(biphasic = TRUE,
                                                                param.name = "psi.i", 
                                                                rj.obj = rj,
                                                                model = rj$mlist[[rj$current.model]],
                                                                values = NULL,
                                                                included.cov = NULL,
                                                                RJ = FALSE,
                                                                lprod = FALSE)
                                   
                                   accept.prob <- 
                                     runif(rj$dat$whales$n) < exp(loglik.proposed - loglik.current)
                                   
                                   # rj$psi.i[i, ] <- rj$psi.i[i - 1, ]
                                   rj$psi.i[i, accept.prob] <- proposed.psi.i[[1]][accept.prob]

                                   if(rj$dat$covariates$n > 0){
                                     rj$pi.ij[i, ] <- 
                                       pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id] + 
                                               cov_effects(rj.obj = rj))
                                   } else {
                                     rj$pi.ij[i, ] <- pnorm(rj$psi.i[rj$iter["psi.i"], rj$dat$whales$id])
                                   }
                                   
                                   if(i > rj$mcmc$n.burn) {
                                     rj$accept[["psi.i"]][accept.prob] <- 
                                       rj$accept[["psi.i"]][accept.prob] + 1
                                   }
                                   
                                   # rj$iter["psi.i"] <- rj$iter["psi.i"] + 1

                                   #'------------------------------
                                   # // k.ij ----
                                   #'------------------------------
                                   
                                   proposed.k.ij <- proposal_mh(rj.obj = rj, param.name = "k.ij")
                                   
                                   loglik.proposed <- likelihood(biphasic = TRUE,
                                                                 param.name = "k.ij",
                                                                 rj.obj = rj,
                                                                 model = rj$mlist[[rj$current.model]],
                                                                 values = proposed.k.ij,
                                                                 included.cov = NULL,
                                                                 RJ = FALSE,
                                                                 lprod = FALSE)
                                   
                                   loglik.current <- likelihood(biphasic = TRUE,
                                                                param.name = "k.ij",
                                                                rj.obj = rj,
                                                                model = rj$mlist[[rj$current.model]],
                                                                values = NULL,
                                                                included.cov = NULL,
                                                                RJ = FALSE,
                                                                lprod = FALSE)
                                   
                                   logprop.forward <- d_binom(x = 2 - proposed.k.ij[[1]], 
                                                              size = 1, 
                                                              prob = rj$config$prop$mh$k.ij, 
                                                              log = TRUE)
                                   
                                   logprop.backward <- d_binom(x = 2 - rj$k.ij[i - 1, ], 
                                                               size = 1, 
                                                               prob = rj$config$prop$mh$k.ij, 
                                                               log = TRUE)
                                   
                                   accept.prob <- 
                                     runif(rj$dat$trials$n) < exp((loglik.proposed + logprop.backward) - 
                                                                    (loglik.current + logprop.forward))
                                   
                                   # Set acceptance probability to 0 if proposed values
                                   # do not align with constraints imposed on censored data
                                   rc.check <- rj$mu.ij[rj$iter["mu.ij"], ,][cbind(seq_along(proposed.k.ij[[1]]), proposed.k.ij[[1]])][rj$dat$obs$censored == 1] < rj$dat$obs$Rc[rj$dat$obs$censored == 1]
                                   
                                   lc.check <- rj$mu.ij[rj$iter["mu.ij"], ,][cbind(seq_along(proposed.k.ij[[1]]), proposed.k.ij[[1]])][rj$dat$obs$censored == -1] > rj$dat$obs$Lc[rj$dat$obs$censored == -1]
                                   
                                   accept.prob[rj$dat$obs$censored == 1][rc.check] <- FALSE
                                   accept.prob[rj$dat$obs$censored == -1][lc.check] <- FALSE
                                   
                                   # rj$k.ij[i, ] <- rj$k.ij[i - 1, ]
                                   rj$k.ij[i, accept.prob] <- proposed.k.ij[[1]][accept.prob]
                                   
                                   if(i > rj$mcmc$n.burn) 
                                     rj$accept[["k.ij"]][accept.prob] <- 
                                     rj$accept[["k.ij"]][accept.prob] + 1
                                   
                                   rj$t.ij[i, ] <-  ((2 - rj$k.ij[i, ]) * rj$mu.ij[i, , 1]) +
                                     ((1 - (2 - rj$k.ij[i, ])) * rj$mu.ij[i, , 2])
                                   
                                   # rj$iter["k.ij"] <- rj$iter["k.ij"] + 1
                                   
                                 } # End if biphasic
                                   
                                 if(any(!rj$iter == i)) stop("Mismatched iteration count")
                                 if(rj$config$function.select & sum(rj$alpha[i, ]) == 0 |
                                    rj$config$biphasic & sum(rj$alpha[i, ]) == 0) stop("Zeroes in alpha")
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
