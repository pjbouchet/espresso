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
    last.iter <- purrr::map(.x = seq_along(dat), 
                            .f = ~glance(dat = dat, which.chain = .x, f = "update"))
  } else {
    last.iter <- vector(mode = "list", length = n.chains)
  }
  
  if(do.update){
    dat <- append(dat[[1]][["dat"]], dat[[1]]["config"])
  }
  
  # Parallel computing
  if(n.chains > parallel::detectCores(all.tests = FALSE, logical = TRUE) -1) 
    stop("Insufficient cores.")
  
  isDarwin <- Sys.info()[['sysname']] == "Darwin"
  isWindows <- Sys.info()[['sysname']] == "Windows"
  
  # Create a parallel cluster and register the backend
  cl <- suppressMessages(parallel::makeCluster(n.chains, outfile = ""))
  doParallel::registerDoParallel(cl)
 
  # Progress bar
  if(isDarwin){
  sink(tempfile())
  pb <- utils::txtProgressBar(min = 0, max = (n.iter + n.burn), style = 3)
  sink()
  }
  
  # Simple progress bar
  pb.int <- c(10, 20, 40, 60, 80, 100) * (n.iter + n.burn) / 100

  # Launch the loop on multiple cores if needed
  rj.res <- foreach::foreach(nc = seq_len(n.chains), .packages = "tcltk") %dopar% {
    
    if (isWindows) {
      pb <- tcltk::tkProgressBar(
        title = "rjMCMC", label = "", min = 0,
        max = (n.iter + n.burn),
        initial = 0, width = 300)
    }
    
    rj <- setup_rjMCMC(
        rj.input = dat,
        n.chains = n.chains,
        n.burn = n.burn,
        n.iter = n.iter,
        p.split = dat$config$move$prob[1],
        p.merge = dat$config$move$prob[2],
        do.update = do.update,
        start.values = last.iter[[nc]]
      )
    
    # Set function environments
    environment(cov_effects) <- rlang::current_env()
    environment(normal_prior) <- rlang::current_env()
    environment(uniform_prior) <- rlang::current_env()
    environment(proposal_t.ij) <- rlang::current_env()
    environment(proposal_mu.i) <- rlang::current_env()
    environment(proposal_mu) <- rlang::current_env()
    environment(proposal_alpha) <- rlang::current_env()
    environment(proposal_nu) <- rlang::current_env()
    environment(proposal_mu.ij) <- rlang::current_env()
    environment(propdens_mh) <- rlang::current_env()
    environment(likelihood) <- rlang::current_env()
    environment(propose_jump) <- rlang::current_env()
    environment(proposal_rj) <- rlang::current_env()
    environment(propdens_rj) <- rlang::current_env()
    environment(proposal_ff) <- rlang::current_env()
    environment(propdens_ff) <- rlang::current_env()
    environment(ff_prior) <- rlang::current_env()
    
    # Define all the parameters need for the run
    model.select <- rj$config$model.select
    function.select <- rj$config$function.select
    covariate.select <- rj$config$covariate.select
    biphasic <- rj$config$biphasic
    n.species <- rj$dat$species$n
    n.per.species <- rj$dat$species$nper
    species.names <- rj$dat$species$names
    species.id <- rj$dat$species$id
    n.whales <- rj$dat$whales$n
    whale.id <- rj$dat$whales$id
    y.ij <- rj$dat$obs$y_ij
    obs.sd <- rj$dat$obs$sd
    dummy <- rj$dat$covariates$dummy
    dummy.df <- rj$dat$covariates$dummy.df
    dose.range <- rj$dat$param$dose.range
    species.trials <- rj$dat$species$trials
    n.trials <- rj$dat$trials$n
    n.covariates <- rj$dat$covariates$n
    covariate.names <- rj$dat$covariates$names
    priors <- rj$config$priors
    priors.cov.sd <- rj$config$priors[covariate.names[1], 2]
    prop.covariates <- rj$config$prop$cov
    prop.t.ij <- rj$config$prop$mh$t.ij
    prop.mu.i <- rj$config$prop$mh$mu.i
    prop.mu <- rj$config$prop$mh$mu
    prop.phi <- rj$config$prop$mh$phi
    prop.sigma <- rj$config$prop$mh$sigma
    prop.nu <- rj$config$prop$mh$nu
    prop.tau <- rj$config$prop$mh$tau
    prop.alpha <- rj$config$prop$mh$alpha
    prop.mu.ij <- rj$config$prop$mh$mu.ij
    prop.omega <- rj$config$prop$mh$omega
    prop.psi.i <- rj$config$prop$mh$psi.i
    prop.k.ij <- rj$config$prop$mh$k.ij
    prop.mh <- rj$config$prop$mh
    prop.RJ <- rj$config$prop$rj
    prob.splitmerge <- rj$mcmc$move$prob
    fL <- rj$dat$covariates$fL
    is.censored <- rj$dat$obs$censored
    is.RightCensored <- rj$dat$obs$censored == 1
    is.LeftCensored <- rj$dat$obs$censored == -1
    Rc <- rj$dat$obs$Rc
    Lc <- rj$dat$obs$Lc
    .p.bounds <- matrix(data = dose.range, ncol = 2, nrow = n.trials, byrow = TRUE)
    tot.iter <- rj$mcmc$tot.iter
    dummy.names <- unlist(purrr::map2(.x = covariate.names, 
                                   .y = purrr::map_dbl(.x = fL, "nL"),
                                   .f = ~rep(.x, each = ifelse(.y == 0, 1, .y))))
    one.group <- paste0("(", paste0(species.names, collapse = ","), ")")
    individual.groups <- 
      paste0(paste0("(", species.names, ")"), collapse = "+")
    var.est <- rj$config$var
    
    # Extract values to speed up calculations
    .t.ij <- rj$t.ij[1, ]
    
    if(!biphasic | function.select){
      .sigma <- rj$sigma[1]
      .phi <- rj$phi[1]
      .mu <- rj$mu[1, ]
      .mu.i <- rj$mu.i[1, ]
    }
    
    if(biphasic | function.select){
      .alpha <- rj$alpha[1, ]
      .nu1 <- rj$nu[1, , 1]
      .nu2 <- rj$nu[1, , 2]
      .tau1 <- rj$tau[1, 1]
      .tau2 <- rj$tau[1, 2]
      .omega <- rj$omega[1]
      .psi <- rj$psi[1]
      .mu.ij1 <- rj$mu.ij[1, , 1]
      .mu.ij2 <- rj$mu.ij[1, , 2]
      .psi.i <- rj$psi.i[1, ]
      .k.ij <- rj$k.ij[1, ]
      .pi.ij <- rj$pi.ij[1, ]
    }
    
    .phase <- rj$phase[1]
    .model <- rj$model[1]
    .mlist <- rj$mlist
    .nglist <- rj$nglist
    
    if(n.covariates > 0){
      .include.covariates <- rj$include.covariates[1, ]
      if("exposed" %in% covariate.names){
        .exposed <- rj$exposed[1, ]}
      if("sonar" %in% covariate.names) {
        .sonar <- rj$sonar[1, ]}
      if("behaviour" %in% covariate.names){
        .behaviour <- rj$behaviour[1, ]}
      if("range" %in% covariate.names){
        .range <- rj$range[1, ] }
    }
    
    # Start timer
    start.time <- Sys.time()
    
    for (i in 2:tot.iter) {
      
      if(i %in% pb.int){
        if(isDarwin) utils::setTxtProgressBar(pb, i)
        if(isWindows) tcltk::setTkProgressBar(pb, i)}
      
      #' -----------------------------------------------------
      # Step 1: Split / merge ----
      #' -----------------------------------------------------
      
      if(model.select & n.species > 1) { 
        
        # Propose jump and new parameters
        proposed.jump <- propose_jump(rj.obj = rj)
        new.params <- proposal_rj(jump = proposed.jump, huelsenbeck = FALSE)
        
        # Priors
        logprior.new <-
          uniform_prior(param.name = if (.phase == 1) "mu" else c("nu", "alpha"),
                        param = new.params)
        
        logprior.cur <-
          uniform_prior(param.name = if (.phase == 1) "mu" else c("nu", "alpha"),
            param = if (.phase == 1) {
              list(out = list(mu = .mu))
            } else {
              list(out = list(nu = cbind(.nu1, .nu2), alpha = .alpha))
            }
          )
        
        # Likelihoods
        loglik.new <-
          likelihood(biphasic = ifelse(.phase == 1, FALSE, TRUE),
                     param.name = if(.phase == 1) "mu" else c("nu", "alpha"),
                     rj.obj = rj,
                     values = new.params$out,
                     included.cov = NULL,
                     RJ = TRUE,
                     lprod = TRUE)
        
        loglik.cur <- 
          likelihood(biphasic = ifelse(.phase == 1, FALSE, TRUE),
                     param.name = if(.phase == 1) "mu" else c("nu", "alpha"),
                     rj.obj = rj,
                     values = NULL,
                     included.cov = NULL,
                     RJ = TRUE, 
                     lprod = TRUE)
        
        # Proposal densities
        logpropdens <- 
          propdens_rj(rj.obj = rj,
                      param = new.params,
                      jump = proposed.jump,
                      phase = .phase)
        
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
        
        # Evaluate acceptance ratio
        if (runif(1) < exp(lognum - logden)) {
          
          rj$accept.model[i] <- 1
          
          new.model <- vec_to_model(
            input.vector = proposed.jump$model$id,
            sp.names = species.names)
          
          .model <- new.model
          
          # Add model to list
          if(!new.model %in% names(.mlist)){
            .mlist[[new.model]] <- proposed.jump$model$id
            .nglist[[new.model]] <- dplyr::n_distinct(proposed.jump$model$id)
          }

          # Save new values
          if(.phase == 1){
            .mu <- new.params$out$mu
          } else {
            .nu1 <- new.params$out$nu[, 1]
            .nu2 <- new.params$out$nu[, 2]
            .alpha <- new.params$out$alpha
          }
          
        } 
       
      }
      
      #' -----------------------------------------------------
      # Step 2: Add/remove covariates ----
      #' -----------------------------------------------------
      
      if(n.covariates > 0){
        
        if(covariate.select){
          
          # We consider all models to be equally likely 
          # when it comes to covariates.
          # Therefore, the priors on models cancel out.
          
          # Select a covariate at random and switch it on/off
          chosen.covariate <- sample(x = 1:n.covariates, size = 1)
          proposed.covariates <- .include.covariates
          proposed.covariates[chosen.covariate] <- add.remove <- 
            1 - .include.covariates[chosen.covariate]
          
          if(add.remove){
            
            # Coefficients for covariate terms
            beta.params <- rnorm(n = fL[[chosen.covariate]]$nparam,
                                 mean = 0, sd = prop.covariates)
            
            # Prior
            logprior <- 
              normal_prior(rj.obj = rj,  
                           param.name = covariate.names[chosen.covariate],
                           param = beta.params)
            
            # Proposal density
            logpropdens <- 
              sum(dnorm(x = beta.params, mean = 0, 
                        sd = prop.covariates, log = TRUE))
            
            # Ratio of prior / proposal density
            R <- logprior - logpropdens
            
            if(fL[[chosen.covariate]]$nL > 1) 
              beta.params <- c(0, beta.params)
            
            beta.params <- stats::setNames(list(beta.params), 
                                           covariate.names[chosen.covariate])
            
          } else {
            
            beta.params <- list()
            
            # Prior
            current.value <- 
              get(paste0(".", covariate.names[chosen.covariate]))[fL[[chosen.covariate]]$index]
            
            logprior <- normal_prior(rj.obj = rj,
                                     param.name = covariate.names[chosen.covariate],
                                     param = current.value)
            
            # Proposal density
            logpropdens <- 
              sum(dnorm(x = current.value, mean = 0, 
                        sd = prop.covariates, 
                        log = TRUE))
            
            # Ratio of prior / proposal density
            R <- logpropdens - logprior
            
          }
          
          # Probabilities of addition / removal
          if(add.remove) pmove <- 1/sum(.include.covariates == 0) else 
            pmove <- 1/sum(.include.covariates == 1)
          pmove <- log(c(pmove, 1/sum(proposed.covariates == add.remove)))
          
          # Likelihoods
          loglik.new <- 
            likelihood(biphasic = ifelse(.phase == 1, FALSE, TRUE),
                       param.name = "covariates",
                       rj.obj = rj,
                       values = beta.params,
                       included.cov = proposed.covariates,
                       RJ = TRUE,
                       lprod = TRUE)
          
          loglik.cur <- 
            likelihood(biphasic = ifelse(.phase == 1, FALSE, TRUE),
                       param.name = "covariates",
                       rj.obj = rj,
                       values = NULL,
                       included.cov = .include.covariates, 
                       RJ = TRUE,
                       lprod = TRUE)
          
          # Posterior ratio
          lognum <- loglik.new + pmove[2] + R
          logden <- loglik.cur + pmove[1]
          
          if (runif(1) < exp(lognum - logden)) {
            rj$accept.covariates[i] <- 1
            .include.covariates <- proposed.covariates
            if(add.remove){
              assign(x = paste0(".", covariate.names[chosen.covariate]),
                     value = unlist(beta.params), envir = rlang::current_env())}
          }
        }  # End if (covariate.select)
      } # End if (n.covariates > 0)
      
      #' --------------------------------------
      # Step 3: Functional form ----
      #' --------------------------------------
      
      if(function.select){
        
        # We can generalize the standard 'monophasic' dose-response
        # function by allowing the threshold to come from one of two
        # truncated normal distributions, one with lower exposure 
        # values than the other. This gives a 'biphasic' 
        # function with a context-dependent part and a 
        # dose-dependent part.
        
        ff.prop <- proposal_ff()
        
        # Likelihoods
        loglik.new <-
          likelihood(biphasic = ifelse(.phase == 1, TRUE, FALSE),
                     param.name = if(.phase == 2) c("mu", "mu.i", "phi", "sigma") else 
                       c("nu", "alpha", "tau1", "tau2", "omega", "psi", "mu.ij", "psi.i", "k.ij"),
                     rj.obj = rj,
                     values = ff.prop,
                     included.cov = NULL,
                     RJ = TRUE,
                     lprod = TRUE)
        
        loglik.cur <-
          likelihood(biphasic = ifelse(.phase == 1, FALSE, TRUE),
                     param.name = NULL,
                     rj.obj = rj,
                     values = NULL,
                     included.cov = NULL,
                     RJ = TRUE,
                     lprod = TRUE)
        
        # Priors
        logprior.new <- ff_prior(param = ff.prop, phase = ff.prop$to.phase)
        logprior.cur <- ff_prior(param = NULL, phase = .phase)
        
        # Proposal densities
        logpropdens <- propdens_ff(param = ff.prop)
        
        # Posterior ratio (Jacobian is 1 and p is 1)
        lognum <- loglik.new + logprior.new + logpropdens[[.phase]]
        logden <- loglik.cur + logprior.cur + logpropdens[[ff.prop$to.phase]]
        
        if (runif(1) < exp(lognum - logden)) {
          
          rj$accept.phase[i] <- 1
          .phase <- ff.prop$to.phase
          
          if(.phase == 1){
            
            .sigma <- ff.prop$sigma
            .mu.i <- ff.prop$mu.i
            .mu <- ff.prop$mu
            .phi <- ff.prop$phi
  
          } else {
            
            .alpha <- ff.prop$alpha
            .nu1 <- ff.prop$nu[, 1]
            .nu2 <- ff.prop$nu[, 2]
            .tau1 <- ff.prop$tau[1]
            .tau2 <- ff.prop$tau[2]
            .omega <- ff.prop$omega
            .psi <- ff.prop$psi
            .mu.ij1 <- ff.prop$mu.ij[, 1]
            .mu.ij2 <- ff.prop$mu.ij[, 2]
            .psi.i <- ff.prop$psi.i
            .k.ij <- ff.prop$k.ij
            .pi.ij <- ff.prop$pi.ij
            
          }
        }
      }
      
      #'---------------------------------------------------
      # Step 4: Update parameters ----
      #'---------------------------------------------------
      
      # if(any(!rj$iter == i)) stop("Mismatched iteration count")
      # if(any(rj$t.ij[i, rj$k.ij[i, ] == 1] > rj$alpha[i, species.trials][rj$k.ij[i, ] == 1])) stop("t.ij inconsistent with alpha / k.ij (1)")
      # if(any(rj$t.ij[i, rj$k.ij[i, ] == 2] < rj$alpha[i, species.trials][rj$k.ij[i, ] == 2])) stop("t.ij inconsistent with alpha / k.ij (2)")
      
      #'------------------------------
      ### covariates ----
      #'------------------------------
      
      cov.effects <- 0
      
      if(n.covariates > 0){
        
        #'------------------------------
        # // exposure history ----
        #'------------------------------
        
        if("exposed" %in% covariate.names &
           .include.covariates["exposed"] == 1){
          
          prop.exposed <-
            c(0, rnorm(n = 1,
                       mean = .exposed[2],
                       sd = prop.mh[["exposed"]]))
          
          if(.phase == 1){
            
            loglik.proposed <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects("exposed", 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
            loglik.current <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects(phase = 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
          } else {
            
            pi.ij.proposed <- 
              pnorm(.psi.i[whale.id] + cov_effects("exposed", 2))
            
            pi.ij.current <- 
              pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
            
            loglik.proposed <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.proposed, 
                                       do_log = TRUE))
            
            loglik.current <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.current, 
                                       do_log = TRUE))
            
          }
          
          
          logprior.proposed <-
            dnorm(x = prop.exposed[2], mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          logprior.current <-
            dnorm(x = .exposed[2], mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          accept.prob <- 
            runif(1) < exp((loglik.proposed + logprior.proposed) - 
                             (loglik.current + logprior.current))
          
          if (accept.prob) .exposed <- prop.exposed
          
        }
        
        #'------------------------------
        # // sonar signal ----
        #'------------------------------
        
        if("sonar" %in% covariate.names &
           .include.covariates["sonar"] == 1) {
          
          prop.sonar <- 
            c(0, rnorm(n = 1,
                       mean = .sonar[2],
                       sd = prop.mh[["sonar"]]))
          
          if(.phase == 1){
            
            loglik.proposed <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects("sonar", 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
            loglik.current <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects(phase = 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
          } else {
            
            pi.ij.proposed <- 
              pnorm(.psi.i[whale.id] + cov_effects("sonar", 2))
            
            pi.ij.current <- 
              pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
            
            loglik.proposed <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.proposed, 
                                       do_log = TRUE))
            
            loglik.current <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.current, 
                                       do_log = TRUE))
            
          }
          
          logprior.proposed <- 
            dnorm(x = prop.sonar[2], mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          logprior.current <- 
            dnorm(x = .sonar[2], mean = 0,
                  sd = priors["sonar", 2], 
                  log = TRUE)
          
          accept.prob <- 
            runif(1) < exp((loglik.proposed + logprior.proposed) - 
                             (loglik.current + logprior.current))
          
          if (accept.prob) .sonar <- prop.sonar
          
        }
        
        #'------------------------------
        # // behavioural mode ----
        #'------------------------------
        
        if("behaviour" %in% covariate.names &
           .include.covariates["behaviour"] == 1) {
          
          prop.behaviour <- 
            c(0, rnorm(n = 1,
                       mean = .behaviour[2],
                       sd = prop.mh[["behaviour"]]))
          
          if(.phase == 1){
            
            loglik.proposed <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects("behaviour", 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
            loglik.current <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects(phase = 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
          } else {
            
            pi.ij.proposed <- 
              pnorm(.psi.i[whale.id] + cov_effects("behaviour", 2))
            
            pi.ij.current <- 
              pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
            
            loglik.proposed <- 
              sum(d_Binom(x = 2 - .k.ij, 
                          size = 1, 
                          prob = pi.ij.proposed, 
                          do_log = TRUE))
            
            loglik.current <- 
              sum(d_Binom(x = 2 - .k.ij, 
                          size = 1, 
                          prob = pi.ij.current, 
                          do_log = TRUE))
            
          }
          
          logprior.proposed <- 
            dnorm(x = prop.behaviour[2], mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          logprior.current <- 
            dnorm(x = .behaviour[2], mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          accept.prob <- 
            runif(1) < exp((loglik.proposed + logprior.proposed) - 
                             (loglik.current + logprior.current))
          
          if (accept.prob) .behaviour <- prop.behaviour
        }
        
        #'------------------------------
        # // Source-whale range ----
        #'------------------------------
        
        if("range" %in% covariate.names & 
           .include.covariates["range"] == 1) {
          
          prop.range <- 
            rnorm(n = 1, mean = .range, sd = prop.mh[["range"]])
          
          if(.phase == 1){
            
            loglik.proposed <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects("range", 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
            loglik.current <- 
              sum(dt_norm(
                x = .t.ij,
                location = .mu.i[whale.id] + cov_effects(phase = 1),
                scale = .sigma,
                L = dose.range[1],
                U = dose.range[2],
                do_log = TRUE))
            
          } else {
            
            pi.ij.proposed <- 
              pnorm(.psi.i[whale.id] + cov_effects("range", 2))
            
            pi.ij.current <- 
              pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
            
            loglik.proposed <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.proposed, 
                                       do_log = TRUE))
            
            loglik.current <- 
              sum(d_Binom(x = 2 - .k.ij, 
                                       size = 1, 
                                       prob = pi.ij.current, 
                                       do_log = TRUE))
          }
          
          logprior.proposed <- 
            dnorm(x = prop.range, mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          logprior.current <- 
            dnorm(x = .range, mean = 0, sd = priors.cov.sd,
                  log = TRUE)
          
          accept.prob <- runif(1) < 
            exp((loglik.proposed + logprior.proposed) - 
                  (loglik.current + logprior.current))
          
          if (accept.prob) .range <- prop.range
          
        }
        
      } # End covariates.n > 0
      
      #'------------------------------
      ### Monophasic ----
      #'------------------------------
      
      if(.phase == 1){
        
        #'------------------------------
        ### t.ij (monophasic) ----
        #'------------------------------
        
        proposed.t.ij <- proposal_t.ij()
        
        loglik.proposed <- 
          dnorm(x = y.ij, mean = proposed.t.ij[[1]], 
                sd = obs.sd, log = TRUE)
        
        loglik.proposed[is.na(loglik.proposed)] <- 0 
        
        loglik.proposed <- loglik.proposed + 
          dt_norm(
            x = proposed.t.ij[[1]],
            location = .mu.i[whale.id] + cov_effects(phase = 1),
            scale = .sigma,
            L = dose.range[1],
            U = dose.range[2],
            do_log = TRUE)
        
        loglik.current <- 
          dnorm(x = y.ij, mean = .t.ij, sd = obs.sd, log = TRUE)
        
        loglik.current[is.na(loglik.current)] <- 0 
        
        loglik.current <- loglik.current + 
          dt_norm(
            x = .t.ij,
            location = .mu.i[whale.id] + cov_effects(phase = 1),
            scale = .sigma,
            L = dose.range[1],
            U = dose.range[2],
            do_log = TRUE)
        
        prop.forward <- 
          propdens_mh(rj.obj = rj,
                      param.name = "t.ij", 
                      dest = proposed.t.ij[[1]], 
                      orig = .t.ij,
                      bounds = proposed.t.ij[[2]])
        
        prop.backward <- 
          propdens_mh(rj.obj = rj,
                      param.name = "t.ij", 
                      dest = .t.ij,
                      orig = proposed.t.ij[[1]],
                      bounds = proposed.t.ij[[2]])
        
        accept.prob <- runif(n.trials) < 
          exp((loglik.proposed + prop.backward) - 
                (loglik.current + prop.forward))
        
        .t.ij[accept.prob] <- proposed.t.ij[[1]][accept.prob]
        
        #'------------------------------
        # // sigma ----
        #'------------------------------
        
        proposed.sigma <- rnorm(n = 1, mean = .sigma, sd = prop.sigma)
        
        if(proposed.sigma < priors["sigma", 1] | 
           proposed.sigma > priors["sigma", 2]){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .t.ij,
              location = .mu.i[whale.id] + cov_effects(phase = 1),
              scale = proposed.sigma,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          loglik.current <- 
            sum(dt_norm(
              x = .t.ij,
              location = .mu.i[whale.id] + cov_effects(phase = 1),
              scale = .sigma,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          accept.prob <- exp(loglik.proposed - loglik.current)
        }
        
        if (runif(1) < accept.prob) .sigma <- proposed.sigma
        
        #'------------------------------
        # // mu.i ----
        #'------------------------------
        
        proposed.mu.i <- proposal_mu.i()
        
        loglik.proposed <- dt_norm(
          x = .t.ij,
          location = proposed.mu.i[whale.id] + cov_effects(phase = 1),
          scale = .sigma,
          L = dose.range[1],
          U = dose.range[2],
          do_log = TRUE)
        
        loglik.proposed <- 
          purrr::map_dbl(.x = seq_len(n.whales), 
                         .f = ~ sum(loglik.proposed[whale.id == .x]))
        
        loglik.proposed <- loglik.proposed + 
          dt_norm(
            x = proposed.mu.i,
            location = .mu[species.id],
            scale = .phi,
            L = dose.range[1],
            U = dose.range[2],
            do_log = TRUE)
        
        loglik.current <- dt_norm(
          x = .t.ij,
          location = .mu.i[whale.id] + cov_effects(phase = 1),
          scale = .sigma,
          L = dose.range[1],
          U = dose.range[2],
          do_log = TRUE)
        
        loglik.current <- 
          purrr::map_dbl(.x = seq_len(n.whales), 
                         .f = ~ sum(loglik.current[whale.id == .x]))
        
        loglik.current <- loglik.current + 
          dt_norm(
            x = .mu.i,
            location = .mu[species.id],
            scale = .phi,
            L = dose.range[1],
            U = dose.range[2],
            do_log = TRUE)
        
        prop.forward <- 
          propdens_mh(rj.obj = rj,
                      param.name = "mu.i", 
                      dest = proposed.mu.i,
                      orig = .mu.i)
        
        prop.backward <- 
          propdens_mh(rj.obj = rj,
                      param.name = "mu.i", 
                      dest = .mu.i,
                      orig = proposed.mu.i)
        
        accept.prob <- 
          runif(n.whales) < 
          exp((loglik.proposed + prop.backward) -
                (loglik.current + prop.forward))
        
        .mu.i[accept.prob] <- proposed.mu.i[accept.prob]
        
        #'------------------------------
        # // mu ----
        #'------------------------------
        
        proposed.mu <- proposal_mu()
        
        if(any(proposed.mu < priors["mu", 1]) | 
           any(proposed.mu > priors["mu", 2])) {
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.i,
              location = proposed.mu[species.id],
              scale = .phi,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.i,
              location = .mu[species.id],
              scale = .phi,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        
        if (runif(1) < accept.prob) .mu <- proposed.mu
        
        #'------------------------------
        # // phi ----
        #'------------------------------
        
        proposed.phi <- rnorm(n = 1, mean = .phi, sd = prop.phi)
        
        if(proposed.phi < priors["phi", 1] |
           proposed.phi > priors["phi", 2]){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.i,
              location = .mu[species.id],
              scale = proposed.phi,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.i,
              location = .mu[species.id],
              scale = .phi,
              L = dose.range[1],
              U = dose.range[2],
              do_log = TRUE))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        
        if (runif(1) < accept.prob) .phi <- proposed.phi
        
      } 
      
      #'------------------------------
      ### Biphasic ----
      #'------------------------------
      
      if(.phase == 2){
        
        if(n.covariates > 0){
          .pi.ij <- pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
        } else {
          .pi.ij <- pnorm(.psi.i[whale.id])
        }
        
        #'------------------------------
        # // alpha ----
        #'------------------------------
        
        proposed.alpha <- proposal_alpha()
        
        if(any(proposed.alpha < priors["alpha", 1]) |
           any(proposed.alpha > priors["alpha", 2]) |
           any(proposed.alpha < .nu1) |
           any(proposed.alpha > .nu2)){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <-
            sum(dt_norm(
              x = .mu.ij1,
              location = .nu1[species.trials],
              scale = .tau1,
              L = dose.range[1],
              U = proposed.alpha[species.trials], 
              do_log = TRUE
            )) +
            sum(dt_norm(
              x = .mu.ij2,
              location = .nu2[species.trials],
              scale = .tau2,
              L = proposed.alpha[species.trials], 
              U = dose.range[2],
              do_log = TRUE
            ))
          
          loglik.current <-
            sum(dt_norm(
              x = .mu.ij1,
              location = .nu1[species.trials],
              scale = .tau1,
              L = dose.range[1],
              U = .alpha[species.trials],
              do_log = TRUE
            )) +
            sum(dt_norm(
              x = .mu.ij2,
              location = .nu2[species.trials],
              scale = .tau2,
              U = dose.range[2],
              L = .alpha[species.trials],
              do_log = TRUE
            ))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .alpha <- proposed.alpha
        
        #'------------------------------
        # // nu1 ----
        #'------------------------------
        
        proposed.nu <- proposal_nu()
        
        if(any(proposed.nu[ ,1] < priors["nu", 1]) | 
           any(proposed.nu[ ,1] > .alpha)){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.ij1,
              location = proposed.nu[species.trials, 1],
              scale = .tau1,
              L = dose.range[1],
              U = .alpha[species.trials], 
              do_log = TRUE
            ))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.ij1,
              location = .nu1[species.trials],
              scale = .tau1,
              L = dose.range[1],
              U = .alpha[species.trials], 
              do_log = TRUE
            ))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .nu1 <- proposed.nu[, 1]
        
        
        #'------------------------------
        # // nu2 ----
        #'------------------------------
        
        if(any(proposed.nu[ ,2] > priors["nu", 2]) | 
           any(proposed.nu[ ,2] < .alpha)){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.ij2,
              location = proposed.nu[species.trials, 2],
              scale = .tau2,
              U = dose.range[2],
              L = .alpha[species.trials],
              do_log = TRUE
            ))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.ij2,
              location = .nu2[species.trials],
              scale = .tau2,
              U = dose.range[2],
              L = .alpha[species.trials],
              do_log = TRUE
            ))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .nu2 <- proposed.nu[, 2]
        
        #'------------------------------
        # // tau1 ----
        #'------------------------------
        
        proposed.tau <- rnorm(n = 2, mean = c(.tau1, .tau2), sd = prop.tau)
        
        if(proposed.tau[1] < priors["tau", 1] |
           proposed.tau[1] > priors["tau", 2]) {
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.ij1,
              location = .nu1[species.trials],
              scale = proposed.tau[1],
              L = dose.range[1],
              U = .alpha[species.trials], 
              do_log = TRUE
            ))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.ij1,
              location = .nu1[species.trials],
              scale = .tau1,
              L = dose.range[1],
              U = .alpha[species.trials], 
              do_log = TRUE
            ))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .tau1 <- proposed.tau[1]
        
        #'------------------------------
        # // tau2 ----
        #'------------------------------
        
        if(proposed.tau[2] < priors["tau", 1] |
           proposed.tau[2] > priors["tau", 2]) {
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dt_norm(
              x = .mu.ij2,
              location = .nu2[species.trials],
              scale = proposed.tau[2],
              U = dose.range[2],
              L = .alpha[species.trials],
              do_log = TRUE
            ))
          
          loglik.current <- 
            sum(dt_norm(
              x = .mu.ij2,
              location = .nu2[species.trials],
              scale = .tau2,
              U = dose.range[2],
              L = .alpha[species.trials],
              do_log = TRUE
            ))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .tau2 <- proposed.tau[2]
        
        #'------------------------------
        # // omega ----
        #'------------------------------
        
        proposed.omega <- rnorm(n = 1, mean = .omega, sd = prop.omega)
        
        if(proposed.omega < priors["omega", 1] | 
           proposed.omega > priors["omega", 2]){
          
          accept.prob <- 0
          
        } else {
          
          loglik.proposed <- 
            sum(dnorm(x = .psi.i, 
                      mean = .psi, 
                      sd = proposed.omega, 
                      log = TRUE))
          
          loglik.current <- 
            sum(dnorm(x = .psi.i, 
                      mean = .psi, 
                      sd = .omega, 
                      log = TRUE))
          
          accept.prob <- exp(loglik.proposed - loglik.current)}
        if (runif(1) < accept.prob) .omega <- proposed.omega
        
        #'------------------------------
        # // psi ----
        #'------------------------------
        
        # Gibbs sampler as we have a Normal prior and a Normal likelihood
        
        sigma2.post <- 1 / (rj$config$psi.gibbs[1] + (n.whales / .omega^2)) 
        mu.post <- sigma2.post * (rj$config$psi.gibbs[2] + 
                                    mean(.psi.i) * (n.whales / .omega ^ 2))
        .psi <- rnorm(1, mean = mu.post, sd = sqrt(sigma2.post))
        
        #'------------------------------
        # // mu.ij ----
        #'------------------------------
        
        proposed.mu.ij <- proposal_mu.ij()
        
        loglik.proposed <- 
          likelihood(biphasic = TRUE,
                     param.name = "mu.ij",
                     rj.obj = rj,
                     values = proposed.mu.ij["mu.ij"],
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        loglik.current <- 
          likelihood(biphasic = TRUE,
                     param.name = "mu.ij",
                     rj.obj = rj,
                     values = NULL,
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        logprop.forward <- 
          propdens_mh(rj.obj = rj, 
                      param.name = "mu.ij",
                      dest = proposed.mu.ij[[1]],
                      orig = cbind(.mu.ij1, .mu.ij2),
                      bounds = proposed.mu.ij$p.bounds)
        
        logprop.backward <- 
          propdens_mh(rj.obj = rj, 
                      param.name = "mu.ij",
                      dest = cbind(.mu.ij1, .mu.ij2),
                      orig = proposed.mu.ij[[1]],
                      bounds = proposed.mu.ij$p.bounds)
        
        accept.prob <-  runif(2 * n.trials) < 
          exp((loglik.proposed + logprop.backward) - 
                (loglik.current + logprop.forward))
        
        .mu.ij1[accept.prob[, 1]] <- 
          proposed.mu.ij[[1]][accept.prob[, 1], 1]
        
        .mu.ij2[accept.prob[, 2]] <- 
          proposed.mu.ij[[1]][accept.prob[, 2], 2]
        
        .t.ij <- ((2 - .k.ij) * .mu.ij1) + ((1 - (2 - .k.ij)) * .mu.ij2)
        
        #'------------------------------
        # // psi.i ----
        #'------------------------------
        
        proposed.psi.i <- rnorm(n = n.whales, mean = .psi.i, sd = prop.psi.i)
        
        loglik.proposed <- 
          likelihood(biphasic = TRUE,
                     param.name = "psi.i",
                     rj.obj = rj,
                     values = proposed.psi.i,
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        loglik.current <- 
          likelihood(biphasic = TRUE,
                     param.name = "psi.i", 
                     rj.obj = rj,
                     values = NULL,
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        accept.prob <- 
          runif(n.whales) < exp(loglik.proposed - loglik.current)
        
        .psi.i[accept.prob] <- proposed.psi.i[accept.prob]
        
        if(n.covariates > 0){
          .pi.ij <- 
            pnorm(.psi.i[whale.id] + cov_effects(phase = 2))
        } else {
          .pi.ij <- pnorm(.psi.i[whale.id])
        }
        
        #'------------------------------
        # // k.ij ----
        #'------------------------------
        
        proposed.k.ij <- (1 - rbinom(n = n.trials, size = 1, prob = prop.k.ij)) + 1
        
        loglik.proposed <- 
          likelihood(biphasic = TRUE,
                     param.name = "k.ij",
                     rj.obj = rj,
                     values = proposed.k.ij,
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        loglik.current <- 
          likelihood(biphasic = TRUE,
                     param.name = "k.ij",
                     rj.obj = rj,
                     values = NULL,
                     included.cov = NULL,
                     RJ = FALSE,
                     lprod = FALSE)
        
        logprop.forward <- 
          d_Binom(x = 2 - proposed.k.ij, 
                               size = 1, 
                               prob = prop.k.ij, 
                               do_log = TRUE)
        
        logprop.backward <-
          d_Binom(x = 2 - .k.ij, 
                               size = 1, 
                               prob = prop.k.ij, 
                               do_log = TRUE)
        
        accept.prob <-
          runif(n.trials) <
          exp((loglik.proposed + logprop.backward) -
                (loglik.current + logprop.forward))
        
        # Set acceptance probability to 0 if proposed values
        # do not align with constraints imposed on censored data
        rc.check <- cbind(.mu.ij1, .mu.ij2)[cbind(seq_along(proposed.k.ij), proposed.k.ij)][is.RightCensored] < Rc[is.RightCensored]
        
        lc.check <- cbind(.mu.ij1, .mu.ij2)[cbind(seq_along(proposed.k.ij), proposed.k.ij)][is.LeftCensored] > Lc[is.LeftCensored]
        
        accept.prob[is.RightCensored][rc.check] <- FALSE
        accept.prob[is.LeftCensored][lc.check] <- FALSE
        
        .k.ij[accept.prob] <- proposed.k.ij[accept.prob]
        
        .t.ij <- ((2 - .k.ij) * .mu.ij1) + ((1 - (2 - .k.ij)) * .mu.ij2)
        
      } # End if biphasic
      
      #'------------------------------
      ### Save values ----
      #'------------------------------
      
      rj$t.ij[i, ] <- .t.ij
      
      if(function.select | !biphasic){
        rj$sigma[i] <- .sigma
        rj$phi[i] <- .phi
        rj$mu[i, ] <- .mu
        rj$mu.i[i, ] <- .mu.i
      }
      
      if(function.select | biphasic){
        rj$alpha[i, ] <- .alpha
        rj$nu[i, , 1] <- .nu1
        rj$nu[i, , 2] <- .nu2
        rj$tau[i, 1] <- .tau1
        rj$tau[i, 2] <- .tau2
        rj$omega[i] <- .omega
        rj$psi[i] <- .psi
        rj$mu.ij[i, , 1] <- .mu.ij1
        rj$mu.ij[i, , 2] <- .mu.ij2
        rj$psi.i[i, ] <- .psi.i
        rj$k.ij[i, ] <- .k.ij
        rj$pi.ij[i, ] <- .pi.ij
      }
      
      if(n.covariates > 0){
        for (k in covariate.names) rj[[k]][i, ] <- get(paste0(".", k))
        rj$include.covariates[i, ] <- .include.covariates
      }
      
      rj$phase[i] <- .phase
      rj$model[i] <- .model
      
    } # End RJMCMC
    
    rj$mlist <- .mlist
    rj$nglist <- .nglist
    
    end.time <- Sys.time()
    rj$run_time <- 
      hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
    
    rj
    
  } # End foreach
  
  # When the function terminates
  on.exit({
    parallel::stopCluster(cl = cl) # Stop cluster
  })
  
  class(rj.res) <- c("rjmcmc", class(rj.res))
  return(rj.res)
  
}