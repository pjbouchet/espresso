#' Posterior inference for dose-response models
#'
#' Fits the Bayesian hierarchical dose-response model of \insertCite{Miller2014;textual}{espresso} to multiple species using \code{\link[rjags:rjags-package]{rjags}} \insertCite{Plummer2019}{espresso}, and estimates posterior model probabilities using a Gibbs Variable Selection (GVS) approach \insertCite{OHara2009}{espresso}.
#' 
#' @details Adapted from original code developed by Dina Sadykova as part of the \href{https://synergy.st-andrews.ac.uk/mocha/}{Mocha} project. The function can accommodate species/species groups either as a fixed or a random effect.
#' 
#' @export
#' @importFrom Rdpack reprompt
#' @param dat Input data. Must be an object of class \code{rjtrace} or \code{brsdata}.
#' @param random.effects Logical. When \code{TRUE}, uses a random effect model formulation.
#'@param pseudo.n Number of iterations for the pseudo-priors.
#' @param mcmc.n Number of posterior samples.
#' @param n.chains Number of MCMC chains.
#' @param thin Thinning interval.
#' @param epsilon.upper Upper bound on the ε parameter used in the random effect model formulation.
#' 
#' @return A list object of class \code{gvs}.
#' @references
#' \insertAllCited{}
#' @author Phil J. Bouchet
#' @seealso \code{\link{summary.gvs}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Simulate data for two species
#' mydat <- simulate_data(n.species = 2, 
#'                        n.whales = 16, 
#'                        max.trials = 3, 
#'                        covariates = list(exposed = c(0, 5), range = 0.5),
#'                        mu = c(101, 158), 
#'                        phi = 20, 
#'                        sigma = 20, 
#'                        Rc = c(210, 211), 
#'                        seed = 58697)
#' summary(mydat)
#'      
#' # Model selection by GVS                        
#' gvs <- gibbs(dat = mydat, 
#'              random.effects = FALSE, 
#'              include.covariates = FALSE, 
#'              mcmc.n = 1000, 
#'              burnin = 500)
#' }
#' @keywords brs gvs dose-response

gibbs <- function(dat,
                  random.effects = FALSE,
                  pseudo.n = 10000,
                  mcmc.n = 1000,
                  burnin = 1000,
                  n.chains = 1,
                  thin = 1,
                  epsilon.upper = 30)
{
  
  
  #' ---------------------------------------------
  # Start timer
  #' ---------------------------------------------
  tictoc::tic()
  
  #' ---------------------------------------------
  # Perform function checks and initial setup
  #' ---------------------------------------------
  if(dat$covariates$n > 0) include.covariates <- TRUE else include.covariates <- FALSE
  fL <- dat$covariates$fL
  if("config" %in% class(dat)) config <- dat$config else config <- NULL
  
  # Covariate selection not implemented
  covariate.select <- FALSE
  
  if("rjtrace" %in% class(dat)) {
    
    mcmc.n <- dat$mcmc$n.iter
    burnin <- dat$mcmc$n.burn
    n.chains <- dat$mcmc$n.chains
    mcmc.params <- dat$mcmc[c("n.chains", "n.burn", "n.iter", "tot.iter", "iter.rge", "thin")]
    dat <- dat$dat
    
  } else {
    
    mcmc.params <- list(n.chains = n.chains,
                        n.burn = burnin,
                        n.iter = mcmc.n,
                        tot.iter = burnin + mcmc.n,
                        iter.rge = paste0(burnin + 1, ":", burnin + mcmc.n),
                        thin = thin)
    
    
  }
  
  if(!"rjtrace" %in% class(dat) & !"brsdata" %in% class(dat)) 
    stop("Input data must be of class <brsdata> or <rjtrace>")
  
  covariates.names <- dat$covariates$names
  N.species <- dat$species$n
  N.whales <- dat$whales$n
  N.trials <- dat$trials$n
  species.names <- dat$species$names
  species.ID <- dat$species$id
  whale.ID <- dat$whales$id
  y <- dat$obs$y_ij
  y.precision <- dat$obs$prec
  lowerbound <- dat$param$bounds["mu", 1]
  upperbound <- dat$param$bounds["mu", 2]
  phi.upper <- dat$param$bounds["phi", 2]
  sigma.upper <- dat$param$bounds["sigma", 2]
  I.censored <- dat$obs$censored
  Rc <- dat$obs$Rc
  simulation <- dat$param$sim
  
  Rc[I.censored == 0] <- upperbound
  
  if(include.covariates){
    
    I.covariates <- purrr::map(.x = covariates.names, 
                               .f = ~ dat$covariates$dummy[[.x]][, fL[[.x]]$index, drop = FALSE])
    
    I.covariates <- purrr::map(.x = I.covariates, 
                               .f = ~{
                                 tmp.names <- colnames(.x)
                                 tmp.res <- split(.x, rep(1:ncol(.x), each = nrow(.x)))
                                 names(tmp.res) <- tmp.names
                                 tmp.res}) %>% purrr::flatten(.)
    
    # JAGS does not like dashes
    names(I.covariates) <- gsub(pattern = "-", replacement = "_", x = names(I.covariates))
  }
  
  
  #' ---------------------------------------------
  # Species indicator
  #' ---------------------------------------------
  I.species0 <- diag(N.species)
  if (N.species == 2) I.species <- I.species0
  
  if(random.effects) n.comb <- N.species - 1
  
  #' ---------------------------------------------
  # Model index
  #' ---------------------------------------------
  n.mod0 <- rep(1, N.species)
  if (N.species == 2) n.mod <- c(n.mod0, 2)
  
  if(random.effects) n.comb.run.species <- floor(N.species / 2)
  
  #' ---------------------------------------------
  # Index martrices for species
  #' ---------------------------------------------
  
  # The next code section ﬁnds a complete set of competing models for the given dataset
  # using the ind.m.fun function. I.species is an index matrix for species of interest.
  # N.modelspecies is the associated number of competing models. pi.species is a list of
  # probabilities for the categorical distribution. Gamma.n and theta.n are 
  # categorical variables.
  
  # Beta_3_i are parameters governing the effects of different species/species group.
  # Here, N.beta_3 = 3, as there is one beta per species (when treated separately) + 1 when grouped together.
  # I.species is an indicator function that takes the value 1 for species/species group and 0 otherwise.
  
  # ind.m.fun returns a list with the following elements:
  # [[1]]: Species groupings. 
  # Species are shown as columns, and rows represents combination of species/species groups.
  # Within a row, all species given a value of 1 are grouped together.
  # [[2]]: Model number. 
  
  I.species <- ind.m.fun(N.species)[[1]]
  N.beta_3 <- dim(I.species)[1]
  
  # Index has model index in first column and beta index in second
  index <- matrix(NA, ncol = 2, nrow = N.beta_3)
  index[,1] <- ind.m.fun(N.species)[[2]]
  index[,2] <- c(1:N.beta_3)
  
  if(random.effects){
    
    species_id <- matrix(0, ncol = length(species.ID), nrow = max(ind.m.fun(N.species)[[2]]))
    species_id[1, ] <- species.ID
    
    for (i in 2:max(index[, 1])) {
      if (all(I.species[index[, 1] == i, ] == 1)) {
        species_id[i, ] <- 1
      } else {
        species.id.temp <- species.ID
        for (j in 1:dim(I.species[index[, 1] == i, ])[1]) {
          ind.temp <- I.species[index[, 1] == i, ][j, ]
          if (sum(ind.temp) > 1) {
            ind.w <- which(ind.temp == 1)
            for (jj in 2:length(ind.w)) {
              species.id.temp[species.id.temp == ind.w[jj]] <-
                species.id.temp[species.id.temp == ind.w[1]][1]
            }
          }
        }
        species_id[i, ] <- species.id.temp
      }
    }
  }
  
  #'--------------------------------------------------------------------
  # Number of competing models + species groups
  #'--------------------------------------------------------------------
  
  N.modelspecies <- max(ind.m.fun(N.species)[[2]])
  species.Groups <- get_groupings(species.matrix = I.species, 
                                  model.index = index, 
                                  species.names = species.names,
                                  simulation = simulation)
  
  #'--------------------------------------------------------------------
  # Prior probabilities for each model
  #'--------------------------------------------------------------------
  # Each model is believed to be equally likely a priori
  pi.species <- rep(1 / N.modelspecies, N.modelspecies)
  
  #'--------------------------------------------------------------------
  # Matrix of indices (species)
  #'--------------------------------------------------------------------
  
  # This shows:
  # Column 1 - ID of competing model
  # Column 2 to 3 - Row numbers (min to max) from index corresponding to competing model
  ind.matr <- matrix(NA, ncol = 3, nrow = N.modelspecies)
  ind.matr[, 1] <- c(1:N.modelspecies)
  for (k in 1:N.modelspecies){
    ind.matr[k, 2] <- min(which(index[, 1]==k))
    ind.matr[k, 3] <- max(which(index[, 1]==k))
  }
  
  #'--------------------------------------------------------------------
  # Parameterisation
  #'--------------------------------------------------------------------
  
  # Theta is a categorical variable for species
  theta.n <- c(1:N.modelspecies)
  
  # || Pseudo-priors ----
  
  # In Bayesian model/variable selection (i.e. when using an indicator variable to compute the 
  # posterior selection probability of a given model/variable), it is common for the MCMC chain 
  # for that index to be highly autocorrelated, as it gets stuck on a particular candidate model 
  # or candidate variable. This is because values for all parameters are drawn at each step of the
  # MCMC chain, yet only some of these parameters are used to describe the data, depending on the 
  # value of the index. The remaining parameters thus float 'unconstrained' by the data and are 
  # sampled randomly from their priors. The resulting proposed values might be far away from posterior
  # credible values, such that the chain rarely jumps to them, instead lingering on the current model/variable. Pseudopriors afford a solution to this problem, by ensuring that draws are made within a zone of 
  # posterior credibility. In essence, pseudopriors are used to mimic the posterior when the parameter 
  # is not being used to describe the data. In this sense, a pseudoprior is not really a prior, but only 
  # a conveniently chosen linking density, required to completely define the joint model specification. 
  
  # Create a results matrix to hold posterior estimates for beta pseudo-priors
  
  if(random.effects) { results.species <- array(0, dim = c(N.species, 2, N.modelspecies)) 
  } else { results.species <- matrix(0, nrow = N.beta_3, ncol = N.species) }
  
  
  if(random.effects){
    
    jags.pseudo <- "model{

  # ------------------------------
  # Priors
  # ------------------------------
  
  # On average whale threshold
  mu ~ dunif(lowerbound, upperbound)
  
  # Between-whale variation (sd) 
  phi ~ dunif(0, phi.upper)
  inv_phi2 <- pow(phi, -2)
  
  # Within-whale variation (sd)
  sigma ~ dunif(0, sigma.upper)
  inv_sigma2 <- pow(sigma, -2)
  
  # Between-species variation (sd)
  epsilon ~ dunif(0, epsilon.upper) 
  inv_epsilon2 <- pow(epsilon, -2)
  
  # Priors on covariates
  
  # Index for competing species models 
  for (l in 1:N.modelspecies) {theta_1[l] <- step(l-theta)*step(theta-l)}

  # ------------------------------
  # Process model
  # ------------------------------

  # Between species
  for(i in 1:N.species){
    for (kk in 1:N.modelspecies){
      mu_i[kk,i] ~ dnorm(mu, inv_epsilon2) T(lowerbound, upperbound)    
    }
  }

  for(j in 1:N.whales){
	for (kk in 1:N.modelspecies){
	  v[kk,j] <- mu_i[kk, species_ID[kk, j]]*theta_1[kk]
	}
    
  mu_ij[j] <- sum(v[1:N.modelspecies,j])

  # Between whales within species
    t_ij[j] ~ dnorm(mu_ij[j], inv_phi2) T(lowerbound, upperbound)   
  }

  # Between trials within whale
  for(k in 1:N.trials){
  mu_ijk[k]<- t_ij[whale.ID[k]] #/#

  # Realized threshold
    t[k] ~ dnorm(mu_ijk[k], inv_sigma2) T(lowerbound, upperbound)    
  }

  # ------------------------------
  # Observation model
  # ------------------------------
  for(k in 1:N.trials){
    y[k] ~ dnorm(t[k], y.precision)
  }
}"
    
  } else {
    
    jags.pseudo <- "model{

  # ------------------------------
  # Priors
  # ------------------------------
  
  # On average whale threshold
  mu ~ dunif(lowerbound, upperbound)
  
  # Between-whale variation (sd)
  phi ~ dunif(0, phi.upper)
  inv_phi2 <- pow(phi, -2)
  
  # Within-whale variation (sd)
  sigma ~ dunif(0, sigma.upper)
  inv_sigma2 <- pow(sigma, -2)
  
  # Priors on covariates
  
  # Index for competing species models
  # step(e) ...... 1 if e >= 0; 0 otherwise
  for (l in 1:N.modelspecies) {theta_1[l] <- step(l-theta)*step(theta-l)}

  # Beta species priors 
  for (i in 1:N.beta_3) {beta_3[i] ~ dnorm(0,(1/30^2))} 

  # To simplify the process model (between species)
  # create an index of competing species models (v).
  for(i in 1:N.species){
    for (j in 1:N.beta_3){
      v[j,i] <- beta_3[j] * I.species[j,i]
    }
  }
  
  # ind.matr is a matrix of indices
  # Column 1 - ID of competing model
  # Column 2 to 3 - Row numbers (min [col2] to max [col3]) from index corresponding to competing model
  
  for (i in 1:N.species){
    for (j in 1:(N.modelspecies-1)){ 
	  sumtheta[j,i]<-theta_1[j]*(sum(v[ind.matr[j,2]:ind.matr[j,3],i]))
    }
  }

  # ------------------------------
  # Process model
  # ------------------------------

  # Between species
  for(i in 1:N.species){
    mu_i[i] <- mu + sum(sumtheta[1:(N.modelspecies-1),i]) +
	             theta_1[N.modelspecies]*(beta_3[N.beta_3])	
  }

  # Between whales within species
  for(j in 1:N.whales){
    mu_ij[j] <- mu_i[species.ID[j]]
    t_ij[j] ~ dnorm(mu_ij[j], inv_phi2) T(lowerbound, upperbound)   
  }

  # Between trials within whale
  for(k in 1:N.trials){
  
  # Expected threshold
  mu_ijk[k] <- t_ij[whale.ID[k]] #/#

  # Realized threshold
  t[k] ~ dnorm(mu_ijk[k], inv_sigma2) T(lowerbound, upperbound)    
  }

  # ------------------------------
  # Observation model
  # ------------------------------
  
  for(k in 1:N.trials){
    y[k] ~ dnorm(t[k], y.precision)
  }
}"
    
  }
  
  if(include.covariates){
    
    jags.pseudo <- gsub(
      pattern = "# Priors on covariates",
      replacement = paste0(names(I.covariates),
                           " ~ dnorm(", paste0(names(I.covariates), ".prior.mean"),
                           ", 1/(", paste0(names(I.covariates), ".prior.sd"), "^2))",
                           collapse = "\n "), x = jags.pseudo)
    
    # k.index <- purrr::map_chr(.x = fL, .f =~ifelse(.x$nL == 0, "[k]", "[k,]"))
    
    jags.pseudo <- gsub(
      pattern = "#/#",
      replacement = paste0("+ ", paste0(names(I.covariates), " * I.", 
                                        names(I.covariates), 
                                        "[k, ]",
                                        collapse = " + ")), x = jags.pseudo)
    
  }
  
  if(sum(I.censored) > 0)
    jags.pseudo <- gsub(pattern = "}\n}",
                        replacement = " I.censored[k] ~ dinterval(t[k], Rc[k]) \n }\n}",
                        x = jags.pseudo)
  
  
  cat("\n--------------------------------------------------\n")
  cat("PSEUDO-PRIORS\n")
  cat("--------------------------------------------------\n")
  
  # Loop over models
  for (i in 1:N.modelspecies){
    
    cat("\n Candidate model (", i," out of ", N.modelspecies, "): ", species.Groups[i], "\n", sep = "")
    
    #'-------------------------------------------------
    # Write the data
    #'-------------------------------------------------
    
    pseudo.data <- list(N.trials = N.trials,
                        N.whales = N.whales,
                        N.species = N.species,
                        lowerbound = lowerbound,
                        upperbound = upperbound,
                        phi.upper = phi.upper,
                        sigma.upper = sigma.upper,
                        theta = theta.n[i],
                        N.modelspecies = N.modelspecies,
                        whale.ID = whale.ID,
                        y = y,
                        y.precision = y.precision)
    
    if(random.effects) {
      pseudo.data <- append(pseudo.data, 
                            list(epsilon.upper = epsilon.upper,
                                 species_ID = species_id)) 
    } else {
      pseudo.data <- append(pseudo.data,
                            list(N.beta_3 = N.beta_3,
                                 I.species = I.species,
                                 ind.matr = ind.matr,
                                 species.ID = species.ID)) }
    
    if(include.covariates){
      for(nc in names(I.covariates)){
        pseudo.data[[paste0(nc, ".prior.mean")]] <- 0
        pseudo.data[[paste0(nc, ".prior.sd")]] <- 30
        pseudo.data[[paste0("I.", nc)]] <- I.covariates[[nc]]}}
    
    if(sum(I.censored) > 0) {
      pseudo.data$I.censored <- I.censored
      pseudo.data$Rc <- Rc}
    
    #'-------------------------------------------------
    # Set up initial values
    #'-------------------------------------------------
    
    # beta_3: species
    
    pseudo.inits <- list(mu = 150, phi = 5, sigma = 5)
    
    if(include.covariates){
      for(nc in names(I.covariates)) pseudo.inits[[nc]] <- 5 }
    
    if(random.effects) { pseudo.inits <- append(pseudo.inits, list(epsilon = 5)) 
    } else { pseudo.inits <- append(pseudo.inits, list(beta_3 = rep(0, N.beta_3))) }
    
    # if(sum(I.censored) > 0){
    #   t.inits <- numeric(N.trials)
    #   t.inits[I.censored == 0] <- runif(n = sum(I.censored == 0), min = lowerbound, upperbound)
    #   t.inits[I.censored == 1] <- runif(n = sum(I.censored == 1), min = Rc[I.censored == 1], upperbound)
    #   pseudo.inits$t <- t.inits}
    
    pseudo.inits$t <- runif(n = N.trials, min = Rc, max = upperbound)
    
    #'-------------------------------------------------
    # Run the model in JAGS
    #'-------------------------------------------------
    
    m.pseudo <- rjags::jags.model(file = textConnection(jags.pseudo), 
                                  n.chains = n.chains, 
                                  data = pseudo.data, 
                                  inits = pseudo.inits, 
                                  quiet = TRUE)
    
    #'-------------------------------------------------
    # Burn-in
    #'-------------------------------------------------
    
    rjags:::update.jags(object = m.pseudo, n.iter = burnin)
    
    #'-------------------------------------------------
    # Draw samples from the posterior
    #'-------------------------------------------------
    
    if(random.effects) var.monitor <- paste0("mu_i[", i, "," , paste0(1:N.species),"]") else 
      var.monitor <- paste0("beta_3[", paste0(index[, 2][index[, 1] == i]), "]")
    
    pseudo.beta3 <- rjags::coda.samples(
      model = m.pseudo,
      variable.names = var.monitor,
      n.iter = pseudo.n)
    
    if(random.effects){
      
      if (i < N.modelspecies) results.species[1:N.species, 1:2, i] <- summary(pseudo.beta3)$statistics[,1:2]
      if (i == N.modelspecies) results.species[1:N.species, 1:2, i] <- summary(pseudo.beta3)$statistics[1:2]  
      
    } else  {
      
      if (i < N.modelspecies) results.species[index[, 1] == i, 1:2] <- summary(pseudo.beta3)$statistics[, 1:2]
      if (i == N.modelspecies) results.species[index[, 1] == i, 1:2] <- summary(pseudo.beta3)$statistics[1:2]
    }
    
  }
  
  # || Main model ----
  
  cat("\n--------------------------------------------------\n")
  cat("MAIN MODEL\n")
  cat("--------------------------------------------------\n\n")
  
  if(random.effects){
    
    jags.main <- "model{

  # ------------------------------
  # Priors
  # ------------------------------
  
  # Average whale threshold
  mu ~ dunif(lowerbound, upperbound)
  
  # Between-whale variation (sd)
  phi ~ dunif(0, phi.upper)
  inv_phi2 <- pow(phi, -2)
  
  # Within-whale variation (sd) 
  sigma ~ dunif(0, sigma.upper)
  inv_sigma2 <- pow(sigma, -2)
  
  # Between-species variation (sd)
  epsilon ~ dunif(0, epsilon.upper) 
  inv_epsilon2 <- pow(epsilon, -2)
  
  # Priors on covariates
  
  # Index for competing species models
  theta ~ dcat(pi.species)
  for (l in 1:N.modelspecies) {theta_1[l] <- step(l-theta)*step(theta-l)}

  # ------------------------------
  # Process model
  # ------------------------------

  # Between species
  for(i in 1:N.species){
    for (kk in 1:N.modelspecies){
      mu_i_prior[kk, i] ~ dnorm(mu, inv_epsilon2) T(lowerbound, upperbound)    
      mu_i_pseudoprior[kk, i] ~ dnorm(results.species[i, 1, kk], 1/results.species[i, 2, kk]^2) T(lowerbound, upperbound)   
      mu_i[kk, i] <- theta_1[kk] * mu_i_prior[kk, i] + (1 - theta_1[kk]) * mu_i_pseudoprior[kk, i] 
    }
  }

  for(j in 1:N.whales){
	for (kk in 1:N.modelspecies){
	  v[kk, j] <- mu_i[kk, species_ID[kk, j]] * theta_1[kk]
	}
    mu_ij[j] <- sum(v[1:N.modelspecies, j])

    # Between whales within species
    t_ij[j] ~ dnorm(mu_ij[j], inv_phi2) T(lowerbound, upperbound)
  }

  # Between trials within whale
  for(k in 1:N.trials){
  
  # Expected threshold
    mu_ijk[k] <- t_ij[whale.ID[k]] #/#

  # Realized threshold
    t[k] ~ dnorm(mu_ijk[k], inv_sigma2) T(lowerbound, upperbound)
  }

  # ------------------------------
  # Observation model
  # ------------------------------
  for(k in 1:N.trials){
    y[k] ~ dnorm(t[k], y.precision)
  }

}"
    
  } else {
    
    jags.main <- "model{

  # ------------------------------
  # Priors
  # ------------------------------
  
  # On average whale threshold
  mu ~ dunif(lowerbound, upperbound)
  
  # Between-whale variation (sd)
  phi ~ dunif(0, phi.upper)
  inv_phi2 <- pow(phi, -2) 
  
  # Within-whale variation (sd)
  sigma ~ dunif(0, sigma.upper)
  inv_sigma2 <- pow(sigma, -2)
  
  # Priors on covariates
  
  # Index for competing species models 
  theta ~ dcat(pi.species) 
  for (l in 1:N.modelspecies) {theta_1[l] <- step(l-theta)*step(theta-l)}

  # Beta species priors and pseudopriors
  for (j in 1:N.beta_3) {
    beta_prior_3[j] ~ dnorm(0, 1/30^2)
    beta_pseudoprior_3[j] ~ dnorm(results.species[j,1], 1/results.species[j,2]^2)
    beta_3[j] <- theta_1[index[j]]*beta_prior_3[j] + (1-theta_1[index[j]])*beta_pseudoprior_3[j]
  }

  # To simplify the process model (between species)
  for (i in 1:N.species){
    for (j in 1:N.beta_3){
      v[j,i] <- beta_3[j]*I.species[j,i]
    }
  }
  
  for (i in 1:N.species){
    for (j in 1:(N.modelspecies-1)){ 
	       sumtheta[j,i] <- theta_1[j]*(sum(v[ind.matr[j,2]:ind.matr[j,3],i]))
    }
  }

  # ------------------------------
  # Process model
  # ------------------------------

  # Between species
  for(i in 1:N.species){
    mu_i[i] <- mu + sum(sumtheta[1:(N.modelspecies-1),i]) +
	             theta_1[N.modelspecies]*(beta_3[N.beta_3])	
  }

  # Between whales within species
  for(j in 1:N.whales){
    mu_ij[j] <- mu_i[species.ID[j]]
    t_ij[j] ~ dnorm(mu_ij[j], inv_phi2) T(lowerbound, upperbound) 
  }

  # Between trials within whale
  for(k in 1:N.trials){
    
  # Expected threshold
  mu_ijk[k] <- t_ij[whale.ID[k]] #/#

  # Realized threshold
    t[k] ~ dnorm(mu_ijk[k], inv_sigma2) T(lowerbound, upperbound)    
  }

  # ------------------------------
  # Observation model
  # ------------------------------
  
  for(k in 1:N.trials){
    y[k] ~ dnorm(t[k], y.precision)
  }
}"
    
  }
  
  
  if(include.covariates){
    
    jags.main <- gsub(
      pattern = "# Priors on covariates",
      replacement = paste0(names(I.covariates),
                           " ~ dnorm(", paste0(names(I.covariates), ".prior.mean"),
                           ", 1/(", paste0(names(I.covariates), ".prior.sd"), "^2))",
                           collapse = "\n "), x = jags.main)
    
    # k.index <- purrr::map_chr(.x = fL, .f =~ifelse(.x$nL == 0, "[k]", "[k,]"))
    
    jags.main <- gsub(
      pattern = "#/#",
      replacement = paste0("+ ", paste0(names(I.covariates), " * I.", 
                                        names(I.covariates), 
                                        "[k, ]", 
                                        collapse = " + ")), x = jags.main)
    
  }
  
  if(sum(I.censored) > 0)
    jags.main <- gsub(pattern = "}\n}",
                      replacement = " I.censored[k] ~ dinterval(t[k], Rc[k]) \n }\n}",
                      x = jags.main)
  
  #'-------------------------------------------------
  # Set up initial values for the model
  #'-------------------------------------------------
  
  main.inits <- list(mu = 190, 
                     phi = 5, 
                     sigma = 5)
  
  if(include.covariates){
    for(nc in names(I.covariates)) main.inits[[nc]] <- 5}
  
  if(random.effects) { main.inits <- append(main.inits, list(epsilon = 5)) 
  } else { 
    main.inits <- append(main.inits, list(beta_prior_3 = rep(0, N.beta_3),
                                          beta_pseudoprior_3 = results.species[, 1])) }
  
  main.inits$t <- runif(n = N.trials, min = Rc, max = upperbound)
  
  # if(sum(I.censored) > 0){
  #   t.inits <- numeric(N.trials)
  #   t.inits[I.censored == 0] <- runif(n = sum(I.censored == 0), min = lowerbound, Rc[I.censored == 0])
  #   t.inits[I.censored == 1] <- runif(n = sum(I.censored == 1), min = Rc[I.censored == 1], upperbound)
  #   main.inits$t <- t.inits
  # }
  
  #'-------------------------------------------------
  # Data for the model
  #'-------------------------------------------------
  
  main.data <- list(N.trials = N.trials,
                    N.whales = N.whales,
                    N.species = N.species,
                    N.modelspecies = N.modelspecies,
                    lowerbound = lowerbound,
                    upperbound = upperbound,
                    results.species = results.species,
                    phi.upper = phi.upper,
                    sigma.upper = sigma.upper,
                    pi.species = pi.species,
                    whale.ID = whale.ID,
                    y = y,
                    y.precision = y.precision)
  
  if(random.effects) {
    main.data <- append(main.data, 
                        list(epsilon.upper = epsilon.upper,
                             species_ID = species_id))
  } else { 
    main.data <- append(main.data, 
                        list(N.beta_3 = N.beta_3,
                             index = index[, 1],
                             I.species = I.species,
                             ind.matr = ind.matr,
                             species.ID = species.ID)) }
  
  if(include.covariates){
    for(nc in names(I.covariates)){
      main.data[[paste0(nc, ".prior.mean")]] <- 0
      main.data[[paste0(nc, ".prior.sd")]] <- 30
      main.data[[paste0("I.", nc)]] <- I.covariates[[nc]]}}
  
  if(sum(I.censored) > 0) {
    main.data$I.censored <- I.censored
    main.data$Rc <- Rc}
  
  #'-------------------------------------------------
  # Run the model in JAGS
  #'-------------------------------------------------
  
  m.main <- rjags::jags.model(file = textConnection(jags.main), 
                              data = main.data, 
                              n.chains = n.chains,
                              inits = main.inits, 
                              quiet = TRUE)
  
  #'-------------------------------------------------
  # Burn-in
  #'-------------------------------------------------
  
  rjags:::update.jags(object = m.main, n.iter = burnin)
  
  #'-------------------------------------------------
  # Sample from the posterior
  #'-------------------------------------------------
  
  if(random.effects){
    cov.monitor <- c("theta", "theta_1", "mu_i", "phi", "sigma", "epsilon")
  } else {
    cov.monitor <- c("theta", "mu_i", "phi", "sigma")
  }
  
  if(include.covariates)  cov.monitor <- c(cov.monitor, names(I.covariates))
  
  samples.ms <- rjags::coda.samples(model = m.main,
                                    variable.names = cov.monitor,
                                    n.iter = mcmc.n,
                                    thin = thin)
  
  ## || Results ----
  
  exec.time <- tictoc::toc(quiet = TRUE)
  run_time <- hms::as_hms(as.numeric(round(exec.time$toc - exec.time$tic, 0)))
  
  run_time <- tibble::enframe(run_time) %>% 
    dplyr::rename(Chains = name, Time = value) %>% 
    dplyr::mutate(Chains = paste0(range(seq_len(n.chains)), collapse = ":"))
  
  mcmc.params$run_time <-  run_time
  
  for(nc in seq_len(n.chains)){
    colnames(samples.ms[[nc]])[which(grepl(pattern = "mu", x = colnames(samples.ms[[nc]])))] <- 
      paste0("mu.", species.names)
  }
       
  res <- list(trace = samples.ms,
              mcmc = mcmc.params,
              dat = dat,
              config = config,
              species.Groups = species.Groups,
              type = ifelse(random.effects, "random", "fixed"))
  
  class(res) <- c("gvs", class(res))
  return(res) 
  
  
}