#' Initial configuration
#'
#' Define several parameters required to configure the rjMCMC sampler.
#' 
#' The function performs three actions:
#' \describe{
#'   \item{Variance}{Returns empirical estimates of between-whale (φ) and within-whale between-exposure (σ) variation, which are needed to generate starting values for the MCMC chains}
#'   \item{Proposals}{Defines the means and standard deviations of relevant proposal distributions and priors}
#'   \item{Clustering}{Performs a cluster analysis using \code{n.rep} bootstrap replicates of the input data, in order to parameterise between-model jumps. The analysis is run using \code{\link[mclust]{Mclust}} on both (1) species means and (2) raw observations. These respectively provide estimates of the probability distributions of (1) unique groupings, and (2) numbers of groups, where the latter ranges from \code{n = 1} when all species belong to the same group, to \code{n = n.species} when each species is treated individually.}
#' }
#' @export 
#' @param dat Input BRS data. Must be an object of class \code{brsdata}.
#' @param model.select Logical. If \code{TRUE}, between-model jumps will be allowed in the rjMCMC sampler. This argument can be set to \code{FALSE} to impose a species grouping and constrain parameter estimation accordingly.
#' @param covariate.select Logical. Set to \code{TRUE} to allow contextual covariates to drop out of the model.
#' @param proposal.mh Named list specifying the standard deviations of the proposal distributions used in the Metropolis-Hastings sampler.
#' @param proposal.rj Named list specifying the standard deviations of the proposal distributions used in the reversible jump sampler. Must contain two elements: \code{mu} for response thresholds, and \code{cov} for covariates.
#' @param prior.covariates Mean and standard deviation of the Normal prior placed on covariates. At present, the same prior is applied to all covariates.
#' @param p.split Probability of choosing to split a group of species when initiating a split-merge move. This parameter is constrained to be \code{1} when all species are in a single group. \code{p.split} and \code{p.merge} must sum to \code{1}. 
#' @param p.merge Probability of choosing to merge two groups of species when initiating a split-merge move. This parameter is constrained to be \code{1} when all species are in their own groups. \code{p.split} and \code{p.merge} must sum to \code{1}.
#' @param move.ratio Relative proportion of calls to data-driven (type I and II) and independence samplers.
#' @param m Integer. Frequency (every \code{m} iterations) at which data-driven (type I and II) and independence samplers are triggered.
#' @param bootstrap Logical, defaults to TRUE. Whether to perform bootstrap clustering. As this step can be time-consuming, setting this argument to \code{FALSE} increase efficiency when reconfiguring the sampler. 
#' @param n.rep Number of replicate bootstrap datasets used for clustering (see \code{Details}).
#' 
#' @return A list object of class \code{rjconfig} identical to the input \code{brsdata} object, with an additional \code{config} element composed of:
#' 
#' \tabular{ll}{
#'   \code{var} \tab Empirical estimates of φ and σ \cr
#'   \code{prop} \tab Standard deviations of MCMC proposal distributions \cr
#'   \code{clust} \tab Outputs from the cluster analysis \cr
#'   \code{boot} \tab Species group assignments for each bootstrap dataset \cr
#'   \code{mlist} \tab Unique species group assignments \cr
#'   \code{prior} \tab Mean and standard deviation of the Normal priors placed on covariates \cr
#'   \code{...} \tab Other elements pulled from the \code{brsdata} object \cr
#'  }
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Import the example data, excluding species 
#' # with less than 5 observations
#' mydat <- read_data(file = NULL, min.N = 5) 
#' summary(mydat)
#' 
#' # Configure the rjMCMC sampler
#' mydat.config <- configure_rjMCMC(dat = mydat,
#'                                  model.select = TRUE,
#'                                  covariate.select = FALSE,
#'                                  proposal.mh = list(t.ij = 10, mu.i = 10, 
#'                                                     mu = 7, phi = 10, sigma = 10),
#'                                  proposal.rj = list(dd = 20, cov = 7),
#'                                  prior.covariates = c(0, 30),
#'                                  n.rep = 100)
#' summary(mydat.config)
#' }
#' @keywords brs rjmcmc dose-response

configure_rjMCMC <- function(dat, 
                             model.select = TRUE,
                             covariate.select = FALSE,
                             proposal.mh = list(t.ij = 10, mu.i = 10, mu = 7, phi = 10, sigma = 10),
                             proposal.rj = list(dd = 1, cov = 7),
                             prior.covariates = c(0, 30),
                             p.split = 0.5,
                             p.merge = 0.5,
                             move.ratio = list(dd1 = 3, dd2 = 1, random = 1),
                             m = 100,
                             bootstrap = TRUE,
                             n.rep = 1000){
  
  #' -----------------------------------------------
  # Perform function checks
  #' -----------------------------------------------
  
  if(!"brsdata" %in% class(dat)){
    if(!"rjconfig" %in% class(dat)) stop("Input data must be of class <brsdata> or <rjconfig>.")}
  
  if(!bootstrap & !"rjconfig" %in% class(dat)) stop("<rjconfig> object required when bootstrap = FALSE.")
  
  if(!length(move.ratio) == 3) stop("Length mismatch for <move.ratio>.")
  if(!sum(p.split, p.merge) == 1) stop("Move probabilities do not sum to 1.")
  if(sum(c(p.split, p.merge) < 0) > 0) stop("Move probabilities cannot be negative.")
  
  if(dat$covariates$n == 0) covariate.select <- FALSE
  if(dat$species$n == 1) model.select <- FALSE
  
  #' -----------------------------------------------
  # Estimate variance parameters from the data
  #' -----------------------------------------------
  
  # Estimate the between-whale SD (phi) and the within-whale between-exposure SD (sigma) from the data.
  # This is useful for generating appropriate starting values for the MCMC sampler (see below).
  
  dat.sigma <- mean(unique(sapply(X = dat$whales$id, FUN = function(n) 
    sd(dat$obs$y_ij[dat$whales$id == n], na.rm = TRUE))), na.rm = TRUE)
  
  dat.phi <- unique(sapply(X = dat$whales$id, FUN = function(n) 
    mean(dat$obs$y_ij[dat$whales$id == n], na.rm = TRUE)))
  
  dat.phi <- mean(sapply(X = dat$species$id, FUN = function(n) 
    sd(dat.phi[dat$species$id == n], na.rm = TRUE)), na.rm = TRUE)
  
  #' -----------------------------------------------
  # Define proposal distributions (RJ and MH steps) + priors
  #' -----------------------------------------------
  
  if(dat$covariates$n > 0){
    
    for (j in dat$covariates$names){
      proposal.mh[[j]] <- ifelse(j == "range", 3, 10)
    }
    
    # Priors on covariate effects
    Normal.priors <- lapply(X = 1:dat$covariates$n, FUN = function(x) prior.covariates)
    names(Normal.priors) <- dat$covariates$names
   
  }
  
  #' -----------------------------------------------
  # Clustering
  #' -----------------------------------------------
  
  # Perform a model-based cluster analysis based on the simulated / real data.
  # The analysis is performed on bootstrap replicates of the data, and both on
  # species means and on observations. The former returns probabilities for each candidate model,
  # whilst the former assigns probabilities to specific numbers of clusters (species groups).
  # This is useful for parameterising alternative types of model jumps.
  
  if(bootstrap){
    
  if(model.select){
    
    # NA values cause numerical issues so need to be dealt with,
    # Cannot simply be removed as that would be problematic for species that only have censored data
    
    boot.dat <- tibble::tibble(y = dat$obs$y_ij, 
                               species = dat$species$trials, 
                               censored = dat$obs$censored,
                               Lc = dat$obs$Lc,
                               Rc = dat$obs$Rc) 

    boot.dat <- boot.dat %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(y = ifelse(is.na(y), 
                               ifelse(censored == -1,
                                      runif(n = 1, min = dat$param$bounds["y.ij", 1], max = Lc),
                                      runif(n = 1, min = Rc, max = dat$param$bounds["y.ij", 2])), y)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-Rc, -Lc, -censored)

    # Bootstrap the data by species
    # Taken from https://humanitarian-user-group.github.io/post/tidysampling/
    # Contrary to the bootstraps function from the <rsamples> package (original implementation),
    # this approach ensures that all species are included in all resamples.
    
    cat("Generating Bootstrap resamples ...\n")
    
    boot.tbl <- purrr::map(.x = seq_len(n.rep),
                           .f = ~{
                             boot.dat %>% 
                             # Partition the data into strata (by species)
                             dplyr::group_nest(species) %>% 
                            # For each stratum
                             dplyr::rowwise() %>%
                             dplyr::mutate(
                               # calculate the required sample size
                               Ns = nrow(data),
                               # then draw the sample
                               sample = list(dplyr::sample_n(data, size = Ns, replace = TRUE))) %>% 
                             dplyr::select(-c(data, Ns)) %>% 
                             tidyr::unnest(sample)}) 
    
    # Original implementation – not reliable
    # boot.tbl <- rsample::bootstraps(data = boot.dat, times = n.rep, strata = "species", apparent = TRUE) %>% 
    #   purrr::map(.x = .$splits, .f = ~rsample::analysis(.x))
    
    # Compute average per species
    boot.avg <- purrr::map(.x = boot.tbl, 
                           .f = ~{
                             .x %>% 
                               dplyr::group_by(species) %>% 
                               dplyr::summarise(y = mean(y, na.rm = TRUE), .groups = "keep") %>% 
                               dplyr::ungroup()})
    
    cat("Performing cluster analysis ...\n")
    
    # Model-based clustering (equal variance)
    datClust <- purrr::map(.x = boot.tbl,
                           .f = ~mclust::Mclust(data = as.numeric(na.omit(.x$y)), 
                                                G = seq_len(dat$species$n), 
                                                modelNames = c("E"), verbose = FALSE))
    
    datClust.avg <- purrr::map(.x = boot.avg,
                               .f = ~{
                                 # If all species randomly have the same average, algorithm gets stuck
                                 if(length(unique(as.numeric(na.omit(.x$y)))) == 1) 
                                   G_mix <- 1 else G_mix <- seq_len(dat$species$n)
                                   mclust::Mclust(data = as.numeric(na.omit(.x$y)), 
                                                  G = G_mix, modelNames = c("E"), verbose = FALSE)})
    
    # // Bootstrap moves 
    
    datGroups <- purrr::map(.x = 1:length(datClust.avg), 
                            .f = ~{
                              
                              # Use a Multinomial distribution to classify species into groups 
                              # based on the assignment probabilities returned by the clustering algorithm
                              if(!is.null(datClust.avg[[.x]])){
                                clustClassif <- sapply(X = 1:nrow(datClust.avg[[.x]]$z), 
                                                       FUN = function(x) stats::rmultinom(n = 1, size = 1, prob = datClust.avg[[.x]]$z[x,]))
                                
                                if("matrix" %in% class(clustClassif)) clustClassif <- 
                                    apply(X = clustClassif, MARGIN = 2, FUN = function(x) which(x == 1))
                                
                                # The above groups do not necessarily have the same labels as the ones used
                                # in the other RJ functions - so must relabel them.
                                # First identify which species are in which groups
                                class.list <- list()
                                for(u in unique(clustClassif)) class.list[[u]] <- which(clustClassif == u)
                                class.list <- purrr::compact(class.list)
                                
                                # Relabel and find matching group from master list
                                res.out <- format_group(input.list = class.list)
                                
                                res.out <- res.out$group %>% relabel(.)
                              } else {
                                res.out <- 1 
                              }})
    
    names(datGroups) <- purrr::map(.x = datGroups, 
                                   .f = ~ vec_to_model(input.vector = .x, sp.names = dat$species$names)) %>%
      do.call(c, .)
    
    p.ClustModels <- table(names(datGroups)) / n.rep
    p.ClustModels <- as.data.frame(p.ClustModels) %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(model = Var1, p = Freq) %>% 
      dplyr::mutate(model = as.character(model)) %>% 
      dplyr::arrange(-p) %>% 
      dplyr::mutate(p_scale = rescale_p(p))
    
    # // Cluster moves 
    
    n.Clust <- purrr::map_dbl(.x = datClust, .f = ~which.max(.x$BIC[seq_len(dat$species$n)]))
    p.Clust <- sapply(X = seq_len(dat$species$n), 
                      FUN = function(x) {tmp <- sum(n.Clust == x, na.rm = TRUE) / n.rep
                      names(tmp) <- x; tmp})
    p.Clust <- tibble::tibble(cluster = as.numeric(names(p.Clust)), p = p.Clust) %>% 
      dplyr::mutate(p_scale = rescale_p(p))
    
    mlist <- unique(datGroups)
    names(mlist) <- purrr::map(.x = mlist, 
                               .f = ~ vec_to_model(input.vector = .x, sp.names = dat$species$names)) %>%
      do.call(c, .)
    
  } else {
    
    p.ClustModels <- tibble::tibble(model = NA, p = 1)
    p.Clust <- tibble::tibble(cluster = NA, p = 1)
    if(dat$species$n == 1) mlist <- datGroups <- list(rep(1, dat$species$n)) else 
      mlist <- datGroups <- list(seq_len(dat$species$n)) 
    names(mlist) <- names(datGroups) <- vec_to_model(input.vector = unlist(mlist),
                                                     sp.names = dat$species$names)
    p.ClustModels$model <- names(mlist)
    p.Clust$cluster <- ifelse(is.null(dat$species$groups), dat$species$n, length(dat$species$groups))
    
  }
    
  }
  
  if(bootstrap){
    
    dat$config <- list(var = setNames(c(dat.sigma, dat.phi), c("sigma", "phi")),
                       prop = list(dd = proposal.rj$dd, mh = proposal.mh),
                       move = list(prob = setNames(c(p.split, p.merge), c("split", "merge")),
                                   ratio = move.ratio, freq = m),
                       model.select = model.select,
                       boot = datGroups,
                       mlist = mlist,
                       clust = list(p.ClustModels, p.Clust),
                       covariate.select = covariate.select)

  } else {
    
    dat$config$var <-  setNames(c(dat.sigma, dat.phi), c("sigma", "phi"))
    dat$config$prop <- list(dd = proposal.rj$dd, mh = proposal.mh)
    dat$config$move <- list(prob = setNames(c(p.split, p.merge), c("split", "merge")),
                                   ratio = move.ratio, freq = m)
    dat$config$model.select <-  model.select
    dat$config$covariate.select <- covariate.select
  } 
  
  if(dat$covariates$n > 0){
    dat$config$prop$cov <- proposal.rj$cov
    dat$config$prior <- Normal.priors
  }
  
  cat("Done!")
  class(dat) <- c("rjconfig", class(dat))
  return(dat)
  
}

