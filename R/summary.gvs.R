#' Some diagnostics for fitted dose-response models
#'
#' Summary method for objects of class \code{gvs}, as returned by \code{\link{gibbs}}. Produces a text-based summary of: (1) effective sample sizes, (2) acceptance rates, (3) model convergence, and (4) posterior model probabilities. Note that \code{\link{gibbs}} does not implement covariate selection. As a result, posterior inclusion probabilities (PIPs) are not returned here, contrary to \code{\link{summary.rjtrace}}.
#'
#' @export
#' @param gvs.obj Input trace object, as returned by \code{\link{gibbs}}.
#' @param eff.n Logical. If \code{TRUE}, returns estimates of the effective sample size for each model parameter.
#' @param accept.rate Logical. If \code{TRUE}, returns the acceptance rate (calculated after burn-in) for each model parameter.
#' @param convergence Logical. If \code{TRUE}, assesses convergence using the multivariate potential scale reduction factor (Gelman-Rubin statistic), as implemented in \code{\link[coda]{gelman.diag}}.
#' @param gelman.rubin Threshold for determining convergence based on the Gelman-Rubin statistic. Defaults to \code{1.1}.
#' @param model.ranks Logical. If \code{TRUE}, returns a summary of posterior model probabilities and associated model rankings.
#' @param n.top Number of top-ranking models to display when \code{model.ranks = TRUE}.
#' 
#' @return A detailed summary, printed to the R console.
#' @import mclust
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.brsdata}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Simulate data for two species
#' mydat <- simulate_data(n.species = 2,
#'                        n.whales = 16,
#'                        min.trials = 1,
#'                        max.trials = 3,
#'                        covariates = list(exposed = c(0, 5), range = 0.5),
#'                        mu = c(101, 158),
#'                        phi = 20,
#'                        sigma = 20,
#'                        Lc = c(60, 65),
#'                        Rc = c(210, 211),
#'                        seed = 58697)
#' summary(mydat)
#' 
#' # Model selection by GVS
#' gvs.model <- gibbs(dat = mydat,
#'              random.effects = FALSE,
#'              include.covariates = FALSE,
#'              mcmc.n = 1000,
#'              burnin = 500)
#' 
#' summary(gvs.model)
#' }
#' @keywords gvs dose-response

summary.gvs <- function(gvs.obj, 
                      eff.n = TRUE,
                      accept.rate = TRUE,
                      convergence = TRUE,
                      gelman.rubin = 1.1,
                      model.ranks = TRUE,
                      n.top = NULL 
                      ){
  
  cat("\n======================================================\n")
  cat("GVS SUMMARY\n")
  cat("======================================================\n\n")
  
  cat("Iterations:", format(gvs.obj$mcmc$iter.rge, big.mark = ","), "\n")
  cat("Thinning interval:", gvs.obj$mcmc$thin, "\n")
  cat("Burn-in:", format(gvs.obj$mcmc$n.burn, big.mark = ","), "\n")
  cat("Number of chains:", gvs.obj$mcmc$n.chains, "\n")
  cat("Sample size per chain:", format(gvs.obj$mcmc$n.iter, big.mark = ","), "\n")
  cat("Total sample size:", format(gvs.obj$mcmc$n.iter * gvs.obj$mcmc$n.chains, big.mark = ","), "\n\n")
  cat("Run times:\n")
  print(gvs.obj$mcmc$run_time)
  
  if(eff.n){
    cat("\n--------------------")
    cat("\nEFFECTIVE SAMPLE SIZES\n")
    cat("--------------------\n")
    
    # print(coda::effectiveSize(gvs.obj$trace))
    print(mcmcse::ess(gvs.obj$trace))
    
  }
  
  if(accept.rate){
    
    cat("\n--------------------")
    cat("\nACCEPTANCE RATES\n")
    cat("--------------------\n")
    
    print(round(1 - coda::rejectionRate(gvs.obj$trace), 4))
  }
  
  if(convergence){
    
    if(coda::nchain(gvs.obj$trace) > 1){
      cat("\n--------------------")
      cat("\nCONVERGENCE ASSESSMENT\n")
      cat("--------------------\n")
      
      mctrace <- do.call(rbind, gvs.obj$trace)
      
      # Excludes columns with a less than 10 unique value 
      # (e.g., when a model has been specified a priori, or when there are only 2 species,
      # and hence only two possible models)
      cvg <- coda::gelman.diag(x = gvs.obj$trace[, 
                                                 which(apply(mctrace, 2, FUN = function(df) length(unique(df))) > 10)],
                               autoburnin = FALSE, multivariate = TRUE)
      
      if(cvg$mpsrf < gelman.rubin){ 
        cat("Convergence: TRUE\n") 
      } else {
        cat("Convergence: FALSE\n")}
      
      print(cvg)
    }
  }
  
  if(model.ranks){
    
    cat("\n--------------------")
    cat("\nMODEL RANKINGS\n")
    cat("--------------------\n")
      
      #'-------------------------------------------------
      # Extract posterior samples
      #'-------------------------------------------------
      
      mcmc.values <- as.matrix(gvs.obj$trace) %>%
        tibble::as_tibble() %>%
        janitor::clean_names(.)
      
      #'-------------------------------------------------
      # Calculate posterior model probabilities
      #'-------------------------------------------------
      
      mselect.species <- table(mcmc.values$theta)/nrow(mcmc.values)
      model.ID.posterior <- as.numeric(names(sort(mselect.species, decreasing = TRUE)))
      names(mselect.species) <- gvs.obj$species.Groups[as.numeric(names(mselect.species))]
      mselect.species <- sort(mselect.species, decreasing = TRUE)
      m_prob <- tibble::tibble(as.data.frame(mselect.species))
      
      if(ncol(m_prob) == 1) m_prob <- m_prob %>%
        dplyr::rename(p = mselect.species) %>%
        dplyr::mutate(model = names(mselect.species)) %>%
        dplyr::select(model, p) else m_prob <- m_prob %>%
        dplyr::rename(model = Var1, p = Freq)
      
      #'-------------------------------------------------
      # Monte Carlo error
      #'-------------------------------------------------
      
      mtrace <- mcmc.values[["theta"]]
      
      # Compute Monte Carlo error
      mce <- sapply(X = unique(mtrace), FUN = function(x) as.numeric(mtrace == x))
      mc.error <- apply(X = mce, MARGIN = 2, FUN = mcmcse::mcse) %>% 
        purrr::map_dbl(.x = ., .f = "se")
      
      # Compute lower and upper interval bounds
      m_prob <- m_prob %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(lower = max(0, round(p - 1.98 * mc.error[dplyr::row_number()], 4)),
                      upper = min(1, round(p + 1.98 * mc.error[dplyr::row_number()], 4))) %>% 
        dplyr::ungroup()
      
      chain.rge <- range(seq_len(gvs.obj$mcmc$n.chains))
      
      if(is.null(n.top)) cat("(Chains ", paste0(chain.rge, collapse = ":"), "):\n", sep = "") else
        cat("\nTop ", n.top, " models (Chains ", paste0(chain.rge, collapse = ":"), "):\n", sep = "")
      print(head(m_prob, ifelse(is.null(n.top), 9999, n.top)))
      cat("\n")

  }
}