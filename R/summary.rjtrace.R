#' Some diagnostics for fitted dose-response models
#'
#' 
#' Summary method for objects of class \code{rjtrace}, as returned by \code{\link{trace_rjMCMC}}. Produces a text-based summary of: (1) effective sample sizes, (2) acceptance rates, (3) model convergence, (4) posterior model probabilities, and (5) posterior inclusion probabilities (PIPs) for contextual covariates (where appripriate).
#'
#' @export
#' @param rj.obj Input rjMCMC object, as returned by \code{\link{trace_rjMCMC}}.
#' @param covariate.prob Logical. If \code{TRUE}, returns a summary of posterior inclusion probabilities (PIPs).
#' @param combine.chains Logical. If \code{TRUE}, outputs are returned for each individual MCMC chain.
#' @inheritParams summary.gvs
#' 
#' @return A detailed summary, printed to the R console.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.rjdata}}
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
#' summary(rj)
#' }
#' @keywords brs dose-response rjmcmc 

summary.rjtrace <- function(rj.obj, 
                            eff.n = TRUE,
                            accept.rate = TRUE,
                            convergence = TRUE,
                            gelman.rubin = 1.1,
                            model.ranks = TRUE,
                            n.top = 10, 
                            covariate.prob = TRUE,
                            combine.chains = FALSE){
  
  options(tibble.width = Inf) 
  options(pillar.neg = FALSE) 
  options(pillar.subtle = TRUE)
  options(pillar.sigfig = 4)
  
  cat("\n======================================================\n")
  cat("SUMMARY\n")
  cat("======================================================\n\n")
  
  cat("Iterations:", format(rj.obj$mcmc$iter.rge, big.mark = ","), "\n")
  cat("Thinning interval:", rj.obj$mcmc$thin, "\n")
  cat("Burn-in:", format(rj.obj$mcmc$n.burn, big.mark = ","), "\n")
  cat("Number of chains:", rj.obj$mcmc$n.chains, "\n")
  size.per.chain <- rj.obj$mcmc$n.iter / rj.obj$mcmc$thin
  cat("Sample size per chain:", format(size.per.chain, big.mark = ","), "\n")
  cat("Total sample size:", format(size.per.chain * rj.obj$mcmc$n.chains, big.mark = ","), "\n\n")
  cat("Run times:\n")
  print(rj.obj$mcmc$run_time)
  
  cat("\n--------------------")
  cat("\nMCMC\n")
  cat("--------------------\n")
  
  cat("Priors:\n\n")
  cat("μ: Uniform", paste0("[", paste0(rj.obj$dat$param$bounds["mu", ], collapse = "-"), "]"), "\n")
  cat("φ: Uniform", paste0("[", paste0(rj.obj$dat$param$bounds["phi", ], collapse = "-"), "]"), "\n")
  cat("σ: Uniform", paste0("[", paste0(rj.obj$dat$param$bounds["sigma", ], collapse = "-"), "]"), "\n")
  for(cc in rj.obj$dat$covariates$names){
  cat(paste0(cc, ":"), "Normal", paste0("(", paste0(rj.obj$config$prior[[cc]], collapse = ", "), ")"), "\n")
  }
  cat("\n")
  cat("p(split):", rj.obj$mcmc$move$prob[1], "\n")
  cat("p(merge):", rj.obj$mcmc$move$prob[2], "\n")
  cat("\n")
  print(rj.obj$mcmc$move$tab)
  
  if(eff.n){
    cat("\n--------------------")
    cat("\nEFFECTIVE SAMPLE SIZES\n")
    cat("--------------------\n")
    
    print(rj.obj$ess)
    # print(coda::effectiveSize(rj.obj$trace))
    # print(mcmcse::ess(rj.obj$trace))
  }
  
  if(accept.rate){
    
    cat("\n--------------------")
    cat("\nACCEPTANCE RATES\n")
    cat("--------------------\n")
    
    print(rj.obj$accept)
  }
  
  if(convergence){
    if(coda::nchain(rj.obj$trace) > 1){
      cat("\n--------------------")
      cat("\nCONVERGENCE ASSESSMENT\n")
      cat("--------------------\n")
      
      mctrace <- do.call(rbind, rj.obj$trace)
      
      # Excludes columns with a less than 10 unique value 
      # (e.g., when a model has been specified a priori, or when there are only 2 species,
      # and hence only two possible models)
      coda.out <- tryCatch(exp = {
        cvg <- coda::gelman.diag(x = rj.obj$trace,
                                 autoburnin = FALSE, 
                                 multivariate = TRUE)},
        error = function(e){ NA })
      
      if(all(is.na(coda.out))){
        cat("Convergence: Not assessed (insufficient MCMC samples)\n")
      } else {
        if(cvg$mpsrf < gelman.rubin){ 
          cat("Convergence: TRUE\n") 
        } else {
          cat("Convergence: FALSE\n")}
        print(cvg)
      }
    }
  }
  
  if(covariate.prob){
    
    cat("\n--------------------")
    cat("\nCOVARIATES\n")
    cat("--------------------\n")
    
    if(rj.obj$config$covariate.select){
      
      tb <- prob_covariates(obj = rj.obj, do.combine = combine.chains)
      
      if(combine.chains){
        
        cat("All chains (n = ", coda::nchain(rj.obj$trace), "):\n", sep = "")
        print(tb)
        cat("\n")
        
      } else {
        
        for(nc in seq_len(coda::nchain(rj.obj$trace))){
          cat("Chain ", nc, ":\n")
          print(tb[[nc]])
          cat("\n")
        }
      }
      
    } else {
      
      cat("Covariate selection: FALSE\n")
    }
    
  }
  
  if(model.ranks){
    
    cat("\n--------------------")
    cat("\nMODEL RANKINGS\n")
    cat("--------------------\n")
    
    if(rj.obj$config$model.select){
      
      res <- prob_models(input.obj = rj.obj, 
                         n.top = n.top,
                         mlist = rj.obj$mlist, 
                         select = rj.obj$config$model.select,
                         do.combine = combine.chains,
                         gvs = "gvs" %in% class(rj.obj))
      
      if(combine.chains){
        
        cat("All chains (n = ", coda::nchain(rj.obj$trace), ")\n", sep = "")
        if(!is.null(n.top)) cat("\nTop ", n.top, " models:\n")
        print(head(res$model$m_prob, ifelse(is.null(n.top), 9999, n.top)))
        cat("\n")
      
        ggres <- dplyr::left_join(x = res$model$m_prob, rj.obj$mlist[, c("model", "group")], by = "model")
        gg.cols <- pals::brewer.paired(max(unlist(ggres$group))) # Brewer paired
        m.matrix <- do.call(rbind, ggres$group)
        
        gg.matrix <- tidyr::expand_grid(x = seq_len(rj.obj$dat$species$n), y = seq_len(n.top)) %>% 
          dplyr::rowwise() %>% 
          dplyr::mutate(grouping = m.matrix[y, x]) %>% 
          dplyr::mutate(species = rj.obj$dat$species$names[y]) %>% 
          dplyr::ungroup() %>% 
          dplyr::mutate(y = n.top - y + 1)
        
        gg_model(dat = gg.matrix, rj.obj = rj.obj, colours = gg.cols, n.top = n.top, combine = TRUE)
        
      } else {
        
        for(nc in seq_len(coda::nchain(rj.obj$trace))){
          if(is.null(n.top)) cat("Chain ", nc, ":\n", sep = "") else
            cat("\nTop ", n.top, " models (chain ", nc, "):\n", sep = "")
          print(head(res[[nc]]$model$m_prob, ifelse(is.null(n.top), 9999, n.top)))
          cat("\n")
        } # End for loop
          
          gg.res <- purrr::map(.x = res, 
             .f = ~{rtib <- dplyr::left_join(x = .x$model$m_prob,  
                            rj.obj$mlist[, c("model", "group")], by = "model")
          do.call(rbind, rtib$group)})
          gg.cols <- pals::brewer.paired(max(unlist(gg.res)))
          
          m.matrix <- lapply(X = seq_len(n.top), 
                           FUN = function(x) do.call(rbind, purrr::map(.x = gg.res, .f = ~.x[x, ])))

          gg.matrix <- purrr::map(.x = m.matrix,
                                  .f = ~{
                                    tidyr::expand_grid(x = seq_len(rj.obj$dat$species$n),
                                                       y = seq_len(rj.obj$mcmc$n.chains)) %>% 
                                      dplyr::rowwise() %>% 
                                      dplyr::mutate(grouping = .x[y, x]) %>% 
                                      dplyr::mutate(species = rj.obj$dat$species$names[y]) %>% 
                                      dplyr::ungroup() %>% 
                                      dplyr::mutate(y = n.top - y + 1)})
          
          gg.tiles <- purrr::map(.x = seq_along(gg.matrix), 
                                 .f = ~gg_model(dat = gg.matrix[[.x]], colours = gg.cols, 
                                                n.top = n.top, rj.obj = rj.obj, combine = FALSE, no = .x))
          patchwork::wrap_plots(gg.tiles)

      }
    } else {
      
      cat("Model selection: FALSE\n")
    }
  }
}