#' Some diagnostics for fitted dose-response models
#'
#' 
#' Summary method for objects of class \code{rjtrace}, as returned by \code{\link{trace_rjMCMC}}. Produces a text-based summary of: (1) effective sample sizes, (2) acceptance rates, (3) model convergence, (4) posterior model probabilities, and (5) posterior inclusion probabilities (PIPs) for contextual covariates (where appripriate).
#'
#' @export
#' @param rj.obj Input rjMCMC object, as returned by \code{\link{trace_rjMCMC}}.
#' @param covariate.prob Logical. If \code{TRUE}, returns a summary of posterior inclusion probabilities (PIPs).
#' @param ff.prob Logical. If \code{TRUE}, returns a summary of posterior probabilities for each functional form (monophasic vs. biphasic).
#' @param rmd Logical. This is used to create a different layout of plots when exporting results using \code{\link{create_report}}.
#' @param do.plot Logical. If \code{TRUE}, returns diagnostic plots in addition to text-based summaries.
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
                            ff.prob = TRUE,
                            rmd = FALSE,
                            do.plot = TRUE){
  
  options(tibble.width = Inf) 
  options(pillar.neg = FALSE) 
  options(pillar.subtle = TRUE)
  options(pillar.sigfig = 4)
  
  if(!rj.obj$config$model.select) model.ranks <- FALSE
  if(!rj.obj$config$covariate.select) covariate.prob <- FALSE
  if(!rj.obj$config$function.select) ff.prob <- FALSE
  
  cat("\n======================================================\n")
  cat("SUMMARY\n")
  cat("======================================================\n\n")
  
  cat("Total iterations:", format(rj.obj$mcmc$tot.iter, big.mark = ","), "\n")
  cat("Burn-in:", format(rj.obj$mcmc$n.burn, big.mark = ","), "\n")
  cat("Thinning interval:", rj.obj$mcmc$thin, "\n")
  cat("Iterations:", format(rj.obj$mcmc$iter.rge, big.mark = ","), "\n")
  cat("Number of chains:", rj.obj$mcmc$n.chains, "\n")
  cat("Sample size per chain:", format(rj.obj$mcmc$n.iter, big.mark = ","), "\n")
  cat("Total sample size:", format(rj.obj$mcmc$n.iter * rj.obj$mcmc$n.chains, big.mark = ","), "\n\n")
  cat("Run times:\n")
  print(rj.obj$mcmc$run_time)
  
  cat("\n--------------------")
  cat("\nMCMC\n")
  cat("--------------------\n")
  
  cat("Priors:\n\n")
  for(pp in 1:nrow(rj.obj$config$priors)){
    cat(rownames(rj.obj$config$priors)[pp], ": ", rj.obj$config$priors$prior[pp], " (", rj.obj$config$priors$lower_or_mean[pp], ", ", rj.obj$config$priors$upper_or_SD[pp], ")\n", sep = "")
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
    
    for(pp in c("mu", "phi", "sigma", "alpha", "nu.lower", "nu.upper", "tau", "omega", "psi", rj.obj$dat$covariates$names, "model", "phase")){
      cat("--", pp, "-- \n\n")
      my.ess <- rj.obj$ess %>%
        dplyr::filter(grepl(pattern = pp, x = parameter))
      if(pp %in% c("mu", "alpha", "nu.lower", "nu.upper")){
        my.ess <- my.ess %>% dplyr::mutate(parameter = rj.obj$dat$species$names) %>%
          tidyr::pivot_wider(., names_from = parameter, values_from = ESS)}
      print(my.ess)
      cat("\n")
    }

  }
  
  if(accept.rate){
    
    cat("\n--------------------")
    cat("\nACCEPTANCE RATES\n")
    cat("--------------------\n")
    
    AR <- rj.obj$accept %>% 
      dplyr::select(chain, param, AR, label) %>% 
      dplyr::arrange(label, param) %>% 
      dplyr::select(-label) %>% 
      tidyr::pivot_wider(names_from = param, values_from = AR)
    
    all.NA <- apply(X = AR, MARGIN = 2, FUN = function(x) all(is.na(x)))
    
    print(AR[, !all.NA])
  }
  
  if(convergence){
    
    if(coda::nchain(rj.obj$trace) > 1){
      
      cat("\n--------------------")
      cat("\nCONVERGENCE ASSESSMENT\n")
      cat("--------------------\n")
      
      # Gelman-Rubin diagnostic
      mctrace <- do.call(rbind, rj.obj$trace)
      
      # Excludes columns with a less than 10 unique value 
      # (e.g., when a model has been specified a priori, or when there are only 2 species,
      # and hence only two possible models)
      coda.out <- tryCatch(expr = {
        cvg <- coda::gelman.diag(x = rj.obj$trace[, sufficient.ess],
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
      
      if(do.plot){
      if(rj.obj$config$model.select){
      # Running means plots for model ID
      mcmcplots::rmeanplot(mcmcout = rj.obj$trace, parms = "model_ID", 
                           plot.title = "", auto.layout = FALSE,
                           col = gg_color_hue(coda::nchain(rj.obj$trace)))
      }}
      
    }
  }
  
  if(covariate.prob){
    
    cat("\n--------------------")
    cat("\nCOVARIATES\n")
    cat("--------------------\n\n")
    
    if(rj.obj$config$covariate.select){
      
      tb.combined <- prob_covariates(obj = rj.obj, do.combine = TRUE)
      tb <- prob_covariates(obj = rj.obj, do.combine = FALSE)

      cat("--- All chains (n = ", coda::nchain(rj.obj$trace), ") --- \n", sep = "")
      print(tb.combined)
      cat("\n")

      cat("--- Individual chains ---\n", sep = "")
      
      for(nc in seq_len(coda::nchain(rj.obj$trace))){
        cat("Chain ", nc, ":\n")
        print(tb[[nc]])
        cat("\n")
      }

    } else {
      
      cat("Covariate selection: FALSE\n")
    }
    
  }
  
  if(ff.prob){
    
    cat("\n--------------------")
    cat("\nFUNCTIONAL FORMS\n")
    cat("--------------------\n\n")
    
    if(rj.obj$config$function.select){
      
      tb.combined <- prob_form(obj = rj.obj, do.combine = TRUE) %>% 
        .[["est"]] %>% 
        dplyr::bind_cols(tibble::tibble(chain = "all"), .)
      tb <- prob_form(obj = rj.obj, do.combine = FALSE)
      
      tb.out <- purrr::map(.x = tb, "est")
      tb.out <- purrr::map(.x = seq_along(tb.out), .f = ~dplyr::mutate(.data = tb.out[[.x]], chain = .x)) %>% 
        do.call(rbind, .) %>% 
        dplyr::relocate(chain) %>% 
        dplyr::mutate(chain = as.character(chain)) %>% 
        dplyr::bind_rows(., tb.combined)
      
      print(tb.out)

    } else {
      cat("Functional form selection: FALSE\n")
    }
  }
  
  if(model.ranks){
    
    cat("\n--------------------")
    cat("\nMODEL RANKINGS\n")
    cat("--------------------\n\n")
    
    if(rj.obj$config$model.select){
      
      res <- prob_models(input.obj = rj.obj, 
                         n.top = n.top,
                         mlist = rj.obj$mlist, 
                         select = rj.obj$config$model.select,
                         do.combine = TRUE,
                         gvs = "gvs" %in% class(rj.obj))
      
      cat("--- All chains (n = ", coda::nchain(rj.obj$trace), ") --- \n", sep = "")
      if(!is.null(n.top)) cat("\nTop ", n.top, " models:\n", sep = "")
      print(head(res$model$m_prob, ifelse(is.null(n.top), 9999, n.top)))
      cat("\n")
      
      if(do.plot){
      ggres <- dplyr::left_join(x = res$model$m_prob, rj.obj$mlist[, c("model", "group")], by = "model")
      if(nrow(ggres) < n.top) n.top <- nrow(ggres)
      gg.cols <- if(max(unlist(ggres$group)) <= 2) pals::parula(max(unlist(ggres$group))) else pals::brewer.paired(max(unlist(ggres$group))) # Brewer paired
      m.matrix <- do.call(rbind, ggres$group)
      
      # Re-assign colours
      m.matrix <- t(apply(X = m.matrix, MARGIN = 1, FUN = function(x) 
       match(seq_along(unique(x)), unique(x))[x]))
      
      gg.matrix <- tidyr::expand_grid(x = seq_len(rj.obj$dat$species$n), y = seq_len(n.top)) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(grouping = m.matrix[y, x]) %>% 
        dplyr::mutate(species = rj.obj$dat$species$names[x]) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(y = n.top - y + 1)
      
      gg.combined <- gg_model(dat = gg.matrix,
                              post.probs = round(res$model$m_prob$p, 3),
                              rj.obj = rj.obj,
                              colours = gg.cols,
                              n.top = n.top,
                              combine = TRUE,
                              x.offset = 0.75,
                              x.margin = 1)
      }
      
      res <- prob_models(input.obj = rj.obj, 
                         n.top = n.top,
                         mlist = rj.obj$mlist, 
                         select = rj.obj$config$model.select,
                         do.combine = FALSE,
                         gvs = "gvs" %in% class(rj.obj))
      
      cat("--- Individual chains ---\n", sep = "")
      
      for(nc in seq_len(coda::nchain(rj.obj$trace))){
        if(is.null(n.top)){
          cat("Chain ", nc, ":\n", sep = "") 
        } else {
          cat("\nTop ", ifelse(nrow(res[[nc]]$model$m_prob) < n.top, nrow(res[[nc]]$model$m_prob), n.top), " models (chain ", nc, "):\n", sep = "")
        }
        print(head(res[[nc]]$model$m_prob, ifelse(is.null(n.top), 9999, n.top)))
        cat("\n")
      } # End for loop
      
      min.n <- min(purrr::map_dbl(.x = res, .f = ~nrow(.x$model$m_prob)))
      if(min.n < n.top) warning(paste0("Fewer than ", n.top, " models available in some MCMC chains. n.top set to ", min.n))
      
      if(do.plot){
      # Posterior probabilities for each ranked model across chains
      pprobs <- purrr::map(.x = seq_len(min(c(min.n, n.top))), 
                           .f = ~sapply(X = seq_len(coda::nchain(rj.obj$trace)),
                                        function(d) res[[d]]$model$m_prob$p[.x]))
      
      # Species groupings
      gg.res <- purrr::map(.x = res, 
                           .f = ~{
                             rtib <- dplyr::left_join(x = .x$model$m_prob[1:min(c(min.n, n.top)), ],  
                                                      rj.obj$mlist[, c("model", "group")], by = "model")
                             do.call(rbind, rtib$group)})
      
      gg.cols <- if(max(unlist(gg.res)) <=2) pals::parula(max(unlist(gg.res))) else pals::brewer.paired(max(unlist(gg.res)))
      
      m.matrix <- lapply(X = seq_len(min(c(min.n, n.top))), 
                         FUN = function(x) do.call(rbind, purrr::map(.x = gg.res, .f = ~.x[x, ])))
      
      # Re-assign colours
      m.matrix <- purrr::map(.x = m.matrix, .f = ~t(apply(X = .x, MARGIN = 1, FUN = function(x) 
        match(seq_along(unique(x)), unique(x))[x])))
      
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
                             .f = ~gg_model(dat = gg.matrix[[.x]], 
                                            post.probs = round(pprobs[[.x]], 3),
                                            colours = gg.cols, 
                                            n.top = min(c(min.n, n.top)), 
                                            rj.obj = rj.obj, 
                                            combine = FALSE, 
                                            p.size = 4,
                                            x.offset = 1,
                                            x.margin = 2,
                                            no = .x))
      
      
      if(!rmd){
        print(gg.combined)
        patchwork::wrap_plots(gg.tiles)
      } else {
        print(gg.combined)
        purrr::walk(.x = gg.tiles, .f = ~print(.x))
        
      }
      }
    } else {
      
      cat("Model selection: FALSE\n")
    }
  }

}
