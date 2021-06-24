#' Summarise dose-response data
#'
#' Summary method for objects of class \code{brsdata}. Print a text-based summary of the data in the R console.
#'
#' @export
#' @param dat.obj Input dataset. Must be an object of class \code{brsdata}, as returned by \code{\link{read_data}} or \code{\link{simulate_data}}.
#' @param print.config Logical. Whether to print a summary of MCMC parameters, when the data have been configured using \code{\link{configure_MCMC}}.
#' 
#' @return A text-based summary of dataset properties.
#'  
#' @author Phil J. Bouchet
#' @seealso \code{\link{read_data}} \code{\link{simulate_data}}
#' @examples
#' library(espresso)
#' 
#' # Simulate data for two species
#' mydat <- simulate_data(n.species = 2, 
#'                       n.whales = 16, 
#'                       max.trials = 3, 
#'                       covariates = list(exposed = c(0, 5), range = 10)
#'                       mu = c(125, 158), 
#'                       phi = 20, 
#'                       sigma = 20, 
#'                       Rc = c(210, 211), 
#'                       seed = 58697)
#'                       
#' summary(mydat)
#' @keywords brs rjmcmc dose-response

summary.brsdata <- function(dat.obj, print.config = TRUE){
  
  if(!"rjconfig" %in% class(dat.obj)) print.config <- FALSE
  
  cat("\n======================================================\n")
  cat("DATA SUMMARY\n")
  cat("======================================================\n")
  
  cat("\nSimulation:", dat.obj$param$sim)
  if(!dat.obj$param$sim) cat("\nData file:", dat.obj$param$data.file)
  
  cat("\n\n--------------------")
  cat("\nOBSERVATIONS\n")
  cat("--------------------\n")
  
  cat("Right-censoring:", any(dat.obj$obs$censored > 0), "\n")
  cat("Total:", sum(dat.obj$obs$censored), "\n")
  
  cat("\n--------------------")
  cat("\nSPECIES\n")
  cat("--------------------\n")
  
  cat("N = ", dat.obj$species$n, "\n")
  if(is.null(dat.obj$species$exclude) | length(dat.obj$species$exclude) == 0){
    cat("Excluded: None \n")
  } else {
    cat("Excluded:", paste0(unlist(dat.obj$species$exclude)), "\n")
  }
  
  if(is.null(dat.obj$species$groups)){
    cat("Species groupings: FALSE\n\n")
  } else {
    cat("Species groupings:\n\n")
    if(!is.null(dat.obj$species$abbrev)) {
      print(dat.obj$species$abbrev)
    } else { 
      groups.df <- tibble::tibble(name = names(dat.obj$species$groups), 
                                  species = unlist(purrr::map(.x = dat.obj$species$groups,
                                                              .f = ~paste0(.x, collapse = ","))))
      print(groups.df)}
    cat("\n")
  }
  
  # print(dat.obj$species$summary, n = ifelse(print.allrows, 9999, NULL))
  print(dat.obj$species$summary, n = 9999)
  
  cat("\n--------------------")
  cat("\nCOVARIATES\n")
  cat("--------------------\n")
  
  if(dat.obj$covariates$n == 0 ){
    
    cat("No covariates\n")
    
  } else {
    
    allcov <- vector(mode = "list", length = 2)
    
    for(cc in dat.obj$covariates$names){
      
      if(is.factor(dat.obj$covariates$df[, cc]) | all(dat.obj$covariates$df[, cc] %in% 0:1)){
        
        ct <- table(dat.obj$covariates$df[, cc])
        c.miss <- apply(dat.obj$covariates$df[, cc, drop = FALSE], 2, FUN = function(x) sum(is.na(x)))
        
        cat("-- ", cc, " --\n")
        cat("Levels:", length(ct), "\n")
        for(j in seq_along(ct)){
          cat("[", names(ct)[j], "]: n = ", ct[j], "\n", sep = "")
        }
        cat("[NA]: n = ", c.miss, "\n\n", sep = "")
        
      } else {
        
        c.min <- round(apply(dat.obj$covariates$df[, cc, drop = FALSE], 2, 
                             FUN = function(x) min(x, na.rm = TRUE)), 3)
        c.mean <- round(apply(dat.obj$covariates$df[, cc, drop = FALSE], 2, 
                              FUN = function(x) mean(x, na.rm = TRUE)), 3)
        c.max <- round(apply(dat.obj$covariates$df[, cc, drop = FALSE], 2, 
                             FUN = function(x) max(x, na.rm = TRUE)), 3)
        c.miss <- apply(dat.obj$covariates$df[, cc, drop = FALSE], 2, FUN = function(x) sum(is.na(x)))
        
        cat("-- ", cc, " --\n")
        cat("[min]: ", c.min, "\n[mean]: ", c.mean, "\n[max]: ", c.max, "\n", sep = "")
        cat("[NA]: n = ", c.miss, "\n\n", sep = "")
        cat("\n")
      }
    }
    
    
  }
  
  if("sonar" %in% dat.obj$covariates$names){
    
    cat("\n--------------------")
    cat("\nSONAR GROUPINGS\n")
    cat("--------------------\n")
    
    print(dat.obj$covariates$signal.types)
  }
  
  if(sum(is.na(dat.obj$covariates$df)) > 0) {
    cat("\n")
    warning("Some covariates contain missing values")}
  
  if(print.config){
    
    cat("\n--------------------")
    cat("\nMCMC\n")
    cat("--------------------\n\n")
    
    cat("Model selection:", dat.obj$config$model.select, "\n")
    cat("Covariate selection:", dat.obj$config$covariate.select, "\n\n")
    
    cat("Variance estimates:", paste0(round(dat.obj$config$var, 1), " (", names(dat.obj$config$var), ")", collapse = "; "), "\n")
    
    cat("\nProposal standard deviations:\n\n")
    prop.df <- tibble::enframe(dat.obj$config$prop$mh) %>% 
      tidyr::unnest(cols = c(value)) %>% 
      dplyr::rename(param = name, SD = value) %>% 
      dplyr::mutate(step = "MH") %>% 
      dplyr::bind_rows(., tibble::tibble(param = "data-driven", 
                                         SD = dat.obj$config$prop$dd, 
                                         step = "RJ"))
    print(prop.df)
    
    cat("\nClustering:\n\n")
    print(dat.obj$config$clust)
    
    
  }
  
}