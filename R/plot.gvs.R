#' Diagnostic plots
#'
#' Generate trace, density, and autocorrelation plots from a \code{gvs} object.
#'
#' @export
#' @param gvs.dat Gibbs Variable Selection object, as returned by \code{\link{gibbs}}.
#' @param param.name Parameter name(s). Defaults to \code{all}, which returns plots for all parameters in the model.
#' @param autocorr Logical. Whether to output chain autocorrelation plots.
#' @param individual Logical. If \code{TRUE}, separate density lines will be plotted for each chain. If \code{FALSE}, one density line will be plotted for all chains.
#' 
#' @details Adapted from Casey Youngflesh's function \code{\link[MCMCvis]{MCMCtrace}}.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{run_rjMCMC}} \code{\link{trace_rjMCMC}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Import the example data
#' 
#' mydat <- read_data(file = NULL) 
#' 
#' # Import a real dataset with the sonar and range covariates, 
#' # excluding sperm whales and any other species with a sample size
#' # smaller than two
#' 
#' mydat <- read_data(file = "path/to/my/data.csv", 
#'                   exclude.species = "Sperm whale",
#'                   min.N = 2) 
#' }
#' @keywords brs gvs dose-response 

plot.gvs <- function(gvs.dat, 
                     param.name = NULL, 
                     autocorr = FALSE,
                     individual = FALSE){
  
  if(!is.null(param.name)){
    
    param.name <- colnames(gvs.dat$trace[[1]])[grepl(pattern = param.name, x = colnames(gvs.dat$trace[[1]]))]
    
  }
  
  MCMC_trace(gvs.dat$trace, 
             iter = gvs.dat$mcmc$n.iter, 
             pdf = FALSE, 
             ind = individual,
             params = ifelse(is.null(param.name), "all", param.name))
  
  if(is.null(param.name)) bpars <- character() else bpars <- param.name
  if(autocorr) bayesplot::mcmc_acf(x = gvs.dat$trace, pars = bpars)
}