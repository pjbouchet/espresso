#' Bayesian inference for multi-species behavioural dose-response functions
#'
#' The \code{espresso} package provides methods for fitting and selecting among behavioural dose-response models using Bayesian reversible jump Markov Chain Monte Carlo (rjMCMC). 
#' 
#' The package was developed to explore similarities in patterns of responsiveness across a range of cetacean species exposed to military sonar, and forms an output from the Double Mocha project led by the University of St Andrews and Duke University. For more details, see \href{https://synergy.st-andrews.ac.uk/mocha/}{https://synergy.st-andrews.ac.uk/mocha/}.
#'
#' @author Phil J. Bouchet
#' @docType package
#' @name espresso
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "4.0.0")  utils::globalVariables(c("."))