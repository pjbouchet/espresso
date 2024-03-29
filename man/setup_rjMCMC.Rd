% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_rjMCMC.R
\name{setup_rjMCMC}
\alias{setup_rjMCMC}
\title{Set up the rjMCMC sampler}
\usage{
setup_rjMCMC(
  rj.input,
  n.burn,
  n.iter,
  n.chains = 3,
  p.split,
  p.merge,
  do.update = FALSE,
  start.values = NULL
)
}
\arguments{
\item{rj.input}{Input dataset. Must be an object of class \code{rjconfig}.}

\item{n.burn}{Number of MCMC iterations to treat as burn-in.}

\item{n.iter}{Number of posterior samples.}

\item{n.chains}{Number of MCMC chains. Defaults to 3.}

\item{p.split}{Probability of performing a group split when \code{model.select = TRUE}.}

\item{p.merge}{Probability of performing a group merge when \code{model.select = TRUE}.}

\item{do.update}{Logical. Whether to update an existing sampler or set up a new one.}

\item{start.values}{Starting values to use when updating an existing sampler and \code{do.update = TRUE}.}
}
\description{
Set up a list object to hold the data from the MCMC sampler, and generate starting values for all model parameters. This function is called internally by \code{\link{run_rjMCMC}}.
}
\details{
In addition to split/merge moves, three other types of MCMC samplers are implemented in \code{espresso} to facilitate convergence and avoid getting stuck in local maxima: two data-driven samplers (type I and type II), in which proposals are informed by the “cues” present in the original data, and one independent sampler, in which proposals are drawn at random from a Uniform distribution bounded by \code{range.dB} (see \code{\link{read_data}} or \code{\link{simulate_data}}), with no dependence on current values.
}
\seealso{
\code{\link{run_rjMCMC}}
}
\author{
Phil J. Bouchet
}
\keyword{-----------------------------------------------}
\keyword{brs}
\keyword{dose-response}
\keyword{rjmcmc}
