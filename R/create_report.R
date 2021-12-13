#' Generate report
#'
#' Convenience function for creating a summary of results using R Markdown
#'
#' @export
#' @param outdir Output directory. Defaults to the current working directory.
#' @param filename Output file name (without the file extension). Defaults to \code{espresso_report.html}.
#' @param plot.height Height of tile plots used for comparing species groupings.
#' @param model.ranks Used to define which models to get dose-response curves for when \code{by.model = TRUE}.
#' @return An html file showing the outputs of \code{\link{summary.brsdata}}, \code{\link{summary.rjtrace}}, \code{\link{plot.rjtrace}}, and \code{\link{plot.dose_response}}.
#' @author Phil J. Bouchet
#' @seealso \code{\link{summary.brsdata}} \code{\link{summary.rjtrace}} \code{\link{plot.rjtrace}} \code{\link{plot.dose_response}}
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
#'                  
#' # Extract results
#' rj.post <- trace_rjMCMC(rj)
#' 
#' # Generate report
#' create_report()
#' }
#' @keywords brs dose-response rjmcmc

create_report <- function(outdir = getwd(),
                          filename = "espresso_report",
                          plot.height = 3,
                          model.ranks = 1){
  
  obj.n <- sapply(X = ls(envir = .GlobalEnv), FUN = function(x) class(get(x))[1]) %>% unlist()
  obj.ind <- purrr::map(.x = obj.n, .f = ~ any(.x %in% c("brsdata", "brsdata.grp", "rjtrace", "dose_response"))) %>% 
    unlist()
  obj.n <- obj.n[obj.ind]
  
  # Create temp directory
  out.dir <- file.path(getwd(), "tmp")
  dir.create(out.dir, showWarnings = FALSE)
  
  # Save objects
  save(list = names(obj.n), file = file.path(out.dir, "tmp_data.rda"))
  
  sink(file.path(out.dir, paste0(filename, ".Rmd")))
  
  writeLines(text = "---")
  writeLines(text = "title: \"espresso\"")
  writeLines(text = "subtitle: \"rjMCMC report\"")
  writeLines(text = "date: \"`r format(Sys.time(), '%d %B, %Y')`\"")
  writeLines(text = "output: html_document")
  writeLines(text = "---")
  
  # CSS
  writeLines(text = "```{css, echo=FALSE}")
  writeLines(text = ".scroll-100 {")
  writeLines(text = "max-height: 500px;")
  writeLines(text = "overflow-y: auto;")
  writeLines(text = "background-color: inherit;")
  writeLines(text = "}")
  writeLines(text = ".hljs-keyword {color: #63a35c;") 
  writeLines(text = " font-weight: normal;}")
  writeLines(text = ".hljs-string {color: #183691;")
  writeLines(text = ".hljs{color: #929292;}")
  writeLines(text = "```")
  
  writeLines(text = "```{r initial, include = FALSE}")
  writeLines(text = "library(espresso)")
  writeLines(text = "library(patchwork)")
  writeLines(text = paste0("load(\"", file.path(out.dir, "tmp_data.rda"), "\")"))
  writeLines(text = "```")
  
  writeLines(text = "## Results {.tabset}") # Interactive tabs
  
  if(sum(obj.n == "brsdata") > 0){
    
    writeLines(text = "### Data (original)")
    writeLines(text = "Original")
    writeLines(text = "```{r data, class.output=\"scroll-100\"}")
    writeLines(text = paste0("summary(", names(obj.n)[which(obj.n == "brsdata")], ")"))
    writeLines(text = "```")
  }
  
  if(sum(obj.n == "brsdata.grp") > 0){
  
  writeLines(text = "### Data (grouped)")
  writeLines(text = "```{r data_grouped, class.output=\"scroll-100\"}")
  writeLines(text = paste0("summary(", names(obj.n)[which(obj.n == "brsdata.grp")], ")"))
  writeLines(text = "```")
  }
  
  if(sum(obj.n == "rjtrace") > 0){
    
    writeLines(text = "### MCMC")
    writeLines(text = paste0("```{r mcmc, fig.height = ", plot.height, ", class.output=\"scroll-100\"}"))
    writeLines(text = paste0("summary(", names(obj.n)[which(obj.n == "rjtrace")], ", rmd = TRUE)"))
    writeLines(text = "```\n")
    
    writeLines(text = "### Trace")
    writeLines(text = paste0("```{r trace}"))
    writeLines(text = paste0("plot(", names(obj.n)[which(obj.n == "rjtrace")], ", autocorr = TRUE)"))
    writeLines(text = "```\n")
    
  }
  
  if(sum(obj.n == "dose_response") > 0){
    
    writeLines(text = "### Dose-response")
    writeLines(text = paste0("```{r doseResponse, class.output=\"scroll-100\"}"))
    
    for(oo in names(obj.n)[which(obj.n == "dose_response")]){
      tmp <- get(oo)
      if(tmp$by.model){
        if(length(model.ranks) == 1){
          writeLines(text = paste0("plot(", oo, ", model.rank = ", model.ranks, ")"))
        } else {
          for(j in model.ranks){writeLines(text = paste0("plot(", oo, ", model.rank = ", j, ")"))}
        }
      } else {
        writeLines(text = paste0("plot(", oo, ")")) 
      }
    }
    writeLines(text = "```\n")
  }
  
  sink()
  rmarkdown::render(input = file.path(out.dir, paste0(filename, ".Rmd")), output_dir = outdir)
  unlink(out.dir, recursive = TRUE)
  
}