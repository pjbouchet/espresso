#' Import dose-response data
#'
#' Read in and format data from behavioural response studies (BRSs). The function provides options for filtering the data based on species and/or sample sizes. See Details for formatting requirements.
#'
#' The input data file must contain the following fields:
#'
#' \tabular{ll}{
#'   \code{species} \tab Species code, as listed in \code{\link{species_brs}}\cr
#'   \code{tag_id} \tab Unique tag identifier\cr
#'   \code{resp_SPL} \tab Sound pressure level at time of response (in dB re 1μPa) \cr
#'   \code{max_SPL} \tab Maximum sound pressure level reached during the exposure\cr
#'   \code{censored} \tab Binary variable indicating whether an observation is right-censored (1) or not (0) \cr
#'  }
#'  
#' When covariates are specified, the below fields must also be included, as relevant:

#' \tabular{ll}{
#'   \code{exp_order} \tab History of exposure (1 = 1st exposure, 2 = 2nd exposure, etc.)\cr
#'   \code{exp_signal} \tab Type of sonar signal (e.g., MFAS, REAL MFA, PRN, CAS)\cr
#'   \code{pre_feeding} \tab Behavioural mode (TRUE = feeding, FALSE = non-feeding)\cr
#'   \code{min_range} \tab Minimum whale-source range during the exposure\cr
#'   \code{resp_range} \tab Whale-source range at the time of response\cr
#'   \code{inferred_resp_range} \tab Best estimate of whale-source range at the time of response\cr
#'   \code{inferred_min_range} \tab Best estimate of minimum whale-source range during the exposure\cr
#'  } 
#'
#' An example data file is included in the package. See \code{\link{example_brs}} for details.
#' 
#' @export
#' @param file Path to the input CSV file. If set to \code{NULL}, imports the \code{\link{example_brs}} data.
#' @param include.species Character vector specifying which species should be retained. These can be selected by any combination of scientific name, common name, or unique identifier, as listed in \code{\link{species_brs}}. All species are included when this argument is set to the default of \code{NULL}.
#' @param exclude.species Character vector specifying which species should be discarded. These can be selected by any combination of scientific name, common name, or unique identifier, as listed in \code{\link{species_brs}}. No species are excluded when this argument is set to the default of \code{NULL}.
#' @param min.N Minimum number of observations per species. Species with sample sizes smaller than \code{min.N} will be removed.
#' @param covariates Contextual covariates. Must be a character string containing one or more of the following: \code{exposed}, \code{sonar}, \code{behaviour} or \code{range}. No covariates are considered when this argument is set to \code{NULL}.
#' @param sonar.groups Named list detailing which sonar signals should be grouped a priori.
#' @param exclude.sonar Character vector specifying which sonar signals should be excluded.
#' @param range.dB Bounds for the dose-response function. Must be a vector of length 2. Defaults to: (1) a lower bound of 60 dB re 1μPa, taken as a conservative lower limit of detectability given hearing sensitivity and the lowest average sea noise conditions; and (2) an upper bound of 215 dB re 1μPa at/above which all animals are expected to respond. This upper bound is consistent with the maximum source levels employed in behavioural response studies (BRSs) to date. 
#' @param range.var Permissible range for the variance term(s) phi (between-whale variance) and sigma (between-exposure variance). Must be a vector of length 2. 
#' @param obs.sd Measurement uncertainty (expressed as a standard deviation in received levels), in dB re 1μPa.
#' @param verbose Logical. Whether to print or suppress warning messages.
#' 
#' @return A list object of class \code{brsdata}, with the following elements:
#' 
#' \tabular{ll}{
#'   \code{dat} \tab Output dataset, after processing. \cr
#'   \code{species} \tab Species data, including species names, groups, sample sizes etc. \cr
#'   \code{whales} \tab Individual data, including animal IDs.  \cr
#'   \code{trials} \tab Exposure data, including exposure IDs. \cr
#'   \code{covariates} \tab Covariate data, including dummy coding, sonar groupings, factor levels etc. \cr
#'   \code{obs} \tab Observations, including right-censoring cutoffs.  \cr  
#'   \code{param} \tab General parameters. \cr
#'  }
#' @import magrittr
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.brsdata}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Import the example data
#' mydat <- read_data(file = NULL) 
#' 
#' # Import a real dataset with the sonar and range covariates, 
#' # whilst excluding sperm whales and any other species with a sample size
#' # smaller than two observations.
#' mydat <- read_data(file = "path/to/my/data.csv", 
#'                   exclude.species = "Sperm whale",
#'                   min.N = 2)
#' }
#' @keywords brs gvs rjmcmc dose-response             

read_data <- function(file = NULL,
                      include.species = NULL, 
                      exclude.species = NULL,
                      min.N = NULL,
                      covariates = NULL,
                      sonar.groups = list(MFAS = c("MFAS", "MFA", "REAL MFA", "MFAS_DS", "MFA HELO",
                                                   "REAL MFA HELO", "SOCAL_d", "REAL MFAS", "MF ALARM"),
                                          LFAS = NULL),
                      exclude.sonar = c("PRN", "CAS", "ALARM", "HF ALARM", "HFAS", 
                                        "HF ALARM CAS", "MFAS CAS"),
                      range.dB = c(60, 215), 
                      range.var = c(0, 45),
                      obs.sd = 2.5,
                      verbose = TRUE){
  
  #' ---------------------------------------------
  # This function is for real data
  #' ---------------------------------------------
  
  simulation <- FALSE
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if(!is.null(covariates)){
    if(any(!covariates %in% c("exposed", "sonar", "behaviour", "range"))) stop("Covariate(s) not recognised. Options include: <exposed>, <sonar>, <behaviour>, and <range>.")
  }
  
  if(!is.null(sonar.groups)){
    if(!is.list(sonar.groups)) stop("signal.types must be a list")
    if(is.null(names(sonar.groups))) stop("Cannot find names for signal type groups")}
  
  signal.null <- purrr::map_lgl(.x = sonar.groups, .f = ~is.null(.x))
  if(sum(signal.null) > 1) stop("Only one NULL element allowed in <signal.types>")
  
  if(is.null(covariates)){
    n.covariates <- 0
    covariate.names <- NULL
  } else {
    n.covariates <- length(covariates)
    covariate.names <- covariates
  }
  
  covariate.types <- list(exposed = "f", sonar = "f", behaviour = "f", range = "d")
  covariate.types <- covariate.types[covariate.names] %>% unlist()
  
  #' ---------------------------------------------
  # Import the species list
  #' ---------------------------------------------
  
  species.list <- espresso::species_brs %>% 
    janitor::clean_names(.)
  
  #' ---------------------------------------------
  # Compile the lists of species to be included / excluded
  #' ---------------------------------------------
  
  if (any(!include.species %in% c(
    tolower(species.list$common_name),
    species.list$common_name,
    species.list$new_code,
    species.list$scientific_name
  ))) {
    
    stop("Species not recognised")
    
  } else {
    
    include.species <- unname(sapply(X = include.species, FUN = function(n) {
      index <- which(purrr::map_lgl(
        .x = seq_len(ncol(species.list)),
        .f = ~ any(n %in% c(
          tolower(unlist(species.list[, .x])),
          unlist(species.list[, .x])
        ))
      ))
      
      species.list$new_code[which(unlist(species.list[, index]) == n |
                                    tolower(unlist(species.list[, index])) == n)]
    }))
    
    exclude.species <- unname(sapply(X = exclude.species, FUN = function(n) {
      index <- which(purrr::map_lgl(
        .x = seq_len(ncol(species.list)),
        .f = ~ any(n %in% c(
          tolower(unlist(species.list[, .x])),
          unlist(species.list[, .x])
        ))
      ))
      
      species.list$new_code[which(unlist(species.list[, index]) == n |
                                    tolower(unlist(species.list[, index])) == n)]
    }))
  }
  
  if("matrix" %in% class(exclude.species)){
    exclude.species <- exclude.species[complete.cases(exclude.species), ]
  }else{
    exclude.species <- exclude.species[!is.na(exclude.species)]
  }
  
  #' ---------------------------------------------
  # Import and process the data
  #' ---------------------------------------------
  
  if(is.null(file)){
    rawdat <- espresso::example_brs
    file <- "example_brs"
  } else {
    rawdat <- readr::read_csv(file = file, na = c(" ", "NA"), col_types = cols())
  }
  
  rawdat <- rawdat %>% janitor::clean_names()
  
  #' ---------------------------------------------
  # Determine whether observations are right-censored
  #' ---------------------------------------------
  # severity.score <- 4
  # rawdat <- rawdat %>% dplyr::mutate(rc = ifelse(is.na(resp_score),
  #                                                ifelse(!is.na(resp_spl), 0, 1), 
  #                                                ifelse(resp_score < severity.score, 1, 0)))
  
  #' ---------------------------------------------
  # Extract relevant columns
  #' ---------------------------------------------
  
  cols.to.extract <- c("project", "species", "tag_id", "resp_spl", "max_spl", "censored")
  
  if(!is.null(covariate.names)){
    if("exposed" %in% covariate.names) cols.to.extract <- c(cols.to.extract, "exp_order")
    if("sonar" %in% covariate.names) cols.to.extract <- c(cols.to.extract, "exp_signal")
    if("behaviour" %in% covariate.names) cols.to.extract <- c(cols.to.extract, "pre_feeding")
    if("range" %in% covariate.names) cols.to.extract <- c(cols.to.extract, "min_range", "resp_range",
                                                          "inferred_resp_range", "inferred_min_range")
  }
  
  rawdat <- rawdat %>% dplyr::select_at(tidyselect::all_of(cols.to.extract))
  
  #' ---------------------------------------------
  # Essential formatting
  #' ---------------------------------------------
  
  # Ensures data are in correct format
  rawdat$max_spl <- as.numeric(rawdat$max_spl)
  rawdat$censored <- as.integer(rawdat$censored)
  
  # Recode species and covariates
  brsdat <- rawdat %>% dplyr::mutate(species = tools::toTitleCase(species)) %>% 
    dplyr::left_join(x = ., y = species.list, by = c("species" = "code")) %>% 
    dplyr::select(-species) %>% 
    dplyr::rename(species = new_code) 
  
  if(!is.null(covariate.names)){
    
    if("exposed" %in% covariate.names) {
      brsdat <- brsdat %>% 
        dplyr::mutate(exp_order = ifelse(exp_order == 1, 0, 1)) %>% 
        dplyr::rename(exposed = exp_order)}
    
    if("sonar" %in% covariate.names){
      
      if(!any(exclude.sonar %in% unique(brsdat$exp_signal))) 
        stop("Unrecognised sonar type(s) in <exclude.signals>")
      
      if(!any(unlist(sonar.groups) %in% unique(brsdat$exp_signal))) 
        stop("Unrecognised sonar type(s) in <exclude.signals>")
      
      brsdat <- brsdat %>% 
        dplyr::filter(!exp_signal %in% exclude.sonar) 
      
      if(sum(signal.null) > 0){
        sonar.groups[[which(signal.null)]] <- 
          unique(unname(brsdat$exp_signal[!brsdat$exp_signal %in% unlist(sonar.groups[-which(signal.null)])]))
      } else {
        if(length(outersect(c(exclude.sonar, unlist(sonar.groups)), unique(brsdat$exp_signal))) > 0){
          for(j in outersect(c(exclude.sonar, unlist(sonar.groups)), unique(brsdat$exp_signal))){
            sonar.groups[j] <- j
          }
        }
      }
      
      signal.df <- sonar.groups %>% 
        tibble::enframe() %>% 
        tidyr::unnest(cols = c(value)) %>% 
        dplyr::rename(sonar = name, signal = value)
      
      brsdat <- brsdat %>% 
        dplyr::left_join(x = ., y = signal.df, by = c("exp_signal" = "signal"))
    }
    
    if("behaviour" %in% covariate.names) brsdat <- brsdat %>% 
        dplyr::mutate(pre_feeding = ifelse(pre_feeding, "feeding", "non-feeding")) %>% 
        dplyr::rename(behaviour = pre_feeding)
    
    #' ---------------------------------------------
    # Extract the correct range information
    #' ---------------------------------------------
    
    if("range" %in% covariate.names) brsdat <- brsdat %>% 
        dplyr::mutate(
          range = dplyr::case_when(
            # Not censored (whale did respond) and response range known
            censored == 0 & !is.na(resp_range) ~ resp_range,
            # Not censored (whale did respond) and response range unknown but estimate available
            censored == 0 & is.na(resp_range) & !is.na(inferred_resp_range) ~ inferred_resp_range,
            # Not censored (whale did respond) and response range unknown but estimate available
            censored == 0 & is.na(resp_range) & !is.na(inferred_resp_range) ~ inferred_resp_range,
            censored == 1 & !is.na(min_range) ~ min_range,

            censored == 1 & is.na(min_range) & !is.na(inferred_min_range) ~ inferred_min_range,
            TRUE ~ NA_real_
          ))
    
    
  }
  
  cols.to.extract <- c("project", "species", "scientific_name", 
                       "common_name", "tag_id", covariate.names, "spl", "Rc", "censored")
  
  brsdat <- brsdat %>% 
    dplyr::rename(spl = resp_spl, Rc = max_spl)%>% 
    dplyr::select_at(., tidyselect::all_of(cols.to.extract)) %>% 
    dplyr::arrange(species, tag_id)
  
  #' ---------------------------------------------
  # Remove empty records
  #' ---------------------------------------------
  
  brsdat <- brsdat %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(na.val = sum(is.na(spl), is.na(Rc))) %>% 
    dplyr::filter(na.val < 2) %>% 
    dplyr::select(-na.val) %>% 
    dplyr::ungroup()
  
  #' ---------------------------------------------
  # Filter out small sample sizes
  #' ---------------------------------------------
  
  if(!is.null(min.N)){
    
    brsdat <- brsdat %>% 
      dplyr::group_by(species) %>% 
      dplyr::mutate(nper = dplyr::n_distinct(tag_id)) %>% 
      dplyr::ungroup()
    
    exclude.species <- brsdat %>% 
      dplyr::filter(nper < min.N) %>% 
      dplyr::pull(species) %>% 
      unique(.) %>%
      c(exclude.species, .)
    
  }
  
  if(length(include.species) == 0) include.species <- unique(brsdat$species)
  if(length(exclude.species) == 0) exclude.species <- outersect(x = include.species, y = unique(brsdat$species))
  
  #' ---------------------------------------------
  # Subset data by species
  #' ---------------------------------------------
  
  brsdat <- brsdat %>% 
    dplyr::filter(species %in% include.species) %>% 
    dplyr::filter(!species %in% exclude.species) %>% 
    dplyr::mutate(r = dplyr::row_number()) %>% 
    dplyr::mutate(sp_orig = species) 
  
  #' ---------------------------------------------
  # Summarise dataset
  #' ---------------------------------------------
  
  # Number of species
  n.species <- length(unique(brsdat$species)) 
  
  # Species names
  sp.names <- sort(unique(brsdat$species))
  
  # Number of whales
  n.whales <- length(unique(brsdat$tag_id))
  
  # Species IDs
  species.id <- purrr::map_dbl(.x = unique(brsdat$tag_id),
                               .f = ~{
                                 tmp <- brsdat %>% 
                                   dplyr::filter(tag_id == .x)
                                 which(sp.names == unique(tmp$species))})
  
  # Number of individuals per species
  n.per.species <- as.numeric(table(species.id))
  
  # Total number of exposures
  n.trials <- nrow(brsdat)
  
  # Number of exposures per animal
  # Use factor trick here to conserve order of tag_id
  n.trials.per.whale <- brsdat %>% 
    dplyr::mutate(tag_f = factor(tag_id, levels = unique(tag_id))) %>% 
    dplyr::group_by(tag_f) %>% 
    dplyr::count(.) %>% 
    dplyr::pull(n)
  
  # Whale IDs
  whale.id <- rep(seq_len(n.whales), n.trials.per.whale)
  
  # Species trials
  species.trials <- sapply(X = seq_along(species.id), 
                           FUN = function(x) rep(species.id[x], n.trials.per.whale[x])) %>% do.call(c, .)
  
  # Bounds for model parameters
  lower.bound <- list(
    y.ij = range.dB[1],
    t.ij = range.dB[1],
    mu = range.dB[1], 
    mu.i = range.dB[1], 
    nu = range.dB[1],
    phi = range.var[1], 
    sigma = range.var[1],
    tau = range.var[1])
  
  upper.bound <- list(
    y.ij = range.dB[2],
    t.ij = range.dB[2],
    mu = range.dB[2], 
    mu.i = range.dB[2],
    nu = range.dB[2],
    phi = range.var[2],
    sigma = range.var[2],
    tau = range.var[2])
  
  for(j in seq_len(n.covariates)) lower.bound[covariate.names[j]] <- -Inf
  for(j in seq_len(n.covariates)) upper.bound[covariate.names[j]] <- Inf
  
  param.bounds <- cbind(purrr::map_df(.x = lower.bound, .f = ~.x) %>% t(.),
                        purrr::map_df(.x = upper.bound, .f = ~.x) %>% t(.))
  
  #' ---------------------------------------------
  # Sampling uncertainty
  #' ---------------------------------------------
  
  measurement.precision <- 1 / (obs.sd^2)
  
  #' ---------------------------------------------
  # Covariates
  #' ---------------------------------------------
  
  if(n.covariates > 0){
    
    covariates.df <- brsdat %>% 
      dplyr::select_at(tidyselect::all_of(covariate.names)) %>% 
      data.frame(.)
    
    for(ii in 1:length(covariate.names)){
      if(covariate.types[ii] == "f") covariates.df[, covariate.names[ii]] <- as.factor(covariates.df[, covariate.names[ii]])
      if(covariate.types[ii] %in% c("d", "i")) covariates.df[, covariate.names[ii]] <- as.numeric(covariates.df[, covariate.names[ii]])}
    
    # Dummy coding
    dummy.cov <- purrr::map(.x = seq_len(n.covariates),
                            .f = ~{
                              if(covariate.types[.x] == "f"){
                                
                                fastDummies::dummy_cols(.data = covariates.df, 
                                                        select_columns = covariate.names[.x],
                                                        remove_first_dummy = FALSE, 
                                                        remove_selected_columns = TRUE) %>% 
                                  dplyr::select(-tidyselect::any_of(covariate.names))
                                
                              } else { covariates.df[, .x, drop = FALSE] }}) %>% 
      purrr::set_names(., covariate.names)
    
  } # End if n.covariates
  
  #' ---------------------------------------------
  # Observations
  #' ---------------------------------------------
  
  y_ij <- brsdat$spl
  is.censored <- brsdat$censored
  
  # Right-censoring
  max.dose <- rep(lower.bound$mu, n.trials)
  y_ij[is.censored == 1] <- NA
  max.dose[is.censored == 1] <- brsdat[is.censored == 1, ]$Rc
  
  #' ---------------------------------------------
  # Create a species summary
  #' ---------------------------------------------
  suppressWarnings(species.summary <- brsdat %>% 
                     dplyr::group_by(common_name) %>% 
                     dplyr::summarise(N_ind = length(unique(tag_id)), 
                                      N_trials = dplyr::n(),
                                      censored = sum(is.na(spl)), 
                                      mean = mean(spl, na.rm = TRUE),
                                      min = min(spl, na.rm = TRUE),
                                      max = max(spl, na.rm = TRUE), .groups = "keep") %>% 
                     dplyr::ungroup())
  
  suppressMessages(species.summary <- 
                     dplyr::left_join(x = species.summary, 
                                      y = dplyr::distinct(brsdat[, c("common_name", "species")])) %>% 
                     dplyr::select(species, common_name, N_ind, N_trials, censored, mean, min, max))
  
  
  #' ---------------------------------------------
  # Store results in list
  #' ---------------------------------------------
  
  
  # Factor levels and other associated info for each covariate
  # nL: Number of levels (will be 0 for integer and continuous covariates)
  # nparam: Number of levels other than the baseline
  # index: column indices for those levels.
  
  res <- list(
    
    # Data
    ddf = brsdat,
    
    # Species
    species = list(names = sp.names,
                   groups = NULL,
                   abbrev = NULL,
                   n = n.species,
                   nper = n.per.species,
                   id = species.id,
                   summary = species.summary,
                   trials = species.trials,
                   exclude = exclude.species),
    
    # Individuals
    whales = list(n = n.whales,
                  id = whale.id),
    # Exposures
    trials = list(n = n.trials,
                  nper = n.trials.per.whale),
    
    # Contextual factors
    covariates = list(n = n.covariates,
                      names = covariate.names,
                      signal.types = append(sonar.groups, list(exclude = exclude.sonar))),
    # Observations
    obs = list(y_ij = y_ij,
               censored = is.censored,
               Rc = max.dose,
               prec = measurement.precision,
               sd = obs.sd),
    # Parameters
    param = list(sim = simulation,
                 data.file = file,
                 bounds = param.bounds,
                 dose.range = range.dB))
  
  if(n.covariates > 0) {
    res$covariates <- append(res$covariates, 
                             list(df = covariates.df,
                                  dummy = dummy.cov))
    
    res$covariates$fL <- sapply(X = covariate.names, 
                                FUN = function(x) factor_levels(covname = x, dat = res$covariates$df), 
                                simplify = FALSE, USE.NAMES = TRUE)
    
  }
  
  if(verbose){
  if(any(res$species$summary$N_ind < 5)) 
    warning("N < 5 for some species. Species with low sample sizes may need to be grouped a priori.")}
  
  class(res) <- c("brsdata", class(res))
  return(res) 
  
}