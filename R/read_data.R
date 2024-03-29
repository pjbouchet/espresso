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
#'   \code{censored} \tab Integer variable indicating whether an observation is left-censored (-1), right-censored (1), or not censored (0) \cr
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
#' @importFrom Rdpack reprompt
#' @export
#' @param file Path to the input CSV file. If set to \code{NULL}, imports the \code{\link{example_brs}} data.
#' @param risk.functions Logical vector of length 4. If \code{TRUE}, indicates whether to include the risk functions derived in \insertCite{Moretti2014;textual}{espresso}, \insertCite{Jacobson2022;textual}{espresso}, \insertCite{Houser2013a;textual}{espresso}a or \insertCite{Houser2013b;textual}{espresso}b.
#' @param n.risk Vector of length 4. Number of samples to draw for simulating exposures from the dose-response curves selected using \code{risk.functions}. Defaults to \code{c(10, 10, 10, 10)}.
#' @param include.species Character vector specifying which species should be retained. These can be selected by any combination of scientific name, common name, or unique identifier, as listed in \code{\link{species_brs}}. All species are included when this argument is set to the default of \code{NULL}.
#' @param exclude.species Character vector specifying which species should be discarded. These can be selected by any combination of scientific name, common name, or unique identifier, as listed in \code{\link{species_brs}}. No species are excluded when this argument is set to the default of \code{NULL}.
#' @param min.N Minimum number of observations per species. Species with sample sizes smaller than \code{min.N} will be removed.
#' @param covariates Contextual covariates. Must be a character string containing one or more of the following: \code{exposed}, \code{sonar}, \code{behaviour} or \code{range}. No covariates are considered when this argument is set to \code{NULL}.
#' @param sonar.groups Named list detailing which sonar signals should be grouped a priori.
#' @param min.spl Minimum SPL value to be considered. Data below this value will be discarded. 
#' @param max.spl Maximum SPL value to be considered. Data above this value will be discarded.
#' @param dose.range Bounds for the dose-response function. Must be a vector of length 2. Defaults to: (1) a lower bound of 60 dB re 1μPa, taken as a conservative lower limit of detectability given hearing sensitivity and the lowest average sea noise conditions; and (2) an upper bound of 215 dB re 1μPa at/above which all animals are expected to respond. This upper bound is consistent with the maximum source levels employed in behavioural response studies (BRSs) to date. 
#' @param obs.sd Measurement uncertainty (expressed as a standard deviation in received levels), in dB re 1μPa.
#' @param verbose Logical. Whether to print or suppress warning messages.
#' 
#' @note Data simulated using \code{risk.functions} are assumed to represent measured responses and do not include left or right-censored observations. Note also that \code{include.species} and \code{exclude.species} override \code{risk.functions}, such that simulated data cannot be produced if the associated species is excluded. 
#' @return A list object of class \code{brsdata}, with the following elements:
#' 
#' \tabular{ll}{
#'   \code{dat} \tab Output dataset, after processing. \cr
#'   \code{species} \tab Species data, including species names, groups, sample sizes etc. \cr
#'   \code{whales} \tab Individual data, including animal IDs.  \cr
#'   \code{trials} \tab Exposure data, including exposure IDs. \cr
#'   \code{covariates} \tab Covariate data, including dummy coding, sonar groupings, factor levels etc. \cr
#'   \code{obs} \tab Observations, including left- and right-censoring cutoffs.  \cr  
#'   \code{param} \tab General parameters. \cr
#'  }
#' @import magrittr
#' @references
#' \insertAllCited{}
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
#'                   risk.functions = c(FALSE, FALSE, FALSE, FALSE),
#'                   exclude.species = "Sperm whale",
#'                   min.N = 2)
#' }

#' @keywords brs gvs rjmcmc dose-response             

read_data <- function(file = NULL,
                      risk.functions = c(FALSE, FALSE, FALSE, FALSE),
                      n.risk = c(10, 10, 10, 10),
                      include.species = NULL,
                      exclude.species = NULL,
                      min.N = NULL,
                      covariates = NULL,
                      sonar.groups = list(MFAS = c("MFAS", "MFA", "REAL MFA", "MFAS_DS", "MFA HELO",
                                                   "REAL MFA HELO", "SOCAL_d", 
                                                   "REAL MFAS", "REAL MFA (USS Dewey)"),
                                          LFAS = c("LFAS", "LFA", "LFAS_DS", "LFAS_LO",
                                                   "LFAS_122", "LFAS_185",
                                                   "HPAS-C", "HPAC-C", "HPAS-D", "HPASF-C",
                                                   "MPAS-C", "MPAS-D", "HPAC-C", "XHPAS-D",
                                                   "XHPAS-C")),
                      min.spl = NULL,
                      max.spl = NULL,
                      dose.range = c(60, 215),
                      obs.sd = 2.5,
                      verbose = TRUE){
  
  #' ---------------------------------------------
  # This function is for real data
  #' ---------------------------------------------
  
  simulation <- FALSE
  n.risk[!risk.functions] <- 0
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if(!is.null(covariates)){
    if(any(!covariates %in% c("exposed", "sonar", "behaviour", "range"))) stop("Covariate(s) not recognised. Options include: <exposed>, <sonar>, <behaviour>, and <range>.")
  }
  
  if(!is.null(sonar.groups)){
    if(!is.list(sonar.groups)) stop("sonar.groups must be a list")
    if(is.null(names(sonar.groups))) stop("Cannot find names for signal type groups")}
  
  if(is.null(covariates)){
    n.covariates <- 0
    covariate.names <- NULL
  } else {
    n.covariates <- length(covariates)
    covariate.names <- covariates
  }
  
  if(any(n.risk < 0)) stop("Number of samples must be positive!")
  
  covariate.types <- list(exposed = "f", sonar = "f", behaviour = "f", range = "d")
  covariate.types <- covariate.types[covariate.names] %>% unlist()
  
  #' ---------------------------------------------
  # Import the species list
  #' ---------------------------------------------
  
  species.list <- species_brs %>% 
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
      
      if(length(index) > 1) index <- max(index)
      
      species.list$new_code[sum(which(unlist(species.list[, index]) == n),
                                which(tolower(unlist(species.list[, index])) == n))]
    }))
    
    exclude.species <- unname(sapply(X = exclude.species, FUN = function(n) {
      index <- which(purrr::map_lgl(
        .x = seq_len(ncol(species.list)),
        .f = ~ any(n %in% c(
          tolower(unlist(species.list[, .x])),
          unlist(species.list[, .x])
        ))
      ))
      
      if(length(index) > 1) index <- max(index)
      
      species.list$new_code[sum(which(unlist(species.list[, index]) == n),
                                which(tolower(unlist(species.list[, index])) == n))]
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
    rawdat <- example_brs
    file <- "example_brs"
  } else {
    rawdat <- readr::read_csv(file = file, na = c(" ", "", "NA"), col_types = readr::cols())
  }
  
  rawdat <- rawdat %>% janitor::clean_names()
  
  #' ---------------------------------------------
  # Add samples from published risk functions, if necessary 
  #' ---------------------------------------------
  
  # Moretti et al. (2014)
  if(risk.functions[1]){
    
    spline.moretti <- stats::splinefun(x = moretti_dose$probability, y = moretti_dose$rl)
    pts.moretti <- spline.moretti(seq(0, 1, 0.001))
    rl.moretti <- sample(x = pts.moretti, size = n.risk[1], replace = TRUE)
    
    tbl.moretti <- tibble::tibble(resp_spl = rl.moretti, project = "Moretti_2014", species = "Md", 
                                  tag_id = paste0("Md_moretti_", sample(1:length(rl.moretti), replace = TRUE)),
                                  start_time = NA, resp_time = NA, resp_type = "foraging",
                                  resp_score = 7, resp_se_lcum = NA, max_spl = rl.moretti, max_se_lcum = NA,
                                  censored = 0, exp_order = 1, exp_duration = 0, exp_signal = "MFA",
                                  pre_feeding = TRUE, inferred_resp_range = NA, inferred_min_range = NA)
    
    tbl.moretti$resp_range <- tbl.moretti$min_range <- purrr::map_dbl(.x = rl.moretti,
                          .f = ~stats::optimize(f = range_finder, interval = c(0, 30000),
                                                SL = 215, target.L = .x)[['minimum']])
    
    tbl.moretti <- dplyr::relocate(tbl.moretti, names(rawdat)) %>% dplyr::arrange(tag_id)
    tbl.moretti <- split(tbl.moretti, tbl.moretti$tag_id) %>% 
      purrr::map(.x = ., .f = ~dplyr::mutate(.x, exp_order = dplyr::row_number())) %>% 
      do.call(rbind, .)

    rawdat <- dplyr::bind_rows(rawdat, tbl.moretti)
    
  }
  
  # Jacobson et al. (2022)
  if(risk.functions[2]){ 
    
    pmrf <- jacobson_dose %>% 
      dplyr::filter(max_rl >= dose.range[1] & max_rl <= dose.range[2]) %>% 
      dplyr::select(max_rl, q50) %>% dplyr::mutate(q50 = -q50)
    spline.pmrf <- stats::splinefun(x = pmrf$q50, y = pmrf$max_rl)
    pts.pmrf <- spline.pmrf(seq(min(pmrf$q50), max(pmrf$q50), 0.001))
    rl.pmrf <- sample(x = pts.pmrf, size = n.risk[2], replace = TRUE)
   
    tbl.pmrf <- tibble::tibble(resp_spl = rl.pmrf, project = "Jacobson_2020", species = "Md", 
                               tag_id = paste0("Md_jacobson_", sample(1:length(rl.pmrf), replace = TRUE)),
                               start_time = NA, resp_time = NA, resp_type = "avoidance",
                               resp_score = 7, resp_se_lcum = NA, max_spl = rl.pmrf, max_se_lcum = NA,
                               censored = 0, exp_order = 1, exp_duration = 0, exp_signal = "MFAS",
                               pre_feeding = FALSE, inferred_resp_range = NA, inferred_min_range = NA)
    
    tbl.pmrf$resp_range <- tbl.pmrf$min_range <- 
      purrr::map_dbl(.x = rl.pmrf, .f = ~stats::optimize(f = range_finder, interval = c(0, 30000), SL = 215, target.L = .x)[['minimum']])
    
    tbl.pmrf <- dplyr::relocate(tbl.pmrf, names(rawdat)) %>% dplyr::arrange(tag_id)
    tbl.pmrf <- split(tbl.pmrf, tbl.pmrf$tag_id) %>% 
      purrr::map(.x = ., .f = ~dplyr::mutate(.x, exp_order = dplyr::row_number())) %>% 
      do.call(rbind, .)
    
    rawdat <- dplyr::bind_rows(rawdat, tbl.pmrf)
    
  }
    
  # Houser et al. (2013a)
  if(risk.functions[3]){
    
    spline.houser.Zca <- stats::splinefun(x = houser_dose$Zca$probability, y = houser_dose$Zca$rl)
    pts.houser.Zca <- spline.houser.Zca(seq(0, 1, 0.001))
    rl.houser.Zca <- sample(x = pts.houser.Zca, size = n.risk[3], replace = TRUE)

    tbl.houser.Zca <- tibble::tibble(
      resp_spl = rl.houser.Zca, 
      project = "Houser_2013a",
      species = "Zca",
      tag_id = paste0("Zca_houser2013a_", sample(1:length(rl.houser.Zca), replace = TRUE)),
      start_time = NA, resp_time = NA, 
      resp_type = sample(x = c("haul out", "refusal to participate"), size = n.risk[3], replace = TRUE),
      resp_score = 7,
      resp_se_lcum = NA, max_spl = rl.houser.Zca, max_se_lcum = NA, censored = 0, 
      exp_order = 1, exp_duration = 0, exp_signal = "MFA",
      pre_feeding = FALSE, resp_range = 0.001, min_range = 0.001, 
      inferred_resp_range = 0.001, inferred_min_range = 0.001)
    
    tbl.houser.Zca <- dplyr::relocate(tbl.houser.Zca, names(rawdat)) %>% dplyr::arrange(tag_id)
    tbl.houser.Zca <- split(tbl.houser.Zca, tbl.houser.Zca$tag_id) %>% 
      purrr::map(.x = ., .f = ~dplyr::mutate(.x, exp_order = dplyr::row_number())) %>% 
      do.call(rbind, .)
    
    rawdat <- dplyr::bind_rows(rawdat, tbl.houser.Zca)
    
  }
  
  # Houser et al. (2013b)
  if(risk.functions[4]){
    
    spline.houser.Tt <- stats::splinefun(x = houser_dose$Tt$probability, y = houser_dose$Tt$rl)
    pts.houser.Tt <- spline.houser.Tt(seq(0, 1, 0.001))
    rl.houser.Tt <- sample(x = pts.houser.Tt, size = n.risk[4], replace = TRUE)
    
    tbl.houser.Tt <- tibble::tibble(
      resp_spl = rl.houser.Tt, 
      project = "Houser_2013b",
      species = "Tt",
      tag_id = paste0("Tt_houser2013b_", sample(1:length(rl.houser.Tt), replace = TRUE)),
      start_time = NA, 
      resp_time = NA, 
      resp_type = sample(x = c("tail slap", "refusal to participate"), size = n.risk[4], replace = TRUE),
      resp_score = 7,
      resp_se_lcum = NA, max_spl = rl.houser.Tt, max_se_lcum = NA, censored = 0, 
      exp_order = 1, exp_duration = 0, exp_signal = "MFA",
      pre_feeding = FALSE, resp_range = 0.001, min_range = 0.001, 
      inferred_resp_range = 0.001, inferred_min_range = 0.001)
    
    tbl.houser.Tt <- dplyr::relocate(tbl.houser.Tt, names(rawdat)) %>% dplyr::arrange(tag_id)
    tbl.houser.Tt <- split(tbl.houser.Tt, tbl.houser.Tt$tag_id) %>% 
      purrr::map(.x = ., .f = ~dplyr::mutate(.x, exp_order = dplyr::row_number())) %>% 
      do.call(rbind, .)
    
    rawdat <- dplyr::bind_rows(rawdat, tbl.houser.Tt)
    
  }
  
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
  brsdat <- rawdat %>% 
    dplyr::mutate(species = tools::toTitleCase(species)) %>% 
    dplyr::left_join(x = ., y = species.list, by = c("species" = "code")) %>% 
    dplyr::select(-species) %>% 
    dplyr::rename(species = new_code) 
  
  if(!is.null(covariate.names)){
    
    if("exposed" %in% covariate.names) {
      brsdat <- brsdat %>% 
        dplyr::mutate(exp_order = ifelse(exp_order == 1, 0, 1)) %>% 
        dplyr::rename(exposed = exp_order)}
    
    if("sonar" %in% covariate.names){
      
      # Determine which sonar signal types to exclude
      exclude.sonar <- unique(brsdat$exp_signal)[!unique(brsdat$exp_signal) %in% unname(unlist(sonar.groups))]
      null.sonar <- unname(unlist(sonar.groups))[!unname(unlist(sonar.groups)) %in% unique(brsdat$exp_signal)]
      
      
      brsdat <- brsdat %>% 
        dplyr::filter(!exp_signal %in% c(exclude.sonar, null.sonar)) 
      
      signal.df <- sonar.groups %>% 
        tibble::enframe() %>% 
        tidyr::unnest(cols = c(value)) %>% 
        dplyr::rename(sonar = name, signal = value)
      
      sonar.groups <- purrr::map(.x = sonar.groups, .f = ~.x[!.x %in% null.sonar])
      
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
            # Not censored (whale did respond), response range unknown, estimate available
            censored == 0 & is.na(resp_range) & !is.na(inferred_resp_range) ~ inferred_resp_range,
            # Not censored (whale did respond), response range unknown, estimate not available
            censored == 0 & is.na(resp_range) & is.na(inferred_resp_range) ~ min_range,
            censored == 0 & is.na(resp_range) & is.na(inferred_resp_range) & is.na(min_range) ~ inferred_min_range,
            
            # Left-censored
            censored == -1 & !is.na(resp_range) ~ resp_range,
            censored == -1 & is.na(resp_range) & !is.na(inferred_resp_range) ~ inferred_resp_range,
            censored == -1 & is.na(resp_range) & is.na(inferred_resp_range) ~ min_range,
            censored == -1 & is.na(resp_range) & is.na(inferred_resp_range) & is.na(min_range) ~ inferred_min_range,
            
            # Right -censored
            censored == 1 & !is.na(min_range) ~ min_range,
            censored == 1 & is.na(min_range) & !is.na(inferred_min_range) ~ inferred_min_range,
            TRUE ~ NA_real_
          ))
    
  }
  
  cols.to.extract <- c("project", "species", "scientific_name", 
                       "common_name", "tag_id", covariate.names, "spl", "Lc", "Rc", "censored")
  
  brsdat <- brsdat %>% 
    dplyr::mutate(Lc = ifelse(censored == -1, resp_spl, NA)) %>% 
    dplyr::rename(spl = resp_spl, Rc = max_spl) %>% 
    dplyr::select_at(., tidyselect::all_of(cols.to.extract)) %>% 
    dplyr::arrange(species, tag_id)
  
  #' ---------------------------------------------
  # Covariate levels
  #' ---------------------------------------------
  
  if(n.covariates > 0){
    
    covariates.df <- brsdat %>% 
      dplyr::select_at(tidyselect::all_of(covariate.names)) %>% 
      data.frame() %>%
      dplyr::mutate(dplyr::across(covariate.names[covariate.types == "f"], ~ as.factor(.x)))
    
    fL <- sapply(X = covariate.names, 
                 FUN = function(x) factor_levels(covname = x, dat = covariates.df), 
                 simplify = FALSE, USE.NAMES = TRUE)
    
    # Dummy coding
    dummy.cov <- purrr::map(
      .x = seq_len(n.covariates),
      .f = ~ {
        if (covariate.types[.x] == "f") {
          fastDummies::dummy_cols(
            .data = covariates.df,
            select_columns = covariate.names[.x],
            remove_first_dummy = FALSE,
            remove_selected_columns = TRUE
          ) %>%
            dplyr::select(-tidyselect::any_of(covariate.names))
        } else {
          covariates.df[, .x, drop = FALSE]
        }
      }
    ) %>%
      purrr::set_names(., covariate.names)
    
    dummy.df <- t(do.call(cbind, dummy.cov))
    rownames(dummy.df) <- unlist(purrr::map(.x = covariate.names, 
                                            .f = ~paste0(.x, if(fL[[.x]]$nL == 0) "" else paste0("_", fL[[.x]][["Lnames"]]))))
    
  }
  
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
    dplyr::mutate(r = dplyr::row_number()) %>%
    dplyr::filter(species %in% include.species) %>% 
    dplyr::filter(!species %in% exclude.species) %>% 
    dplyr::mutate(sp_orig = species) 
  
  #' ---------------------------------------------
  # SPL filter
  #' ---------------------------------------------
  
  if(!is.null(min.spl)){ brsdat <- brsdat %>% dplyr::filter(spl > min.spl | is.na(spl)) }
  if(!is.null(max.spl)){ brsdat <- brsdat %>% dplyr::filter(spl < max.spl | is.na(spl)) }
  
  #' ---------------------------------------------
  # Covariates
  #' ---------------------------------------------
  
  if(n.covariates > 0){
    
    covariates.df <- covariates.df[brsdat$r, , drop = FALSE]
    dummy.df <- dummy.df[, brsdat$r, drop = FALSE]
    dummy.cov <- purrr::map(.x = dummy.cov, .f = ~dplyr::slice(.x, brsdat$r))
    
  }

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
                                 tmp <- brsdat %>% dplyr::filter(tag_id == .x)
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
                           FUN = function(x) rep(species.id[x], n.trials.per.whale[x]))
  if("list" %in% class(species.trials)) species.trials <- do.call(c, species.trials)
  
  #' ---------------------------------------------
  # Sampling uncertainty
  #' ---------------------------------------------
  
  measurement.precision <- 1 / (obs.sd^2)
  
  #' ---------------------------------------------
  # Observations
  #' ---------------------------------------------
  
  y_ij <- brsdat$spl
  is.censored <- brsdat$censored
  
  # Right-censoring
  max.dose <- rep(dose.range[1], n.trials)
  y_ij[is.censored == 1] <- NA
  max.dose[is.censored == 1] <- brsdat[is.censored == 1, ]$Rc
  
  # Left-censoring
  min.dose <- rep(dose.range[2], n.trials)
  y_ij[is.censored == -1] <- NA
  min.dose[is.censored == -1] <- brsdat[is.censored == -1, ]$Lc
  
  #' ---------------------------------------------
  # Create a species summary
  #' ---------------------------------------------
  suppressWarnings(species.summary <- brsdat %>% 
                     dplyr::group_by(common_name) %>% 
                     dplyr::summarise(N_ind = length(unique(tag_id)), 
                                      N_trials = dplyr::n(),
                                      censored.L = sum(censored == -1), 
                                      censored.R = sum(censored == 1), 
                                      mean = mean(spl, na.rm = TRUE),
                                      min = min(spl, na.rm = TRUE),
                                      max = max(spl, na.rm = TRUE), .groups = "keep") %>% 
                     dplyr::ungroup())
  
  suppressMessages(species.summary <- 
       dplyr::left_join(x = species.summary, 
                        y = dplyr::distinct(brsdat[, c("common_name", "species")])) %>% 
                     dplyr::select(species, common_name, N_ind, N_trials, censored.L,
                                   censored.R, mean, min, max))
  
  
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
                      signal.types = if("sonar" %in% covariate.names) append(sonar.groups, list(exclude = exclude.sonar, not_found = null.sonar)) else NULL),
    # Observations
    obs = list(y_ij = y_ij,
               censored = is.censored,
               Lc = min.dose,
               Rc = max.dose,
               prec = measurement.precision,
               sd = obs.sd,
               risk.functions = risk.functions,
               n.risk = n.risk),
    # Parameters
    param = list(sim = simulation,
                 data.file = basename(file),
                 # bounds = param.bounds,
                 dose.range = dose.range))
  
  if(n.covariates > 0) {
    res$covariates <- append(res$covariates, list(df = covariates.df, dummy = dummy.cov, dummy.df = dummy.df))
    res$covariates$fL <- fL
    
  }
  
  if(verbose){
  if(any(res$species$summary$N_ind < 5)) 
    warning("N < 5 for some species. Species with low sample sizes may need to be grouped a priori.")}
  
  class(res) <- c("brsdata", class(res))
  return(res) 
  
}