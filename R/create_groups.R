#' Create groupings
#'
#' Assign species to *a priori* groupings. 
#' 
#' @export
#' @param dat.obj Input dataset. Must be an object of class \code{brsdata}, as returned by \code{\link{read_data}} or \code{\link{simulate_data}}. 
#' @param species.groups List indicating which species should be grouped. These can be given using any combination of scientific name, common name, or unique identifier, as defined in \code{\link{species_brs}}. If \code{species.groups} is given as a named list, the names of each list element will be used as group names. If \code{species.groups} is unnamed, the function will combine species codes from \code{\link{species_brs}} to generate appropriate group names.
#' @param abbrev Logical. If \code{TRUE}, group names will be abbreviated. This is useful for plotting and data summaries. 
#' 
#' @return A list object of class \code{brsdata.grp}.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{undo_groups}} \code{\link{summary.brsdata}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Import the example data, excluding Risso's dolphins and minke whales
#' mydat <- read_data(file = NULL, exclude.species = c("gg", "ba")) 
#' summary(mydat)
#' 
#' # Group all beaked whales together
#' mydat.grouped <- create_groups(mydat, 
#'                   species.groups = list(c("Md", "Cuvier's beaked whale", 
#'                                         "Hyperoodon ampullatus")))
#' 
#' Examine the resulting grouping
#' summary(mydat.grouped)
#' }
#' @keywords brs rjmcmc dose-response

create_groups <- function(dat.obj, 
                          species.groups = NULL, 
                          abbrev = FALSE){
  
  #' ---------------------------------------------
  # Perform function checks
  #' ---------------------------------------------
  
  if(!"brsdata" %in% class(dat.obj)) stop("Input must be of class <brsdata>")
  if(!is.null(species.groups) & !is.list(species.groups)) stop("species.groups must be a list")
  if(!all(unname(unlist(species.groups)) %in% c(unique(dat.obj$ddf$species), 
                                                tolower(sort(unique(dat.obj$ddf$common_name))),
                                                unique(dat.obj$ddf$scientific_name),
                                                unique(dat.obj$ddf$common_name)))){

    stop(paste0("Groups do not match available species.\n", "Must be one of: ",
                paste0(c(sort(unique(dat.obj$ddf$species)),
                         sort(unique(dat.obj$ddf$scientific_name)),
                         sort(unique(dat.obj$ddf$common_name)),
                         tolower(sort(unique(dat.obj$ddf$common_name)))), collapse = ", ")))
  }
  
  #' ---------------------------------------------
  # Retrieve species codes
  #' ---------------------------------------------
  
  species.list <- espresso::species_brs %>% 
    janitor::clean_names(.) %>% 
    dplyr::select(scientific_name, common_name, new_code)
  
  sp.groups <- lapply(X = 1:length(species.groups), FUN = function(j) {
    sapply(X = species.groups[[j]], FUN = function(n) {
      index <- which(purrr::map_lgl(
        .x = seq_len(ncol(species.list)),
        .f = ~ any(n %in% c(
          unlist(species.list[, .x]),
          tolower(unlist(species.list[, .x]))
        ))
      ))
      species.list$new_code[which(species.list[, index] == n | tolower(unlist(species.list[, index])) == n)]
    }) %>% unname()
  })
  
  if(is.null(names(species.groups))) {
    names(sp.groups) <- unlist(purrr::map(.x = sp.groups, .f = ~paste0(.x, collapse = ",")))
  } else { names(sp.groups) <- names(species.groups) }
  
  #' ---------------------------------------------
  # Abbreviate group names if necessary
  #' ---------------------------------------------
  
  if(abbrev){
    abbrev.names <- base::abbreviate(names.arg = names(sp.groups), minlength = 5) %>% 
      sub("_$", "", .)
    abbrev.df <- tibble::tibble(name = names(sp.groups), 
                                abbrev = abbrev.names,
                                species = unlist(purrr::map(.x = sp.groups, .f = ~paste0(.x, collapse = ","))))
    names(sp.groups) <- abbrev.names
  }
  
  #' ---------------------------------------------
  # Create groupings
  #' ---------------------------------------------
  
  if(!is.null(species.groups)){
    group.df <- sp.groups %>% tibble::enframe(.) %>% tidyr::unnest(., cols = c(value)) %>% 
      dplyr::rename(species = value, group_name = name)
    brsdat <- dplyr::left_join(x = dat.obj$ddf, y = group.df, by = "species")
    
    # If species.groups contains a subset of all available species
    brsdat <- brsdat %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(group_name = ifelse(is.na(group_name), species, group_name)) %>% 
      dplyr::mutate(species = group_name) %>% 
      dplyr::select(-group_name)
  } else {
    brsdat <- dat.obj$ddf
  }
  
  #' ---------------------------------------------
  # Update dataset and relevant parameters
  #' ---------------------------------------------
  
  n.species <- length(unique(brsdat$species))
  
  species.id <- purrr::map_dbl(.x = unique(brsdat$tag_id),
                               .f = ~{tmp <- brsdat %>% 
                                 dplyr::filter(tag_id == .x)
                               which(unique(unlist(brsdat[, "species"])) == unique(unlist(tmp[, "species"])))})
  
  n.per.species <- as.numeric(table(species.id))
  n.trials <- nrow(brsdat)
  
  # Number of exposures per animal
  # Use factor trick here to conserve order of tag_id
  n.trials.per.whale <- brsdat %>% 
    dplyr::mutate(tag_f = factor(tag_id, levels = unique(tag_id))) %>% 
    dplyr::group_by(tag_f) %>% 
    dplyr::count(.) %>% 
    dplyr::pull(n)
  
  species.trials <- sapply(X = seq_along(species.id), 
                           FUN = function(x) rep(species.id[x], n.trials.per.whale[x])) %>% do.call(c, .)
  
  suppressWarnings(species.summary <- brsdat %>% 
                     dplyr::group_by(species) %>% 
                     dplyr::summarise(N_ind = length(unique(tag_id)), 
                                      N_trials = dplyr::n(), 
                                      censored = sum(is.na(spl)), 
                                      mean = mean(spl, na.rm = TRUE),
                                      min = min(spl, na.rm = TRUE),
                                      max = max(spl, na.rm = TRUE), .groups = "keep") %>% 
                     dplyr::ungroup())
  
  #' ---------------------------------------------
  # Return results
  #' ---------------------------------------------
  
  dat.obj$ddf <- brsdat
  if(!is.null(species.groups)) dat.obj$species$groups <- append(dat.obj$species$groups, sp.groups) else
    dat.obj$species$groups <- NULL
  dat.obj$species$names <- unique(brsdat$species)
  dat.obj$species$n <- n.species
  if(abbrev) dat.obj$species$abbrev <- abbrev.df
  dat.obj$species$nper <- n.per.species
  dat.obj$species$id <- species.id
  dat.obj$species$summary <- species.summary
  dat.obj$species$trials <- species.trials
  
  class(dat.obj) <- c("brsdata.grp", class(dat.obj))
  return(dat.obj)
  
}