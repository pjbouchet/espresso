#' Dose-response curves
#'
#' Plot method for \code{dose_response} objects produced by \code{\link{compile_rjMCMC}}.
#'
#' @param dr.object Input object of class \code{dose_response}.
#' @param which.model Integer indicating the ID of the model for which dose-response curves are to be plotted. The \code{\link{summary.dose_response}} function can be used to list model IDs, if necessary.
#' @param overlay Logical. If \code{TRUE}, plots will be faceted by species (or species group). If \code{FALSE}, all curves will be overlaid on the same plot.
#' @param colour.by By default, all dose-response curves follow the same colour scheme (\code{colour.by = NULL}). However, curves can also be coloured by \code{group} or by \code{species}.
#' @param colour Used to overwrite the \code{colour.by} argument, if desired.
#' @param colour.median Colour to use for the posterior median line. Overwrites the \code{colour.by} argument when specified.
#' @param order.by How should plots be ordered? When \code{order.by = "species"}, plots are ordered alphabetically. \code{order.by = "response"}, plots are ordered by response thresholds (from most to least sensitive species, as determined by posterior medians for μ. Only relevant when \code{overlay = FALSE}.
#' @param show.p0_5 Logical. If \code{TRUE}, posterior medians are added to facet labels. Only relevant when \code{overlay = FALSE}.
#' @param rotate.y Logical. If \code{TRUE}, y-axis labels are rotated clockwise by 90 degrees.
#' @param all.quants 
#' @param scientific.name Logical. If \code{TRUE}, use species' scientific names in plot titles.
#' @param common.name Logical. If \code{TRUE}, use species' common names in plot titles.
#' @param overwrite.abbrev Logical. If \code{TRUE}, overwrites abbreviated names for species grouped a priori using \code{\link{cerate_groups}}.
#' @param outline.95 Logical. If \code{TRUE}, outlines the 95% posterior credible interval bounds.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.rjdata}}
#' @examples
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
#' 
#' @export
#' @keywords brs rjmcmc 

plot.dose_response <- function(dr.object,
                               which.model = 1,
                               overlay = FALSE,
                               colour.by = NULL,
                               colour = NULL,
                               colour.median = NULL,
                               order.by = "species", # or "response"
                               show.p0_5 = TRUE,
                               rotate.y = FALSE,
                               all.quants = FALSE,
                               scientific.name = TRUE,
                               common.name = FALSE,
                               overwrite.abbrev = TRUE,
                               outline.95 = FALSE){
  
  covariate <- dr.object$covariate
  covariate.values <- dr.object$covariate.values
  species <- dr.object$species
  
  word <- NULL
  if(!is.null(covariate)) if(covariate == "range") word <- "km"
  
  if(!"dose_response" %in% class(dr.object)) stop("dr.object must be of class <dose_response>")
  if(scientific.name & common.name) stop("Only one name allowed.")
  
  if(!dr.object$by.model) which.model <- 1
  
  if(!is.null(colour.by)){
    if(!colour.by %in% c("group", "species")) stop("Can only colour by <species> or <group>")}
  
  if(!order.by %in% c("response", "species")) stop("Can only order by <species> or <response>")
  
  if(!is.null(colour.by)){
    if(!dr.object$by.model){
      colour.by <- "species"
      warning("No group-level information available. Dose-response curves will be coloured by species.")}}
  
  if(!is.null(colour)) warning("A colour was specified. <colour.by> argument ignored.")
  
  n.colours <- dr.object$dat[[which.model]]$posterior %>% 
    dplyr::filter(param == "median") %>% 
    nrow()
  
  if (is.null(colour)) { # No custom colour specified
    if (is.null(colour.by)) { # No colouring by group or species (all same colour)
      mycols <- rep("#0098a0", n.colours)
    } else {
      if (n.colours == 1) {
        mycols <- "#0098a0"
      } else if (n.colours == 2) {
        mycols <- c("#0098a0", "#E0B200")
      } else if (n.colours > 2) {
        mycols <- pals::viridis(n = n.colours)
      }
    }
  } else {
    if (is.null(colour.by)) {
      mycols <- rep(colour, n.colours)
    } else {
      mycols <- rep("#0098a0", n.colours)
    }
  }
  
  if(is.null(colour.by)) colour.by <- "species"
  
  if(!is.null(dr.object$sp.groups) & n.colours > 1){
    scientific.name <- common.name <- FALSE
  }
  
  if(!dr.object$by.model){
    which.model <- 1
    if(!is.null(covariate)) model.name <- paste0("Species: ", dr.object$species) else model.name <- ""
    model.groupings <- 1
    
  } else {
    
    model.name <- names(dr.object$dat)[which.model]
    model.groupings <- dr.object$mlist %>% 
      dplyr::filter(model == model.name) %>% 
      dplyr::pull(group) %>% unlist()}
  
  is.sim <- dr.object$sim
  
  dr.obj <- dr.object$dat
  x.breaks <- pretty(x = dr.obj$dose.range)
  
  if(all.quants) plot.quants <- dr.obj$quants else plot.quants <- 95
  
  if(!is.null(covariate)){
    dr.obj[[which.model]]$posterior <- 
      dr.obj[[which.model]]$posterior %>% 
      dplyr::mutate(parcov = names(value))}
  
  
  if(!is.null(covariate.values)){
    dr.obj[[which.model]]$posterior <- 
      dr.obj[[which.model]]$posterior %>% 
      dplyr::mutate(parcov = rep(covariate.values, length(covariate.values)))
  }
  
  median.list <- dr.obj[[which.model]]$posterior %>% 
    dplyr::filter(param == "median") %>% 
    dplyr::mutate(value = purrr::map(.x = value, .f = ~data.frame(x = dr.obj$dose.range, y = .x))) %>% 
    dplyr::mutate(grp = species)
  
  if(colour.by == "group") median.list$colgrp <- as.character(model.groupings[species]) else median.list$colgrp <- median.list$species
  
  posterior.medians.grp <- dr.object$p.med %>% 
    dplyr::filter(stringr::str_detect(param, "mu")) %>% 
    dplyr::rename(grp = param) %>% 
    dplyr::mutate(grp = gsub(pattern = "mu.", replacement = "", grp))
  
  if(!is.null(covariate)){
    
    posterior.medians.grp <- posterior.medians.grp %>% 
      dplyr::filter(grp == dr.object$species)
    
    posterior.medians.covar <- dr.object$p.med %>% 
      dplyr::filter(stringr::str_detect(param, covariate)) %>% 
      dplyr::filter(!stringr::str_detect(param, "incl.")) %>% 
      dplyr::mutate(param = gsub(pattern = paste0(covariate, "_"), replacement = "", param)) %>% 
      dplyr::rename(pmed.cov = pmed)
    
    posterior.medians <- posterior.medians.covar %>% 
      dplyr::mutate(pmed = posterior.medians.grp$pmed) %>% 
      dplyr::arrange(pmed.cov + pmed) %>% 
      dplyr::mutate(rank = 1:dplyr::n()) %>% 
      dplyr::mutate(grp = dr.object$species) %>% 
      dplyr::rename(parcov = param)
    
    if(!is.null(covariate.values)) posterior.medians <- posterior.medians %>% 
      slice(rep(1:n(), each = length(covariate.values))) %>% 
      dplyr::mutate(parcov = covariate.values)
    
  } else {
    
    posterior.medians <- posterior.medians.grp %>% 
      dplyr::arrange(pmed) %>% 
      dplyr::mutate(rank = 1:dplyr::n())
  }
  
  if(!is.null(covariate.values)) posterior.medians <- posterior.medians %>% 
    dplyr::mutate(rank = dplyr::dense_rank(parcov))
  
  # is.sim <- is.numeric(suppressWarnings(as.numeric(posterior.medians$grp)))
  
  low <- dr.obj[[which.model]]$posterior %>% 
    dplyr::filter(param == "lower") %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., parcov = covariate.values) else dplyr::mutate(.)} %>% 
    tidyr::unnest(cols = c(value)) %>% 
    dplyr::group_by(species) %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., quant = rep(dr.obj$quants, 
                                                                  length(covariate.values)))
      else dplyr::mutate(., quant = rep(dr.obj$quants, 
                                        ifelse(!is.null(covariate), 
                                               fL[[covariate]]$nL, 1)))} %>% 
    
    dplyr::ungroup() %>% 
    dplyr::rename(lower = value) %>% 
    dplyr::select(-param)
  
  up <- dr.obj[[which.model]]$posterior %>% 
    dplyr::filter(param == "upper") %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., parcov = covariate.values) else dplyr::mutate(.)} %>% 
    tidyr::unnest(cols = c(value)) %>% 
    dplyr::group_by(species) %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., quant = rep(dr.obj$quants, 
                                                                  length(covariate.values)))
      else dplyr::mutate(., quant = rep(dr.obj$quants, 
                                        ifelse(!is.null(covariate), 
                                               fL[[covariate]]$nL, 1)))} %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(upper = value) %>% 
    dplyr::select(-param)
  
  if(!is.null(covariate)){
    interval.list <- dplyr::left_join(x = low, y = up, by = c("species", "quant", "parcov"))
  } else {
    interval.list <- dplyr::left_join(x = low, y = up, by = c("species", "quant"))
  }
  
  interval.list <- interval.list %>% 
    dplyr::mutate(value = purrr::map2(.x = lower, 
                                      .y = upper,
                                      .f = ~data.frame(x = c(dr.obj$dose.range, rev(dr.obj$dose.range)),
                                                       y = c(.x, rev(.y))))) %>% 
    dplyr::mutate(grp = species)
  
  if(colour.by == "group") interval.list$colgrp <- as.character(model.groupings[species]) else interval.list$colgrp <- interval.list$species
  
  
  if(!is.null(covariate)){
    median.list <- median.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = c("parcov", "grp"))
    interval.list <- interval.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = c("parcov", "grp"))
  } else {
    median.list <- median.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, c("colgrp" = "grp"))
    interval.list <- interval.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, c("colgrp" = "grp"))}
  
  
  if(!is.sim){
    
    which.name <- c("scientific_name", "common_name")[which(c(scientific.name, common.name))]
    
    if(length(which.name) == 0){
      
      if(!is.null(dr.object$abbrev)){
        
        if(overwrite.abbrev){
          
          median.list <- dr.object$abbrev %>% 
            dplyr::rename(label = species) %>% 
            dplyr::select(abbrev, label) %>% 
            dplyr::left_join(x = median.list,
                             y = ., 
                             by = c("species" = "abbrev"))
          
          interval.list <- dr.object$abbrev %>% 
            dplyr::rename(label = species) %>% 
            dplyr::select(abbrev, label) %>% 
            dplyr::left_join(x = interval.list,
                             y = ., 
                             by = c("species" = "abbrev"))
          
        } else {
          
          median.list$label <- median.list$species
          interval.list$label <- interval.list$species
          
        }
        
        
      } else {
        
        median.list$label <- median.list$species
        interval.list$label <- interval.list$species
        
      }
      
    } else {
      
      species.list <- espresso::species_brs %>% 
        janitor::clean_names(.) %>% 
        dplyr::filter(new_code %in% unique(dr.obj$output$posterior$species))
      
      species.list <- species.list[, c(which.name, "new_code")]
      names(species.list)[which(names(species.list) == which.name)] <- "label"
      
      median.list <- dplyr::left_join(x = median.list, y = species.list, by = c("species" = "new_code"))
      interval.list <- dplyr::left_join(x = interval.list, y = species.list, by = c("species" = "new_code"))
    }
    
  }
  
  if(show.p0_5){
    
    if(!is.null(covariate)){
      
      median.list <- median.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(colgrp = paste0(parcov, " ", word, " (", 
                                      round(pmed + ifelse(!is.null(covariate.values), pmed.cov * parcov, pmed.cov), 1), " dB)")) %>% 
        dplyr::ungroup()
      
      interval.list <- interval.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(colgrp = paste0(parcov, " ", word, " (", round(pmed + ifelse(!is.null(covariate.values), pmed.cov * parcov, pmed.cov), 1), " dB)")) %>% 
        dplyr::ungroup()
      
    } else {
      
      median.list <- median.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(grp = paste0(ifelse(is.sim, "Species ", ""),
                                   ifelse(!is.null(covariate), species, label),
                                   " (", round(pmed, 1), " dB)"))  %>% 
        dplyr::ungroup()
      
      interval.list <- interval.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(grp = paste0(ifelse(is.sim, "Species ", ""),
                                   ifelse(!is.null(covariate), species, label),
                                   " (", round(pmed, 1), " dB)")) %>% 
        dplyr::ungroup()
    }
    
  } else {
    
    median.list <- median.list %>% dplyr::mutate(colgrp = parcov)
    interval.list <- interval.list %>% dplyr::mutate(colgrp = parcov)
  }
  
  if(order.by == "response") {
    
    median.list <- median.list %>% 
      dplyr::mutate(grp = factor(grp)) %>% 
      dplyr::mutate(grp = forcats::fct_reorder(grp, rank))
    
    interval.list <- interval.list %>% 
      dplyr::mutate(grp = factor(grp)) %>%
      dplyr::mutate(grp = forcats::fct_reorder(grp, rank))
    
  }
  
  if(!is.null(covariate.values)) {
    
    median.list <- median.list %>% 
      dplyr::mutate(colgrp = factor(colgrp)) %>% 
      dplyr::mutate(colgrp = forcats::fct_reorder(colgrp, rank))
    
    interval.list <- interval.list %>% 
      dplyr::mutate(colgrp = factor(colgrp)) %>%
      dplyr::mutate(colgrp = forcats::fct_reorder(colgrp, rank))
    
  }
  
  gplot <- ggplot() +
    theme(axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(angle = ifelse(rotate.y, 90, 0), hjust = 0.5, vjust = 0.5),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 00, l = 0)),
          plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, "cm"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) + 
    {if(!overlay) theme(legend.position = "none", strip.text.x = element_text(size = 11))}
  
  
  for(qq in seq_along(plot.quants)){
    
    gdat <- interval.list %>% dplyr::filter(quant == plot.quants[qq]) %>% 
      tidyr::unnest(cols = c(value))
    
    if(!is.null(covariate)){
      
      gplot <- gplot +
        geom_polygon(data = gdat, aes(x = x, y = y, group = colgrp), fill = "white") +
        geom_polygon(data = gdat, aes(x = x, y = y, group = colgrp, fill = colgrp), 
                     alpha = seq(0.3, 1, length = length(plot.quants))[qq]) 
    } else {
      
      gplot <- gplot +
        geom_polygon(data = gdat, aes(x = x, y = y), fill = "white") +
        {if(colour.by == "group") geom_polygon(data = gdat,
                                               aes(x = x, y = y, fill = colgrp),
                                               alpha = seq(0.3, 1, length = length(plot.quants))[qq]) } + 
        {if(colour.by == "species") geom_polygon(data = gdat,
                                                 aes(x = x, y = y, fill = grp),
                                                 alpha = seq(0.3, 1, length = length(plot.quants))[qq]) } 
    }
  } 
  
  mdat <- median.list %>% tidyr::unnest(cols = c(value))
  
  if(!is.null(covariate)){
    
    gplot <- gplot +
      geom_line(data = mdat, aes(x = x, y = y, group = parcov, colour = colgrp)) +
      {if(!is.null(colour.median)) geom_line(data = mdat, 
                                             aes(x = x, y = y, group = colgrp), 
                                             colour = colour.median) } +
      { if(!overlay) facet_wrap(~ colgrp) }
    
  } else {
    gplot <- gplot +
      {if(colour.by == "group") geom_line(data = mdat, aes(x = x, y = y, colour = colgrp)) } +
      {if(colour.by == "species") geom_line(data = mdat, aes(x = x, y = y, colour = grp)) } +
      {if(!is.null(colour.median)) geom_line(data = mdat, aes(x = x, y = y), colour = colour.median) } +
      {if(!overlay) facet_wrap(~ grp) }
  }
  
  gplot <- gplot +
    scale_fill_manual(values = mycols) +
    scale_colour_manual(values = mycols) +
    ylab("p(response)") + 
    xlab(expression(paste("Dose (dB re 1", mu, "Pa)")))
  
  if(outline.95){
    
    if(!is.null(covariate)){
      gplot <- gplot +
        geom_polygon(data = gdat %>% dplyr::filter(quant == 95),
                     aes(x = x, y = y, colour = colgrp), fill = "transparent", linetype = "dashed") +
        {if(!overlay) facet_wrap(~ colgrp) }
    } else {
      gplot <- gplot +
        geom_polygon(data = gdat %>% dplyr::filter(quant == 95),
                     aes(x = x, y = y, colour = grp), fill = "transparent", linetype = "dashed") +
        {if(!overlay) facet_wrap(~ grp) }
    }
  }
  
  gplot <- gplot + 
    ggtitle(ifelse(!is.null(covariate), "Dose-response model", "Multi-species dose-response model"),
            subtitle = model.name) +
    theme(plot.title = element_text(size = 12, vjust = 5, face = "bold"),
          plot.subtitle = element_text(size = 11, vjust = 5))
  
  print(gplot)
  
}