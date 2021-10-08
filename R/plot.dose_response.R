#' Dose-response curves
#'
#' Plot \code{dose_response} curve objects produced by \code{\link{compile_rjMCMC}}.
#'
#' @export
#' @param dr.object Input object of class \code{dose_response}, as returned by \code{\link{compile_rjMCMC}}.
#' @param model.rank Integer indicating the rank of the model for which dose-response curves are to be plotted. The model with highest posterior probability is given a rank of 1.
#' @param overlay Logical. If \code{TRUE}, all curves will be overlaid on the same plot.
#' @param do.facet Logical. If \code{TRUE}, plots will be faceted by species (or species group).
#' @param colour.by Curves can be coloured by species group (\code{colour.by = "group"}) or by individual species (\code{colour.by = "species"}), or can use a single colour scheme (the default, \code{colour.by = NULL}).
#' @param colour Used to overwrite the \code{colour.by} argument, if desired.
#' @param colour.median Colour to use for the posterior median. Overwrites the \code{colour.by} argument.
#' @param fill.alpha Opacity applied to the credible intervals. Ranges from 0 (transparent) to 1 (fully opaque).
#' @param order.by How should plots be arranged? Use \code{order.by = "species"} to arrange plots by species (groups) names in alphabetical order, or \code{order.by = "response"}, to arrange plots by response thresholds (from most to least sensitive species (groups), as determined by the estimated posterior medians for μ). Only relevant when \code{overlay = FALSE}.
#' @param show.pmed Logical. If \code{TRUE}, the values of posterior medians are added to facet labels. Only relevant when \code{overlay = FALSE}.
#' @param rotate.y Logical. If \code{TRUE}, y-axis labels are rotated clockwise by 90 degrees.
#' @param all.credint Logical. If \code{TRUE}, plot all credible intervals. If \code{FALSE}, only plot the outermost interval.
#' @param scientific.name Logical. If \code{TRUE}, use species' scientific names in plot titles.
#' @param common.name Logical. If \code{TRUE}, use species' common names in plot titles.
#' @param overwrite.abbrev Logical. If \code{TRUE}, overwrites abbreviated names for species grouped *a priori* using \code{\link{create_groups}}.
#' @param outline.outer Logical. If \code{TRUE}, outline the outermost credible interval bounds.
#' @param do.save Logical. Whether to save plots on disk.
#' @param file.format Character. File format for output plots on disk. Can be one of \code{png}, \code{jpeg}, \code{tiff} or \code{pdf}. Only used when \code{do.save = TRUE}.
#' @param file.dir Path to the output directory. This defaults to the current working directory. Only used when \code{do.save = TRUE}.
#' @param ... Additional arguments as passed to \code{\link[ggplot2]{gg.save}}.
#' 
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{compile_rjMCMC}}
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
# 
#' # Burn and thin
#' rj.trace <- trace_rjMCMC(rj.dat = rj)
#' 
#' # Get dose-response functions
#' doseR <- compile_rjMCMC(rj.trace)
#' 
#' # Plot dose-response curves
#' plot(dr.object = doseR, colour.by = "group")
#' }
#' @keywords brs rjmcmc 

plot.dose_response <- function(dr.object,
                               model.rank = 1,
                               overlay = FALSE,
                               do.facet = TRUE,
                               colour.by = NULL,
                               colour = NULL,
                               colour.median = NULL,
                               fill.alpha = 0.3,
                               order.by = "response", # or "species"
                               show.pmed = TRUE,
                               rotate.y = FALSE,
                               all.credint = FALSE,
                               scientific.name = FALSE,
                               common.name = FALSE,
                               overwrite.abbrev = TRUE,
                               outline.outer = FALSE,
                               do.save = FALSE,
                               file.format = "pdf",
                               ...){
  
  #' ---------------------------------------------
  # Initialisation
  #' ---------------------------------------------
  covariate <- dr.object$covariate
  covariate.values <- dr.object$covariate.values
  species <- dr.object$species
  species.list <- species_brs %>% janitor::clean_names(.)
  fL <- dr.object$fL
  word <- NULL
  
  #' ---------------------------------------------
  # Function checks
  #' ---------------------------------------------
  
  if(!is.null(covariate)) if(covariate == "range") word <- "km"
  if(is.null(colour.median) & all.credint) colour.median <- "black"
  if(!"dose_response" %in% class(dr.object)) stop("dr.object must be of class <dose_response>")
  if(scientific.name & common.name) stop("Only one name allowed.")
  if(!dr.object$by.model) model.rank <- 1
  
  if(!is.null(colour.by)){
    if(!colour.by %in% c("group", "species")) stop("Can only colour by <species> or <group>")}
  
  if(!order.by %in% c("response", "species")) stop("Can only order by <species> or <response>")
  
  if(!is.null(colour.by)){
    if(colour.by == "group"){
    if(!dr.object$by.model){
      colour.by <- "species"
      warning("Curves are coloured by species when by.model = FALSE.")}}}
  
  if(!is.null(colour)) warning("A colour was specified. <colour.by> argument ignored.")
  if(is.null(colour.by)) colour.by <- "species"
  if(dr.object$by.model) colour.by <- "group"
  
  #' ---------------------------------------------
  # Plot colours
  #' ---------------------------------------------
  
  # Number of colours
  if(dr.object$by.model){
    n.colours <- nb_groups(unlist(dr.object$mlist$group[dr.object$mlist$model == dr.object$ranks$model[model.rank]]))
  } else {
  n.colours <- dr.object$dat[[model.rank]]$posterior %>% 
    dplyr::filter(param == "median") %>% 
    nrow()
}
  
  # Generate colour palette
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

  #' ---------------------------------------------
  # Chosen model
  #' ---------------------------------------------
  
  # Retrieve model name and corresponding species groupings
  if(!dr.object$by.model){
    model.rank <- 1
    if(!is.null(covariate)){
      sn <- species
      if(scientific.name) sn <- species.list %>%
          dplyr::filter(new_code == species) %>% dplyr::pull(scientific_name)
      if(common.name) sn <- species.list %>%
          dplyr::filter(new_code == species) %>% dplyr::pull(common_name)
      model.name <- paste0("Species: ", sn) 
      } else {
        model.name <- ""
      }
    
    model.groupings <- 1
    
  } else {
    
    model.name <- dr.object$ranks[model.rank, ]$model
    model.groupings <- dr.object$mlist %>% 
      dplyr::filter(model == model.name) %>% 
      dplyr::pull(group) %>%
      unlist()
    }
  
  #' ---------------------------------------------
  # Set up
  #' ---------------------------------------------
  
  is.sim <- dr.object$sim
  dr.obj <- dr.object$dat
  x.breaks <- pretty(x = dr.obj$dose.range)
  if(all.credint) plot.quants <- dr.obj$cred.int else plot.quants <- max(dr.obj$cred.int)
  
  if(!is.null(covariate)){
    dr.obj[[model.rank]]$posterior <- 
      dr.obj[[model.rank]]$posterior %>% 
      dplyr::mutate(parcov = names(value))}
  
  if(!is.null(covariate.values)){
    dr.obj[[model.rank]]$posterior <- 
      dr.obj[[model.rank]]$posterior %>% 
      dplyr::mutate(parcov = rep(covariate.values, 
                                 ifelse(length(covariate.values) == 1, 
                                                      nrow(dr.obj[[model.rank]]$posterior), 3)))
  }
  
  #' ---------------------------------------------
  # Posterior medians
  #' ---------------------------------------------
  
  median.list <- dr.obj[[model.rank]]$posterior %>% 
    dplyr::filter(param == "median") %>% 
    dplyr::mutate(value = purrr::map(.x = value, .f = ~data.frame(x = dr.obj$dose.range, y = .x))) %>% 
    dplyr::mutate(grp = species)
  
  if(colour.by == "group"){
    median.list <- median.list %>% 
      dplyr::mutate(colgrp = as.character(model.groupings))
  } else {
    median.list <- median.list %>% 
      dplyr::mutate(colgrp = species)
  }

  # Posterior medians for each group
  if(dr.object$by.model) posterior.medians.grp <- dr.object$p.med.bymodel[[model.name]] else posterior.medians.grp <- dr.object$p.med
  
  if(dr.object$phase == 1){
  posterior.medians.grp <- posterior.medians.grp %>% 
    dplyr::filter(stringr::str_detect(param, "mu")) %>% 
    dplyr::rename(grp = param) %>% 
    dplyr::mutate(grp = gsub(pattern = "mu.", replacement = "", grp))
  } else {
    posterior.medians.grp <- posterior.medians.grp %>%
      dplyr::rowwise() %>%
      dplyr::mutate(word.count = sum(c(grepl("nu", param)))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(word.count > 0) %>%
      dplyr::select(-word.count) %>%
      dplyr::rename(grp = param) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(nu = ifelse(grepl("lower", grp), "lower", "upper")) %>% 
      dplyr::mutate(grp = gsub(pattern = "nu.lower.", replacement = "", grp)) %>% 
      dplyr::mutate(grp = gsub(pattern = "nu.upper.", replacement = "", grp)) %>% 
      dplyr::mutate(grp = as.numeric(grp)) %>% 
      tidyr::pivot_wider(., names_from = nu, values_from = pmed) %>% 
      dplyr::arrange(lower) %>% 
      dplyr::mutate(pmed = paste0(round(lower, 2), " / ", round(upper, 2))) %>% 
      dplyr::select(-lower, -upper)
  }
  
  if(!is.null(covariate)){
    
    posterior.medians.grp <- posterior.medians.grp %>% 
      dplyr::filter(grp == which(dr.object$names == dr.object$species))
    
    posterior.medians.covar <- dr.object$p.med %>% 
      dplyr::filter(stringr::str_detect(param, covariate)) %>% 
      dplyr::filter(!stringr::str_detect(param, "incl.")) %>% 
      dplyr::mutate(param = gsub(pattern = paste0(covariate, "_"), replacement = "", param)) %>% 
      dplyr::rename(pmed.cov = pmed)
    
    posterior.medians <- posterior.medians.covar %>% 
      dplyr::mutate(pmed = posterior.medians.grp$pmed) %>% 
      {if(dr.object$phase == 1) dplyr::arrange(., pmed.cov + pmed) else dplyr::arrange(., pmed.cov)} %>% 
      dplyr::mutate(rank = 1:dplyr::n()) %>% 
      dplyr::mutate(grp = dr.object$species) %>% 
      dplyr::rename(parcov = param) %>% 
      tibble::as_tibble() %>% 
      {if(dr.object$phase == 1) . else dplyr::mutate(., grp = as.numeric(grp))}
    
    if(!is.null(covariate.values)) posterior.medians <- posterior.medians %>% 
      dplyr::slice(rep(1:dplyr::n(), each = length(covariate.values))) %>% 
      dplyr::mutate(parcov = covariate.values)
    
  } else {
    
    posterior.medians <- posterior.medians.grp %>% 
      {if(dr.object$phase == 1) dplyr::arrange(., pmed) else .} %>% 
      dplyr::mutate(rank = 1:dplyr::n()) %>% 
      dplyr::mutate(., grp = dr.object$names[as.numeric(grp)])
    
  }
  
  if(!is.null(covariate.values)) posterior.medians <- posterior.medians %>% 
    dplyr::mutate(rank = dplyr::dense_rank(parcov))
  
  #' ---------------------------------------------
  # Credible intervals
  #' ---------------------------------------------
  
  low <- dr.obj[[model.rank]]$posterior %>% 
    dplyr::filter(param == "lower") %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., parcov = covariate.values) else dplyr::mutate(.)} %>% 
    tidyr::unnest(cols = c(value)) %>% 
    dplyr::group_by(species) %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., quant = rep(dr.obj$cred.int, 
                                                                  length(covariate.values)))
      else dplyr::mutate(., quant = rep(dr.obj$cred.int, ifelse(!is.null(covariate), fL$nL, 1)))} %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(lower = value) %>% 
    dplyr::select(-param)
  
  up <- dr.obj[[model.rank]]$posterior %>% 
    dplyr::filter(param == "upper") %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., parcov = covariate.values) else dplyr::mutate(.)} %>% 
    tidyr::unnest(cols = c(value)) %>% 
    dplyr::group_by(species) %>% 
    {if (!is.null(covariate.values)) dplyr::mutate(., quant = rep(dr.obj$cred.int, 
                                                                  length(covariate.values)))
      else dplyr::mutate(., quant = rep(dr.obj$cred.int, ifelse(!is.null(covariate), fL$nL, 1)))} %>% 
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
  
  if(colour.by == "group"){
    interval.list <- interval.list %>% 
      dplyr::mutate(colgrp = rep(as.character(model.groupings), each = length(dr.obj$cred.int)))
  } else {
    interval.list <- interval.list %>% 
      dplyr::mutate(colgrp = species)
  }
  
  if(!is.null(covariate)){
    median.list <- median.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = c("parcov", "grp"))
    interval.list <- interval.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = c("parcov", "grp"))
  } else {
    median.list <- median.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = "grp")
    interval.list <- interval.list %>% 
      dplyr::left_join(x = ., y = posterior.medians, by = "grp")
    }
  
  #' ---------------------------------------------
  # Species labels
  #' ---------------------------------------------
  
  if(!is.sim){
    
    which.name <- c("scientific_name", "common_name")[which(c(scientific.name, common.name))]
    
    if(length(which.name) == 0){
      
      if(!is.null(dr.object$abbrev[[1]])){
        
        if(overwrite.abbrev){
          
          median.list <- dr.object$abbrev %>% 
            dplyr::rename(label = species) %>% 
            dplyr::select(abbrev, label) %>% 
            dplyr::left_join(x = median.list,
                             y = ., 
                             by = c("species" = "abbrev")) %>% 
            dplyr::mutate(label = ifelse(is.na(label), colgrp, label))
          
          interval.list <- dr.object$abbrev %>% 
            dplyr::rename(label = species) %>% 
            dplyr::select(abbrev, label) %>% 
            dplyr::left_join(x = interval.list,
                             y = ., 
                             by = c("species" = "abbrev")) %>% 
            dplyr::mutate(label = ifelse(is.na(label), colgrp, label))
          
        } else {
          
          median.list$label <- median.list$species
          interval.list$label <- interval.list$species
          
        }
        
        
      } else {
        
        median.list$label <- median.list$species
        interval.list$label <- interval.list$species
        
      }
      
    } else {
      
      species.list <- species.list %>% 
        dplyr::filter(new_code %in% unique(c(dr.obj$output$posterior$species,
                                             unlist(dr.object$sp.groups))))
      
      species.list <- species.list[, c(which.name, "new_code")]
      names(species.list)[which(names(species.list) == which.name)] <- "label"
      
      if(!is.null(dr.object$sp.groups)){
        for(nn in seq_along(dr.object$sp.groups)){
          grp.name <- names(dr.object$sp.groups)[nn]
          grp.label <- species.list %>% 
            dplyr::filter(new_code %in% dr.object$sp.groups[[nn]]) %>% 
            dplyr::arrange(label) %>% 
            dplyr::pull(label) %>% 
            paste0(., collapse = " / ")
        
          species.list <- species.list %>% 
            dplyr::filter(!new_code %in% dr.object$sp.groups[[nn]]) %>% 
            dplyr::bind_rows(., tibble::tibble(label = grp.label, new_code = grp.name))
        }
      }
      
      median.list <- dplyr::left_join(x = median.list, y = species.list, by = c("species" = "new_code"))
      interval.list <- dplyr::left_join(x = interval.list, y = species.list, by = c("species" = "new_code"))
    }
    
  } else {
    
    median.list$label <- median.list$species
    interval.list$label <- interval.list$species
    
  }
  
  if(dr.object$by.model){
    
  median.list$grp <- purrr::map(.x = median.list$colgrp,
               .f = ~median.list %>% dplyr::filter(colgrp %in% .x) %>% dplyr::pull(grp)) %>% 
      purrr::map(.x = ., .f = ~paste0("(", paste0(., collapse = ","), ")")) %>% 
    unlist()
  
  interval.list <- interval.list %>% 
    dplyr::select(-grp) %>% 
    dplyr::left_join(., y = median.list[, c("grp", "colgrp")], by = "colgrp")

  }
  
  #' ---------------------------------------------
  # Labels for posterior medians
  #' ---------------------------------------------
  
  if(show.pmed){
    
    if(!is.null(covariate)){
      
      median.list <- median.list %>% 
        dplyr::rowwise() %>% 
        {if(dr.object$phase == 1) 
        dplyr::mutate(., colgrp = paste0(parcov, " ", word, " (", 
        round(pmed + ifelse(!is.null(covariate.values), pmed.cov * parcov, pmed.cov), 1), " dB)"))
          else dplyr::mutate(., colgrp = paste(parcov, word)) } %>% 
        dplyr::ungroup()
      
      interval.list <- interval.list %>% 
        dplyr::rowwise() %>% 
        {if(dr.object$phase == 1) 
        dplyr::mutate(., colgrp = paste0(parcov, " ", word, " (", round(pmed + ifelse(!is.null(covariate.values), pmed.cov * parcov, pmed.cov), 1), " dB)")) 
          else dplyr::mutate(., colgrp = paste(parcov, word)) } %>% 
        dplyr::ungroup()
      
    } else {
      
      median.list <- median.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(grp = paste0(ifelse(is.sim, "Species ", " "),
                                   ifelse(is.null(covariate), grp, label),
                                   " (", ifelse(is.numeric(pmed), round(pmed, 1), pmed), " dB)"))  %>% 
        dplyr::ungroup()
      
      interval.list <- interval.list %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(grp = paste0(ifelse(is.sim, "Species ", " "),
                                   ifelse(is.null(covariate), grp, label),
                                   " (", ifelse(is.numeric(pmed), round(pmed, 1), pmed), " dB)")) %>% 
        dplyr::ungroup()
    }
    
  } else {
    
    if(!is.null(covariate)){
      median.list <- median.list %>% dplyr::mutate(colgrp = parcov)
      interval.list <- interval.list %>% dplyr::mutate(colgrp = parcov)
    } else {
      if(is.sim){
      median.list <- median.list %>% dplyr::mutate(grp = label)
      interval.list <- interval.list %>% dplyr::mutate(grp = label)
      }
    }
    
    
  }
  
  if(order.by == "response") {
    
    median.list <- median.list %>% 
      dplyr::mutate(grp = factor(grp)) %>% 
      dplyr::mutate(grp = forcats::fct_reorder(grp, rank))
    
    interval.list <- interval.list %>% 
      dplyr::mutate(grp = factor(grp)) %>%
      dplyr::mutate(grp = forcats::fct_reorder(grp, rank))
    
  }
  
  if(!is.null(covariate.values) | !is.null(covariate)) {
    
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
    
    gdat <- interval.list %>% 
      dplyr::filter(quant == plot.quants[qq]) %>% 
      tidyr::unnest(cols = c(value))
    
    if(!is.null(covariate)){
      
      gplot <- gplot +
        geom_polygon(data = gdat, aes(x = x, y = y, group = colgrp), fill = "white") +
        geom_polygon(data = gdat, aes(x = x, y = y, group = colgrp, fill = colgrp), 
                     alpha = seq(fill.alpha, 1, length = length(plot.quants))[qq]) 
    } else {
      
      gplot <- gplot +
        geom_polygon(data = gdat, aes(x = x, y = y), fill = "white") +
        {if(colour.by == "group" & is.null(colour)) geom_polygon(data = gdat,
                                               aes(x = x, y = y, fill = colgrp),
                                               alpha = seq(fill.alpha, 1, length = length(plot.quants))[qq]) } + 
        {if(colour.by == "species" & is.null(colour)) geom_polygon(data = gdat,
                                                 aes(x = x, y = y, fill = grp),
                                                 alpha = seq(fill.alpha, 1, length = length(plot.quants))[qq]) } +
        {if(!is.null(colour)) geom_polygon(data = gdat, aes(x = x, y = y), fill = colour,
                                                 alpha = seq(fill.alpha, 1, length = length(plot.quants))[qq]) }
    }
    
    if(outline.outer){
      
      if(!is.null(covariate)){
        gplot <- gplot +
          geom_polygon(data = gdat %>% dplyr::filter(quant == max(plot.quants)),
                       aes(x = x, y = y, colour = colgrp), fill = "transparent", linetype = "dashed") +
          {if(!overlay & do.facet) facet_wrap(~ colgrp) }
      } else {
        if(all.credint){
          gplot <- gplot +
            geom_polygon(data = gdat %>% dplyr::filter(quant == max(plot.quants)),
                         aes(x = x, y = y), colour = "black", fill = "transparent", linetype = "dashed")
        }else{
          gplot <- gplot +
            geom_polygon(data = gdat %>% dplyr::filter(quant == max(plot.quants)),
                         aes(x = x, y = y, colour = grp), fill = "transparent", linetype = "dashed")
        }
        # Facet plots
        if(!overlay & do.facet) { gplot <- gplot + facet_wrap(~ grp) }
        
        # Output plots separately
        if(!overlay & !do.facet){ gplot <- gplot + ggforce::facet_wrap_paginate(~ grp, ncol = 1, nrow = 1) }
       
      }
    }
    
  } 
  
  mdat <- median.list %>% tidyr::unnest(cols = c(value))
  
  if(!is.null(covariate)){
    
    gplot <- gplot +
      geom_line(data = mdat, aes(x = x, y = y, group = parcov, colour = colgrp)) +
      {if(!is.null(colour.median)) geom_line(data = mdat, 
                                             aes(x = x, y = y, group = colgrp), 
                                             colour = colour.median) } +
      { if(!overlay & do.facet) facet_wrap(~ colgrp) }
    
  } else {
    
    gplot <- gplot +
      {if(colour.by == "group") geom_line(data = mdat, aes(x = x, y = y, colour = colgrp)) } +
      {if(colour.by == "species") geom_line(data = mdat, aes(x = x, y = y, colour = grp)) } +
      {if(!is.null(colour.median)) geom_line(data = mdat, aes(x = x, y = y), colour = colour.median) } +
      {if(!overlay & do.facet) facet_wrap(~ grp) } +
      {if(!overlay & !do.facet) ggforce::facet_wrap_paginate(~ grp, nrow = 1, ncol = 1) }
  }
  
  gplot <- gplot +
    scale_fill_manual(values = mycols) +
    scale_colour_manual(values = mycols) +
    ylab("p(response)") + 
    xlab(expression(paste("Dose (dB re 1", mu, "Pa)")))
  
  
  gplot <- gplot + 
    ggtitle(ifelse(!is.null(covariate), "Dose-response model", "Multi-species dose-response model"),
            subtitle = model.name) +
    theme(plot.title = element_text(size = 12, vjust = 5, face = "bold"),
          plot.subtitle = element_text(size = 11, vjust = 5))
  
  if(!overlay & !do.facet){
    for(page in 1:ggforce::n_pages(gplot)){
      p_save <- gplot +
        facet_wrap_paginate(~ grp, ncol = 1, nrow = 1, page = page)
      
      if(do.save){
        suppressMessages(ggplot2::ggsave(filename = paste0("dose_response (page ", page, ").", file.format), 
                                         device = file.format, plot = p_save, ...))
      } else {
        print(p_save) 
      }
      
    }
  } else {
    if(do.save){
      suppressMessages(ggplot2::ggsave(filename = paste0("dose_response.", file.format), 
                                       device = file.format, plot = gplot, ...))
    } else {
      print(gplot) 
      }
  }
  
}