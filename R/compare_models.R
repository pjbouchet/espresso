#' Benchmark results
#'
#' Evaluate outputs from models implemented using a (1) rjMCMC and (2) Gibbs Variable Selection approach. Returns comparative plots of posterior model rankings and posterior parameter estimates.
#'
#' @export
#' @param rj.dat Input rjMCMC object, as returned by \code{\link{trace_rjMCMC}}.
#' @param gvs.fixed Input Gibbs object, as returned by \code{\link{gibbs}} using a fixed effect implementation of the dose-response model.
#' @param gvs.random Input Gibbs object, as returned by \code{\link{gibbs}} using a random effect implementation of the dose-response model.
#' @param by.model Logical. If \code{TRUE}, the functions subsets posterior estimates by candidate model.
#' @param kernel.adj Bandwidth adjustment. The bandwidth used to create density plots is given by \code{kernel.adj*bw}.
#' @param viridis.col Logical. Whether to use a viridis colour scheme.
#' @param density Logical. If \code{TRUE}, compares density plots for each model parameter.
#' @param prob Logical. If \code{TRUE}, compares posterior rankings for candidate models.
#' 
#' @author Phil J. Bouchet
#' @seealso \code{\link{simulate_data}} \code{\link{example_brs}} \code{\link{summary.rjdata}}
#' @examples
#' \dontrun{
#' library(espresso)
#' 
#' # Simulate data for two species
#' mydat <- simulate_data(n.species = 2, 
#'                        n.whales = 16, 
#'                        max.trials = 3, 
#'                        covariates = list(exposed = c(0, 5), range = 0.5),
#'                        mu = c(101, 158), 
#'                        phi = 20, 
#'                        sigma = 20, 
#'                        Rc = c(210, 211), 
#'                        seed = 58697)
#' summary(mydat)
#'      
#' # Model selection by GVS                        
#' gvs <- gibbs(dat = mydat, 
#'              random.effects = FALSE, 
#'              include.covariates = FALSE, 
#'              mcmc.n = 1000, 
#'              burnin = 500)
#' 
#' # Run the reversible jump MCMC
#' rj <- run_rjMCMC(dat = mydat.config,
#'                  n.chains = 2,
#'                  n.burn = 100,
#'                  n.iter = 100,
#'                  do.update = FALSE)
#'                  
#' # Burn and thin
#' rj.trace <- trace_rjMCMC(rj.dat = rj)
#' 
#' # Compare outputs
#' compare_models(rj.dat = rj.trace, gvs.fixed = gvs)
#' }
#' @keywords brs rjmcmc gvs dose-response

compare_models <- function(rj.dat = NULL,
                           gvs.fixed = NULL,
                           gvs.random = NULL,
                           by.model = FALSE,
                           kernel.adj = 2,
                           viridis.col = FALSE,
                           density = TRUE,
                           prob = TRUE
){
  
  if(!"rjtrace" %in% class(rj.dat)) stop("Input data must be of class <rjtrace>.")
  if(!"gvs" %in% c(class(gvs.fixed), class(gvs.random))) stop("Input data must be of class <gvs>.")
  if(!is.null(gvs.fixed)){
    if(!gvs.fixed$type == "fixed") stop("Erroneous input for <gvs.fixed>.")}
  if(!is.null(gvs.random)){
    if(!gvs.random$type == "random") stop("Erroneous input for <gvs.random>.")}
  
  if(by.model){
    print(data.frame(rj.dat$mlist))
    plot.model <- as.numeric(readline(prompt = "Which model should be plotted? "))
    if(!is.numeric(plot.model)) stop("Unrecognised model.")
    which.model <- rj.dat$mlist$model[plot.model]
  } else {
    plot.model <- NULL
  }
  
  master.list <- list(rj = rj.dat, gvs_fixed = gvs.fixed, gvs_random = gvs.random)
  master.L <- purrr::compact(.x = master.list) %>% length
  if(master.L <= 1) stop("Two or more elements required.")
  
  m.methods <- tibble::tibble(method = c("rjMCMC", "GVS (fixed effects)", "GVS (random effects)"),
                              mcmc = c("rj", "gvs_fixed", "gvs_random"))
  legend.methods <- m.methods$method[!sapply(X = master.list, FUN = is.null)]
  
  # Colours for plotting
  if(viridis.col) plot.colours <- c("#2C708E", "#29AF7F", "#ffce0a") else plot.colours <- c("#00AFBB", "#E7B800", "#FC4E09")
  
  plot.colours <- plot.colours[!sapply(X = master.list, FUN = is.null)]
  master.list <- purrr::compact(master.list)
  
  
  if(prob){
    
    combined.res <- purrr::map2(.x = master.list, 
                                .y = names(master.list),
                                .f = ~prob_models(input.obj = .x, 
                                                  mlist = .x$mlist, 
                                                  select = .x$config$model.select,
                                                  do.combine = TRUE,
                                                  gvs = grepl(pattern = "gvs", .y))) %>% 
      tibble::enframe(.) %>% 
      dplyr::rename(mcmc = name) %>% 
      dplyr::left_join(., m.methods, by = "mcmc") %>% 
      dplyr::mutate(value = purrr::map(.x = value, 
                                       .f = ~{
                                         .x$model$m_prob %>% 
                                           dplyr::arrange(p) %>%
                                           dplyr::mutate(model = factor(model, levels = model))})) %>% 
      tidyr::unnest(cols = c(value))
    
    
    gg.opts <- theme(axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     legend.position = "top")
    
    gPlot <- ggplot2::ggplot(data = combined.res, aes(x = model, y = p, colour = method)) +
      ggplot2::geom_pointrange(aes(ymin = lower, ymax = upper), 
                               position = position_dodge(width = 0.2), 
                               size = 0.5, shape = 16) + 
      ggplot2::geom_point(aes(group = mcmc), 
                          colour = "black", size = 2, shape = 21, 
                          position = position_dodge(width = 0.2)) + 
      ggplot2::scale_color_manual(values = plot.colours, name = "") +
      coord_flip() +
      ylab("Posterior probability") + xlab("") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      gg.opts + theme(plot.margin = unit(c(0.5, 2, 0.5, 0.5), "cm"),
                      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 20, l = 0)))
    
    print(gPlot)
  }
  
  
  if(density){
    
    # Split trace by model ID
    master.trace <- purrr::map(.x = master.list, 
                               .f = ~{.x$trace %>% 
                                   do.call(rbind, .) %>% 
                                   tibble::as_tibble(.) %>% 
                                   dplyr::rename_at(vars(matches("theta")), function(x) "model_ID") %>% 
                                   dplyr::select(-contains("incl.")) %>% 
                                   dplyr::select(-contains("size"))})
    
    master.trace$rj <- master.trace$rj %>% 
      dplyr::left_join(x = ., y = master.list$rj$mlist[, c("ID", "model")], by = c("model_ID" = "ID"))
    
    if(!is.null(gvs.fixed)){
      master.trace$gvs_fixed <- master.trace$gvs_fixed %>% 
        dplyr::left_join(x = ., 
                         y = tibble::tibble(ID = seq_along(master.list$gvs_fixed$species.Groups),
                                            model = master.list$gvs_fixed$species.Groups), 
                         by = c("model_ID" = "ID"))
      
      names(master.trace$gvs_fixed) <- 
        gsub(pattern = "_i\\[", replacement = ".", x = names(master.trace$gvs_fixed)) %>% 
        gsub(pattern = "\\]", replacement = "", x =.)
      
    }
    
    nn <- purrr::map(.x = master.trace, .f = ~colnames(as.data.frame(.x)))
    nn.check <- purrr::map(.x = nn, .f = ~paste0(sort(.x), collapse = "+")) %>% unlist(.) %>% unique(.)
    if(length(nn.check) > 1) stop("Column names must match!") else nn <- unlist(nn) %>% unique(.)
    
    # par(mfrow = c(ifelse(round(length(nn) / 3) == 1, 2, 1), 
    #               ifelse(round(length(nn) / 3) == 1, 3, 4)))
    
    master.plot <- purrr::map(.x = 1:length(master.trace),
                              .f = ~ dplyr::mutate(master.trace[[.x]],
                                                   method = legend.methods[.x])) %>% 
      do.call(rbind, .) %>% 
      dplyr::select(-model_ID)
    
    # When by.model, only makes senses to subset mu as other params are not grouping dependent
    mu.plot <- master.plot %>% 
      dplyr::select_at(vars(contains("mu"), "model", "method"))
    
    master.plot <- master.plot %>% 
      dplyr::select_at(vars(-contains("mu")))
    
    if(!is.null(plot.model))  mu.plot <- dplyr::filter(mu.plot, model == which.model)
    
    mu.plot <- mu.plot %>% 
      dplyr::select(., -model) %>% 
      tidyr::pivot_longer(!method, names_to = "param", values_to = "value") 
    
    master.plot <- master.plot %>% 
      dplyr::select(., -model) %>% 
      tidyr::pivot_longer(!method, names_to = "param", values_to = "value") 
    
    master.plot <- dplyr::bind_rows(master.plot, mu.plot)
    
    gPlot2 <- ggplot(master.plot, aes(x = value, colour = method))  +
      geom_density(adjust = kernel.adj, show.legend = FALSE) +
      stat_density(geom = "line", position = "identity", size = 0) +
      facet_wrap(~ param, scales = "free") +
      scale_colour_manual(values = plot.colours) +
      xlab("") + ylab("Density") +
      gg.opts + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                      strip.text.x = element_text(size = 10, face = "bold", colour = "black"),
                      legend.title = element_blank(),
                      plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, "cm"),
                      strip.background = element_rect(fill = "grey90")) +
      guides(colour = guide_legend(override.aes = list(size = 0.8)))
    
    if(!is.null(plot.model)) gPlot2 <- gPlot2 + 
      ggtitle("Multi-species dose-response model",
              subtitle = paste0("Model: ", which.model)) +
      theme(plot.title = element_text(size = 12, vjust = 5, face = "bold"),
            plot.subtitle = element_text(size = 11, vjust = 5))
    
    print(gPlot2)
  }
  
  
}