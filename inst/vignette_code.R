## OBTAINED WITH KNITR::PURL
## ----eval=FALSE, include=TRUE---------------------------------------------------------------
## # install.packages("remotes")
## remotes::install_github("pjbouchet/espresso", dependencies = TRUE)


## ----echo=T, results='hide', message = F----------------------------------------------------
library(espresso)
library(tidyverse)
library(magrittr)

#'--------------------------------------------------------------------
# Set tibble options
#'--------------------------------------------------------------------
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

Sys.setenv(TZ = "GMT")

#'--------------------------------------------------------------------
# Set knitr options
#'--------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
knitr::kable(head(example_brs), format = "pandoc")
knitr::kable(head(species_brs), format = "pandoc")


## ----eval=FALSE-----------------------------------------------------------------------------
## mydat <- read_data(file = "path/to/my/data.csv")


## ----eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------------------------------
## mydat <- read_data(file = NULL)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
mydat <- read_data(file = NULL, exclude.species = "Risso's dolphin", min.N = 3)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
summary(mydat)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
mydat <- read_data(file = NULL, 
                   exclude.species = c("Risso's dolphin", "Tursiops truncatus"),
                   min.N = 2,
                   covariates = c("behaviour", "range"))


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
summary(mydat)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
mydat.grp <- create_groups(dat.obj = mydat, abbrev = TRUE,
            species.groups = list(Beaked_whales = c("Ha", "Cuvier's beaked whale")))


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
class(mydat.grp)
summary(mydat.grp)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
mydat.ungrp <- undo_groups(mydat.grp)
class(mydat.ungrp)
summary(mydat.ungrp)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
mydat.config <- configure_rjMCMC(dat = mydat.grp,
                                 model.select = TRUE,
                                 covariate.select = TRUE,
                                 n.rep = 100)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
class(mydat.config)
summary(mydat.config, print.config = TRUE)


## ----modelfit, eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE--------------------------
## rj <- run_rjMCMC(dat = mydat.config, n.chains = 3, n.burn = 1000, n.iter = 10000)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
multicool::Bell(n = 5) # Number of candidate models

# List of candidate models
partitions::listParts(x = 5) %>%
  purrr::map_depth(.x = ., .depth = 2, .f = ~mydat.config$species$names[.x]) %>%
  purrr::map(.x = ., 
             .f = ~ lapply(X = .x, FUN = function(x) paste(x, collapse = ","))) %>%
  purrr::map(.x = ., .f = ~ paste0("(", .x, ")")) %>% 
  purrr::map(.x = ., .f = ~paste0(.x, collapse = "+")) %>% 
  tibble::enframe() %>% dplyr::select(-name) %>% 
  dplyr::rename(model = value) %>% data.frame()


## ----message=FALSE, warning=FALSE, include=FALSE--------------------------------------------
rj.posterior <- espresso:::vignette.posterior


## ----eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------------------------------
## rj.posterior <- trace_rjMCMC(rj, thin = 10)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
class(rj.posterior)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
summary(rj.posterior, combine.chains = FALSE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(rj.posterior, individual = T)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
# Plot for the range covariate only
plot(rj.posterior, param.name = "range", dautocorr = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
doseR <- compile_rjMCMC(rj.posterior)
class(doseR)
plot(doseR)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
doseR.beh <- compile_rjMCMC(rj.posterior, covariate = "behaviour", species = "Oo")
plot(doseR.beh)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
doseR.rge <- compile_rjMCMC(rj.posterior, 
                            covariate = "range", 
                            covariate.values = c(5, 10),
                            species = "Oo")
plot(doseR.rge, scientific.name = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(doseR, colour = "gray20", colour.median = "black")


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(doseR, colour.by = "species", order.by = "response")


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(doseR, scientific.name = TRUE, show.p0_5 = FALSE, outline.outer = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(doseR, colour.by = "species", overlay = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
doseR <- compile_rjMCMC(rj.posterior, by.model = TRUE)


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
plot(doseR, model.rank = 1, order.by = "response")


## ----echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------
doseR


## ----echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE------------------------------------
## create_report(outdir = "path/to/directory", filename = "file.name")

