# Generate example dataset
library(tidyverse)
library(magrittr)

set.seed(1587)

# All tibble columns shown, no colouring negative numbers, etc.
options(tibble.width = Inf) 
options(pillar.neg = FALSE) 
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

species.list <- readr::read_csv("../rjmcmc/data/species_list.csv", col_types = cols()) %>% 
  janitor::clean_names(.) %>% 
  dplyr::mutate(code = tolower(code))

rawdat <- readr::read_csv("../rjmcmc/data/brsdata_3s_socal_autec.csv", col_types = cols()) %>% 
  janitor::clean_names(.) %>% 
  dplyr::mutate(species = tolower(species))

dat <- rawdat %>% 
  dplyr::mutate(project = "Example_project") %>% 
  dplyr::left_join(x = ., y = species.list, by = c("species" = "code"))

dat <- dat %>% 
  dplyr::filter(species %in% c("ba", "oo", "bw", "gg", "zc", "ha", "md", "sw")) %>% 
  dplyr::group_by(tag_id) %>% 
  dplyr::mutate(tag_id = paste0(tolower(species), rpois(n =1, lambda = 10), "_", rpois(n =1, lambda = runif(n = 1, min = 1, max = 1000)))) %>% 
  dplyr::ungroup()

dat$exp_signal <- sample(dat$exp_signal, size = nrow(dat), replace = FALSE)
dat$exp_duration <- round(runif(n = nrow(dat), min = 8, max = 70), 0)
dat$resp_type <- sample(dat$resp_type, size = nrow(dat), replace = FALSE)
dat$resp_score <- sample(1:7, size = nrow(dat), replace = TRUE)
dat$resp_score[dat$resp_type == "Any"] <- 0

spl.by.species <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(spl_mean = mean(resp_spl, na.rm = TRUE), spl_sd = sd(resp_spl, na.rm = TRUE)) %>% 
  dplyr::ungroup()

spl.by.species[spl.by.species$species == "gg", ]$spl_mean <- 125
spl.by.species[spl.by.species$species == "gg", ]$spl_sd <- 10
spl.by.species[spl.by.species$species == "ba", ]$spl_sd <- 10
spl.by.species[spl.by.species$species == "md", ]$spl_sd <- 10

spl.by.species$species <- sample(x = spl.by.species$species, size = nrow(spl.by.species), replace = FALSE)

spl.na <- sum(is.na(dat$resp_spl))
dat$resp_spl <- purrr::map_dbl(.x = 1:nrow(dat), 
                               .f = ~{
                                 current.spl <- dat[.x, ]$resp_spl
                                 current.species <- dat[.x, ]$species
                                 rnorm(n = 1, 
                                         mean = spl.by.species[spl.by.species$species == current.species, ]$spl_mean,
                                         sd = spl.by.species[spl.by.species$species == current.species, ]$spl_sd)
                               })

is.censored <- numeric(nrow(dat))
for(i in seq_len(nrow(dat))){
  df <- dat[i, ]
  if(!df$resp_score >=4){
    is.censored[i] <- 1
  } else {
    if(is.na(df$resp_spl)){
      is.censored[i] <- 1
    } else {
      if(df$resp_spl < df$max_spl){
        is.censored[i] <- 0
      } else {
        is.censored[i] <- 1
      }
    }
  }
}

dat$censored <- is.censored
dat$resp_spl[is.censored == 1] <- NA

dat <- dat %>% 
  dplyr::mutate(pre_feeding = sample(x = dat$pre_feeding, size = nrow(dat), replace = FALSE),
                start_time = sample(x = dat$start_time, size = nrow(dat), replace = FALSE),
                max_spl = rnorm(n = nrow(dat), mean = mean(dat$max_spl), sd = sd(dat$max_spl))) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(max_spl = ifelse(max_spl < resp_spl & !is.na(resp_spl), resp_spl, max_spl)) %>% 
  dplyr::ungroup()

dat <- dat %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(resp_se_lcum = resp_spl + abs(rnorm(n = 1, mean = 0, sd = 10)),
                max_se_lcum = max_spl + abs(rnorm(n = 1, mean = 0, sd = 10))) %>% 
  dplyr::ungroup()


# Range

for(i in 1:nrow(dat)){
  if(!is.na(dat$resp_range[i])) dat$resp_range[i] <- rgamma(n = 1, shape = 1, rate = 0.1)
  if(!is.na(dat$min_range[i])) dat$min_range[i] <- dat$resp_range[i] - rgamma(n = 1, shape = 0.01, rate = 2)
  
  if(dat$censored[i] == 0){ # Response range
    if(is.na(dat$resp_range[i])){
      dat$inferred_resp_range[i] <- rgamma(n = 1, shape = 1, rate = 0.1)
    } else {
      dat$inferred_resp_range[i] <- NA
    }
  } else {
    if(is.na(dat$min_range[i])){
      dat$inferred_min_range[i] <- rgamma(n = 1, shape = 1, rate = 0.1)
    } else {
      dat$inferred_min_range[i] <- NA
    }
  }
  
}

example_brs <- dat %>% dplyr::select(-scientific_name, -common_name, -new_code)
names(example_brs) <- names(rawdat)
usethis::use_data(example_brs, overwrite = TRUE)