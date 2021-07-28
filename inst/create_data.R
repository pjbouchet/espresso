# Generate example dataset
library(tidyverse)
library(magrittr)

set.seed(1587)

# All tibble columns shown, no colouring negative numbers, etc.
options(tibble.width = Inf) 
options(pillar.neg = FALSE) 
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

# Import species list and raw data

species.list <- readr::read_csv("../rjmcmc/data/species_list.csv", col_types = cols()) %>% 
  janitor::clean_names(.) %>% 
  dplyr::mutate(code = tolower(code))

rawdat <- readr::read_csv("../rjmcmc/data/brsdata_3s_socal_autec.csv", col_types = cols()) %>% 
  janitor::clean_names(.) %>% 
  dplyr::mutate(species = tolower(species))

# Project 

dat <- rawdat %>% 
  dplyr::mutate(project = "Example_project") %>% 
  dplyr::left_join(x = ., y = species.list, by = c("species" = "code"))

# Filter species and generate random tag IDs

dat <- dat %>% 
  dplyr::filter(species %in% c("ba", "oo", "bw", "gg", "zc", "ha", "md", "sw")) %>% 
  dplyr::group_by(tag_id) %>% 
  dplyr::mutate(tag_id = paste0(tolower(species), rpois(n =1, lambda = 10), "_", rpois(n =1, lambda = runif(n = 1, min = 1, max = 1000)))) %>% 
  dplyr::ungroup()

# Summarise mean SPL by species

spl.by.species <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(spl_mean = mean(resp_spl, na.rm = TRUE), spl_sd = sd(resp_spl, na.rm = TRUE)) %>% 
  dplyr::ungroup()

spl.by.species[spl.by.species$species == "ba", ]$spl_mean <- runif(1, 60, 215)
spl.by.species[spl.by.species$species == "bw", ]$spl_mean <- runif(1, 60, 215)
spl.by.species[spl.by.species$species %in% c("ha", "zc"), ]$spl_mean <- rep(runif(1, 60, 215), 2)
spl.by.species[spl.by.species$species == "md", ]$spl_mean <- runif(1, 60, 215)
spl.by.species[spl.by.species$species == "oo", ]$spl_mean <- runif(1, 60, 215)
spl.by.species[spl.by.species$species == "sw", ]$spl_mean <- runif(1, 60, 215)

dat <- dplyr::left_join(x = dat, y = spl.by.species, by = "species")

dat <- dat %>% dplyr::rowwise() %>% 
  dplyr::mutate(resp_spl = ifelse(is.na(resp_spl), NA, 
                              rtnorm(n = 1, location = spl_mean, scale = ifelse(is.na(spl_sd), 2, spl_sd), L = 60, U = 190))) %>%
  dplyr::mutate(resp_se_lcum = ifelse(is.na(resp_se_lcum), NA, 
                                      resp_spl + rgamma(n = 1, shape = 1, rate = 1))) %>%
  dplyr::mutate(max_spl = ifelse(censored == 0, runif(n = 1, min = resp_spl, max = 200), max_spl)) %>% 
  dplyr::ungroup()

# Permute values
dat$exp_signal <- sample(dat$exp_signal, size = nrow(dat), replace = FALSE)
dat$exp_duration <- round(runif(n = nrow(dat), min = 8, max = 70), 0)
dat$resp_type <- sample(dat$resp_type, size = nrow(dat), replace = FALSE)

for(j in 1:nrow(dat)){
  if(!is.na(dat$resp_spl[j])){
    dat$resp_se_lcum[j] <- dat$resp_spl[j] + rgamma(1, 0.5, 4)
  }
  
  if(!is.na(dat$max_spl[j])){
    dat$max_se_lcum[j] <- dat$max_spl[j] + rgamma(1, 0.5, 4)
  }
}


# dat$resp_score <- sample(1:7, size = nrow(dat), replace = TRUE)

# dat$resp_score[dat$resp_type == "Any"] <- 0
# dat$resp_score[is.na(dat$resp_score) | dat$resp_score == "not scored"] <- 7
# dat$resp_score <- as.numeric(dat$resp_score)

# Right-censoring

# is.censored <- numeric(nrow(dat))
# for(i in seq_len(nrow(dat))){
#   df <- dat[i, ]
#   if(df$resp_score >=4){
#     if(is.na(df$resp_spl)){
#       is.censored[i] <- 1
#     } else {
#       if(df$resp_spl <= df$max_spl){
#         is.censored[i] <- 0
#       } else {
#         is.censored[i] <- 1
#       }
#     }
#   } else {
#     is.censored[i] <- 1
#   }
# }
# 
# dat$censored <- is.censored
# dat$resp_spl[is.censored == 1] <- NA

dat <- dat %>% 
  dplyr::mutate(pre_feeding = sample(x = dat$pre_feeding, size = nrow(dat), replace = FALSE),
                start_time = sample(x = dat$start_time, size = nrow(dat), replace = FALSE)) 

# dat <- dat %>% 
#   dplyr::mutate(pre_feeding = sample(x = dat$pre_feeding, size = nrow(dat), replace = FALSE),
#                 start_time = sample(x = dat$start_time, size = nrow(dat), replace = FALSE),
#                 max_spl = rnorm(n = nrow(dat), mean = mean(dat$max_spl), sd = sd(dat$max_spl))) %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(max_spl = ifelse(max_spl < resp_spl & !is.na(resp_spl), resp_spl, max_spl)) %>% 
#   dplyr::ungroup()

# dat <- dat %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(resp_se_lcum = resp_spl + abs(rnorm(n = 1, mean = 0, sd = 10)),
#                 max_se_lcum = max_spl + abs(rnorm(n = 1, mean = 0, sd = 10))) %>% 
#   dplyr::ungroup()


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

example_brs <- dat %>% 
  dplyr::select(-scientific_name, -common_name, -new_code, -spl_mean, -spl_sd) %>% 
  dplyr::select(project, species, tag_id, start_time, resp_time, resp_type, resp_score, resp_spl, resp_se_lcum, max_spl, max_se_lcum, censored, exp_order, exp_duration, exp_signal, pre_feeding, resp_range, min_range, inferred_resp_range, inferred_min_range)
# names(example_brs) <- names(rawdat)
usethis::use_data(example_brs, overwrite = TRUE)
