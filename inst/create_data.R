# Generate example dataset
library(magrittr)
library(tidyverse)
library(janitor)
library(ds4psy)
library(lubridate)

t_end <- lubridate::ceiling_date(lubridate::now(), "year") 
t_start <- t_end - lubridate::years(5)
set.seed(246)

template_brs <- readr::read_csv(file = "/Users/philippebouchet/Google Drive/Documents/git/rjmcmc/data/brsdata_3s_socal_autec.csv") %>% janitor::clean_names()

# All tibble columns shown, no colouring negative numbers, etc.
options(tibble.width = Inf) 
options(pillar.neg = FALSE) 
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

# dat <- simulate_data(n.species = 10,
#                      n.whales = c(1, 2, 4, 32, 15, 7, 22, 12, 9, 46),
#                      max.trials = 3,
#                      covariates = list(exposed = c(0, 5), sonar = c(0, 8), behaviour = c(0, 2, -12), range = 0.5),
#                      mu = runif(n = 10, min = 60, max = 215),
#                      phi = 20,
#                      sigma = 20,
#                      Lc = c(60, 65),
#                      Rc = c(180, 190),
#                      seed = 58690)

dat <- simulate_data(n.species = 7,
                     n.whales = c(1, 2, 4, 32, 15, 7, 22),
                     max.trials = 3,
                     covariates = list(exposed = c(0, 5), sonar = c(0, 8), behaviour = c(0, 2, -12), range = 0.5),
                     mu = runif(n = 7, min = 60, max = 215),
                     phi = 20,
                     sigma = 20,
                     Lc = c(60, 65),
                     Rc = c(180, 190),
                     seed = 58690)

# summary(dat)

n.obs <- length(dat$obs$y_ij)
dt_org <- sort(ds4psy::sample_time(from = t_start, to = t_end, size = n.obs))

ddf <- tibble::tibble(project = rep("Example BRS", n.obs))
ddf$species <- dat$species$id[dat$whales$id]
ddf <- ddf %>% 
  dplyr::mutate(species = dplyr::case_when(species == 1 ~ "Gg",
                                           species == 2 ~ "Ha",
                                           species == 3 ~ "Zc",
                                           species == 4 ~ "Bw",
                                           species == 5 ~ "Sw",
                                           species == 6 ~ "Oo",
                                           species == 7 ~ "Gm",
                                           TRUE ~ "other")) %>% 
  dplyr::mutate(indiv = dat$whales$id) %>% 
  dplyr::group_by(indiv) %>% 
  dplyr::mutate(tag_id = paste0(tolower(species), rpois(n =1, lambda = 10), "_", 
                                rpois(n = 1, lambda = runif(n = 1, min = 1, max = 1000)))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-indiv) %>% 
  dplyr::mutate(start_time = lubridate::ymd_hms(dt_org),
                resp_type = sample(x = template_brs$resp_type, replace = TRUE, size = n.obs),
                resp_spl = dat$obs$y_ij,
                resp_se_lcum = resp_spl + rgamma(n = 1, shape = 4, rate = 1),
                max_spl = dat$obs$Rc,
                max_se_lcum = max_spl + rgamma(n = 1, shape = 4, rate = 1),
                censored = dat$obs$censored,
                exp_order = dat$covariates$values$exposed,
                exp_duration = round(runif(n.obs, min = min(template_brs$exp_duration, na.rm = TRUE), max = max(template_brs$exp_duration, na.rm = TRUE)), 1),
                exp_signal = dat$covariates$values$sonar,
                pre_feeding = dat$covariates$values$behaviour)


mfas <- c("MFAS", "MFA", "REAL MFA", "MFAS_DS", "MFA HELO", "REAL MFA HELO", "SOCAL_d", "REAL MFAS", "MF ALARM")
lfas <- c("LFAS", "LFAS_DS", "LFAS_LO", "LFA")

ddf <- ddf %>%
  dplyr::mutate(min_dose = dat$obs$Lc) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(exp_signal = ifelse(exp_signal == "MFAS", sample(mfas, size = 1, replace = TRUE), 
                              sample(lfas, size = 1, replace = TRUE)),
                pre_feeding = ifelse(pre_feeding == "Feed", TRUE, FALSE),
                resp_time = start_time + runif(n = 1, min = 1, max = 60 * exp_duration),
                max_spl = ifelse(censored == 0, runif(n = 1, min = resp_spl, max = 215), 
                                 ifelse(censored == -1, min_dose, max_spl)),
                max_se_lcum = ifelse(censored == 0, runif(n = 1, min = resp_se_lcum, max = 215), max_se_lcum)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-min_dose)

range.row <- template_brs[, c("resp_range", "min_range", "inferred_resp_range", "inferred_min_range")] %>% 
  dplyr::filter_all(dplyr::any_vars(!is.na(.)))
range.ind <- sample(1:nrow(range.row), size = nrow(ddf), replace = TRUE)
range.val <- range.row[range.ind, ]

ddf <- ddf %>% dplyr::bind_cols(., range.val)

ddf <- ddf %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(resp_score = ifelse(censored == 0, sample(4:7, 1, TRUE), sample(0:3, 1, TRUE))) %>% 
  dplyr::mutate(resp_spl = ifelse(resp_score > 0 & censored == 1, runif(n = 1, 60, max_spl), resp_spl)) %>% 
  dplyr::mutate(resp_se_lcum = ifelse(resp_score > 0 & censored == 1, 
                                      runif(n = 1, resp_spl, max_se_lcum), resp_se_lcum)) %>% 
  dplyr::mutate(resp_spl = ifelse(censored == -1, max_spl, resp_spl),
                resp_se_lcum = ifelse(censored == -1, max_spl, resp_se_lcum)) %>% 
  dplyr::ungroup()
  
for(j in 1:nrow(ddf)){
  amount <- rgamma(n = 1, shape = 0.5, scale = 2)
  ddf$resp_range[j] <- ddf$resp_range[j] + amount
  ddf$min_range[j] <- ddf$min_range[j] + amount
  ddf$inferred_resp_range[j] <- ddf$inferred_resp_range[j] + amount
  ddf$inferred_min_range[j] <- ddf$inferred_min_range[j] + amount
  
  if(ddf$censored[j] == 1 & is.na(ddf$min_range[j])) ddf$min_range[j] <- ddf$resp_range[j]
  
}



example_brs <- ddf %>% 
  dplyr::select(project, species, tag_id, start_time, resp_time, resp_type, resp_score, resp_spl, resp_se_lcum, max_spl, max_se_lcum, censored, exp_order, exp_duration, exp_signal, pre_feeding, resp_range, min_range, inferred_resp_range, inferred_min_range)

usethis::use_data(example_brs, overwrite = TRUE)
