#!/usr/bin/env Rscript

################################
## fit the specified stan model
## to the specified clean dataset,
## and save the result as an .Rds
## file
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(shinystan))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <- args[1]
datapath <- args[2]
mcmc_output_path <- args[3]

## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 2000
nchains = n_cores
adapt_d= 0.8
max_tree = 15
fixed_seed = 232032

## load data
dat <- read_csv(datapath,
                col_types = cols())

if(grepl('normed', model_src_path)){
    dat <- dat[!is.na(dat$mean_centered_log_epi_size),]
} else {
    dat <- dat[!is.na(dat$new_ag_marker),]
}


## nuisance value for ignored standardized vals
## should screw things up visibly if model wrong
dat$standardized_log_prior_cumulative <-
    dat$standardized_log_prior_cumulative %>% replace_na(9999)
dat$standardized_prior_season_activity <-
    dat$standardized_prior_season_activity %>% replace_na(9999)


## make data into list
data_list <- list(
    n_cities = max(dat$city_id),
    n_epidemics = length(dat$mean_centered_log_epi_size),
    n_possible_epidemics = length(dat$epidemic_flag),
    n_subtypes = max(dat$subtype_id),
    subtype = dat$subtype_id,
    incidences = log(dat$incidence_per_mil),
    normed_metric = dat$mean_centered_log_epi_size,
    epi_occurred = dat$epidemic_flag,
    city = dat$city_id,
    antigenic_change = dat$new_ag_marker,
    cumulative_prior_incidence_std = dat$standardized_log_prior_cumulative,
    abs_humidity = dat$mean_epi_ah,
    temp = dat$mean_epi_temp,
    rainfall = dat$mean_epi_rainfall,
    is_first_of_season = dat$is_first_of_season,
    prior_season_activity_std = dat$standardized_prior_season_activity,
    start_date = dat$start)


hyperparam_list <- list(
    mean_city_reporting_rates_per_mil = 0,
    sd_city_reporting_rates_per_mil = 5000,
    alpha_average_epi_attack_rate = 2,
    beta_average_epi_attack_rate = 10,
    sd_sd_incidences = 1,
    sd_mean_effect_sizes = 1,
    sd_sd_effect_sizes = 0.25,
    sd_mean_intercept = 1,
    sd_sd_intercept = 0.5,
    nu = 3)

stan_data <- c(
    data_list,
    hyperparam_list)


###############################
## Compile, fit, and save model
###############################
fit <- stan(
    model_src_path,
    data=stan_data,
    iter=niter,
    seed=fixed_seed,
    chains=nchains, 
    control = list(max_treedepth = max_tree,
                   adapt_delta = adapt_d))

saveRDS(fit, mcmc_output_path)

warnings()
