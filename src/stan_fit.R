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


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <- "multilevel_incidence_model.stan"
datapath <- "../epi_table.csv"
humid_datapath <- "../mean_fortnightly_climate_30years.csv"

## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 1500
nchains = n_cores
adapt_d= 0.8
max_tree = 15
fixed_seed = 232032

## load data
dat <- read_csv(datapath,
                col_types = cols())
## make data into list
data_list <- list(
    n_cities = max(dat$city_id),
    n_epidemics = length(dat$epi_z_score),
    incidences = dat$epi_z_score,
    city = dat$city_id,
    antigenic_change = dat$new_ag_marker,
    mean_abs_humidity = dat$mean_AH_over_epi,
    prior_activity = dat$prior_everything_scaled,
    cumulative_prior_incidence = dat$frac_max)

hyperparam_list <- list(
    mean_city_reporting_rates_per_hundred=0,
    sd_city_reporting_rates_per_hundred=0.25,
    alpha_average_epi_attack_rate=2,
    beta_average_epi_attack_rate=15,
    sd_sd_incidences=2)

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

warnings()
