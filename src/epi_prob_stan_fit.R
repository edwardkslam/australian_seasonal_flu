#!/usr/bin/env Rscript

################################
## fit the specified epi probability stan model
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
model_src_path <- "multilevel_epi_probability_model.stan"
datapath <- "../dat/cleaned/clean_stan_data.csv"
mcmc_output_path <- '../out/mcmc_chains/prob_stan_fit_output.Rds'

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

dat <- dat[!is.na(dat$new_ag_marker),]


## make data into list
data_list <- list(
    n_cities = max(dat$city_id),
    n_possible_epidemics = length(dat$epidemic_flag),
    epi_occurred = dat$epidemic_flag,
    n_subtypes = max(dat$subtype_id),
    subtype = dat$subtype_id,
    city = dat$city_id,
    antigenic_change = dat$new_ag_marker,
    abs_humidity = dat$mean_epi_ah,
    temperature = dat$mean_epi_temp,
    other_subtype_activity = dat$prior_everything_scaled,
    cumulative_prior_incidence = dat$standardised_prior_cumulative,
    start_date = dat$start)

hyperparam_list <- list(
    sd_mean_effect_sizes = 1,
    sd_sd_effect_sizes = 0.1,
    sd_mean_intercept = 1,
    sd_sd_intercept = 1)

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
