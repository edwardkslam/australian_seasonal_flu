#!/usr/bin/env Rscript

#####################################
## name: extract_reported_quantities.R
## date: 2019-03-19
## version 0.0.1
## author: Dylan Morris <dhmorris@princeton.edu>
##
## read in mcmc chains
## and output quantities that
## we report in the paper
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot)) # for rhat, neff_ratio

## read command line args
args <- commandArgs(trailingOnly=TRUE)
n_models <- length(args) - 2
mcmc_fit_paths <- args[1:n_models]
data_path <- args[n_models + 1]
outpath <- args[n_models + 2]

## read data
dat <- read.csv(data_path)


list_of_quantities <- list()

## do regression model chain quantity extraction
cat("\n Reading in regression model...\n")
fit <- readRDS(mcmc_fit_paths[1])
chains <- extract(fit)
cat(paste0("Model read successfully, ",
           "extracting reported ",
           "quantities...\n"))

for(quantity in c('mean_effect_antigenic_change',
                  'mean_effect_abs_humidity',
                  'mean_effect_cumulative_prior_inc',
                  'mean_effect_prior_season_activity',
                  'mean_effect_start_date',
                  'mean_effect_is_first_of_season',
                  'sd_incidences',
                  'sd_effect_antigenic_change',
                  'sd_effect_abs_humidity',
                  'sd_effect_cumulative_prior_inc',
                  'sd_effect_prior_season_activity',
                  'sd_effect_start_date',
                  'sd_effect_is_first_of_season')){
    print(quantile(chains[[quantity]]))
    list_of_quantities[paste0(quantity, '_q025')] <-
        quantile(chains[[quantity]], 0.025)
    list_of_quantities[paste0(quantity, '_q975')] <-
        quantile(chains[[quantity]], 0.975)
    list_of_quantities[paste0(quantity, '_q50')] <-
        quantile(chains[[quantity]], 0.5)
}

## key-value to dataframe
df <- stack(list_of_quantities)
names(df) <- c("value", "quantity")

## reorder for better display
df <- df[c("quantity", "value")]
## save to outpath

## output to file
cat(sprintf("saving to %s\n",
            outpath))
write.csv(df, outpath)
