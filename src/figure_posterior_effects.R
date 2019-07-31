#!/usr/bin/env Rscript

#####################################
## name: figure_posterior_effects.R
## version 0.0.1
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for effects
## from regression
##
####################################

script_packages <- c(
    'ggplot2',    # general plotting
    'tibble',     # for tibble()
    'dplyr',      # for do()
    'magrittr',   # for pipe operator
    'tidybayes',  # for spread_draws()
    'readr',      # for read_csv()
    'cowplot')    # for publication-ready ggplot

for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}

args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
plotting_params_path <- args[3]
save_path <- args[4]

dat <- read_csv(data_path,
                col_types = cols())

source(plotting_params_path)

fit <- readRDS(mcmc_fit_path)

parameter_names = tibble(

    parameter_name = c(
        'effect_antigenic_change',
        'effect_abs_humidity',
        'effect_is_first_of_season',
        'effect_cumulative_prior_inc',
        'effect_prior_season_activity',
        'effect_start_date',
        'effect_temperature'),

    display_name = c(
        'antigenic change',
        'absolute humidity',
        'first epi. of season',
        'prior variant cases',
        'prior season cases',
        'start date',
        'temperature')
)
    
tidychains <- fit %>% gather_draws(effect_abs_humidity[subtype_id],
                                   effect_cumulative_prior_inc[subtype_id],
                                   effect_antigenic_change[subtype_id],
                                   effect_prior_season_activity[subtype_id],
                                   effect_is_first_of_season[subtype_id],
                                   effect_temperature[subtype_id],
                                   effect_start_date[subtype_id])

quants <- tidychains %>%
    do(tibble(post_quant = quantile(.$.value, ppoints(100))))
quants <- quants %>%
    ungroup() %>%
    rename(parameter_name = .variable)

quants$parameter_name <- factor(quants$parameter_name)

pointcolor <- "#56c9ff"

quants <- quants %>% left_join(parameter_names,
                               by = 'parameter_name')

effect_fig <- quants %>%
    ggplot(aes(x = post_quant)) +
    geom_dotplot(
        fill=pointcolor,
        alpha=1,
        binwidth=0.1) +
    geom_vline(xintercept=0) + 
    scale_x_continuous(limits=c(-2, 2)) + 
    facet_wrap(vars(display_name), ncol=2) +
    xlab("regression coeffecient\n(common effect size scale)") +
    ylab("posterior frequency") +
    theme_classic() +
    theme_cowplot() +
    panel_border()


save_plot(save_path,
          effect_fig,
          base_height=9,
          base_aspect_ratio=0.75)


