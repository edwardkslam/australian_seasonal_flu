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

#######################################
## plot styling 
#######################################
source(plotting_params_path)
bin_width <- 0.1
dot_size <- 1
multilevel_incidence_parameter_names = tibble(

    parameter_name = c(
        'effect_antigenic_change',
        'effect_abs_humidity',
        'effect_is_first_of_season',
        'effect_cumulative_prior_inc',
        'effect_prior_season_activity',
        'effect_start_date'),

    display_name = c(
        'antigenic change',
        'absolute humidity',
        'first epi of season',
        'prior variant cases',
        'prior season cases',
        'start date')
)


fit <- readRDS(mcmc_fit_path)

tidychains <- fit %>% gather_draws(effect_abs_humidity[subtype_id],
                                   effect_cumulative_prior_inc[subtype_id],
                                   effect_antigenic_change[subtype_id],
                                   effect_prior_season_activity[subtype_id],
                                   effect_is_first_of_season[subtype_id],
                                   effect_start_date[subtype_id])

## SQL-ishly reintroduce subtype human-readable names
tidychains <-
    tidychains %>%
    inner_join(
        distinct(dat, subtype_id, subtype),
        by='subtype_id') %>%
    group_by(subtype, .variable)

quants <- tidychains %>%
    do(tibble(post_quant = quantile(.$.value, ppoints(100))))
quants <- quants %>%
    ungroup() %>%
    rename(parameter_name = .variable)

quants <- quants %>%
    left_join(multilevel_incidence_parameter_names,
              by = 'parameter_name')


effect_fig <- quants %>%
    ggplot(aes(x = post_quant,
               fill = subtype)) +
    geom_dotplot(
        alpha = 1,
        binwidth = bin_width,
        dotsize=dot_size,
        method='histodot') +
    geom_vline(xintercept=0) + 
    scale_x_continuous(limits = default_posterior_limits) +
    scale_fill_manual(values = subtype_colors) + 
    facet_grid(rows=vars(subtype), cols=vars(display_name)) +
    xlab("regression coeffecient\n(common effect size scale)") +
    ylab("posterior frequency") +
    theme_cowplot(font_size=20) +
    theme(legend.position = "none") + 
    panel_border()


save_plot(save_path,
          effect_fig,
          base_height=9,
          base_aspect_ratio=1.6)
