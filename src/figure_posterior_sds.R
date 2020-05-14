#!/usr/bin/env Rscript

#####################################
## name: figure_posterior_sds.R
## version 0.0.1
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for sds
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

###################################
## fixed and variable plot style
###################################

source(plotting_params_path)
## includes multilevel_incidence_parameter_names
## and default_posterior_limits

bin_width <- 0.028
dot_size <- 1
fig_n_col <- 3
sd_posterior_limits <- c(0, 1)

fit <- readRDS(mcmc_fit_path)

## determine whether temperature
## or abs humidity is the climate
## variable
if(grepl("temp", mcmc_fit_path)){
    tidychains <- fit %>% gather_draws(
                              sd_intercept,
                              sd_effect_temp,
                              sd_effect_rainfall,
                              sd_effect_temp,
                              sd_effect_cumulative_prior_inc,
                              sd_effect_prior_season_activity,
                              sd_effect_is_first_of_season,
                              sd_effect_start_date)
} else{
    tidychains <- fit %>% gather_draws(
                              sd_intercept,
                              sd_effect_abs_humidity,
                              sd_effect_rainfall,
                              sd_effect_antigenic_change,
                              sd_effect_cumulative_prior_inc,
                              sd_effect_prior_season_activity,
                              sd_effect_is_first_of_season,
                              sd_effect_start_date)
}

quants <- tidychains %>%
    do(tibble(post_quant = quantile(.$.value, ppoints(100))))
quants <- quants %>%
    ungroup() %>%
    rename(parameter_name = .variable)

quants <- quants %>%
    left_join(multilevel_incidence_parameter_names,
              by = 'parameter_name')

effect_fig <- quants %>%
    ggplot(aes(x = post_quant)) +
    geom_dotplot(
        fill=posterior_pointcolor,
        alpha=1,
        binwidth=bin_width,
        dotsize=dot_size,
        method='histodot') +
    scale_x_continuous(limits = sd_posterior_limits) + 
    facet_wrap(vars(display_name), ncol=fig_n_col) +
    xlab("standard deviation") +
    ylab("posterior frequency") +
    theme_cowplot(font_size=25) +
    panel_border()


save_plot(save_path,
          effect_fig,
          base_height=15,
          base_aspect_ratio=1.2)


