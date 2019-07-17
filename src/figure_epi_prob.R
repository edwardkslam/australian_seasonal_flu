#!/usr/bin/env Rscript

#####################################
## name: figure_epi_prob.R
## date: 2019-07-15
## version 0.0.1
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for
## epidemic probability and related
## quantities
##
####################################

script_packages <- c(
    'readr',     # for read_csv()
    'ggplot2',   # for plotting
    'cowplot',   # for publication-ready ggplot
    'magrittr',  # for pipe operator %>%
    'tidybayes', # for spread_draws(),
    'dplyr',     # for sql-style manipulations
    'tidyr')     # for cartesian product via crossing(),

## silently load packages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


################################
## process script call inputs
## and read in data
################################


args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
plotting_style_path <- args[3]
save_path <- args[4]

set.seed(234902385) # reproducible!

source(plotting_style_path) # get plotting params

dat <- read_csv(data_path,
                col_types = cols())
fit <- readRDS(mcmc_fit_path)


## add human-readable parameter names
parameter_names = tibble(

    parameter_name = c(
        'effect_antigenic_change',
        'effect_cumulative_prior_inc'),

    display_name = c(
        'antigenic change',
        'prior variant cases')
)


##############################
## plot posterior dotplots
##############################

tidychains <- fit %>% gather_draws(effect_cumulative_prior_inc[subtype_id],
                                   effect_antigenic_change[subtype_id])

## SQL-ishly reintroduce subtype human-readable names
tidychains <-
    tidychains %>%
    inner_join(
        distinct(dat, subtype_id, subtype),
        by='subtype_id') %>%
    left_join(parameter_names, by=c('.variable'='parameter_name')) %>% 
    group_by(subtype, display_name)


quants <- tidychains %>%
    do(tibble(post_quant = quantile(.$.value, ppoints(100))))
quants <- quants %>% ungroup()

effect_fig <- quants %>%
    ggplot(aes(x = post_quant,
               fill = subtype)) +
    geom_dotplot(
        alpha=1,
        binwidth=.2) +
    geom_vline(xintercept=0) + 
    scale_fill_manual(values = subtype_colors) + 
    facet_grid(rows=vars(subtype), cols=vars(display_name)) +
    xlab("regression coeffecient\n(common effect size scale)") +
    ylab("posterior frequency") +
    theme_classic() +
    theme_cowplot() +
    theme(legend.position = "none") + 
    panel_border()

#############################################
## plot distribution of regression lines
#############################################

save_plot(save_path,
          effect_fig,
          base_height=9,
          base_aspect_ratio=1.61)
