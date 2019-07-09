#!/usr/bin/env Rscript

#####################################
## name: figure_attack_rate.R
## version 0.0.1
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for difference between
## attack rates with and without antigenic
## change
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
save_path <- args[3]

dat <- read_csv(data_path,
                col_types = cols())

fit <- readRDS(mcmc_fit_path)

    
## tidychains <- fit %>% gather_draws(epi_no_change_attack[epi_id],
##                                    epi_with_change_attack[epi_id])

tidychains <- fit %>% gather_draws(epi_no_change[epi_id],
                                   epi_with_change[epi_id])



attack_fig <- tidychains %>%
    ggplot(aes(x = exp(.value),
               fill = .variable)) +
    geom_density(
        alpha=0.5) +
    scale_fill_manual(values=c('gray','red')) + 
    xlab("true cases") +
    ylab("posterior density") +
    scale_x_log10() +
    theme_classic() +
    theme_cowplot() +
    panel_border()


save_plot(save_path,
          attack_fig,
          base_height=9,
          base_aspect_ratio=1)


