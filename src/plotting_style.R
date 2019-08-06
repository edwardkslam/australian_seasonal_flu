suppressPackageStartupMessages(library(tibble)) # for tibble()

subtype_colors <-
    list("B/Yam"="#CC79A7",
         "B/Vic"="#009E73",
         "H1sea"="#56B4E9",
         "H1pdm09"="#999999",
         "H3"="#E69F00")

posterior_pointcolor <- "#56c9ff" # for non-subtype-specific posteriors

default_posterior_limits <- c(-1.5, 1.5)


parameter_basenames <- c(
    'effect_antigenic_change',
    'effect_abs_humidity',
    'effect_is_first_of_season',
    'effect_cumulative_prior_inc',
    'effect_prior_season_activity',
    'effect_start_date')

parameter_display_names <- c(
    'antigenic change',
    'absolute humidity',
    'first epi of season',
    'prior variant cases',
    'prior season cases',
    'start date')


multilevel_incidence_parameter_names <- tibble(

    parameter_name = c(
        parameter_basenames,
        paste0("mean_", parameter_basenames),
        paste0("sd_", parameter_basenames)
    ),

    display_name = rep(parameter_display_names, 3)
)


cat("Loaded plotting style successfully")
