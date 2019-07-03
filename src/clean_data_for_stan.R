#!/usr/bin/env Rscript

#####################################
## clean data for clear, reproducible
## input into Bayesian regression
## analysis
######################################


#######################
# load needed packages
#######################
script_packages <- c(
    'magrittr',   # for (ceci n'est pas une) pipe operator %>%
    'readr',      # for read_csv()
    'dplyr',      # for group_by()
    'tidyr',      # for full_seq()
    'tibble')     # for add_column() 
    
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


## function for running data integrity checks before
## saving
check_data_integrity <- function(dat){
    ag_markers_right <- xor(
        dat$standardised_prior_cumulative > 0.001 &
        (is.na(dat$new_ag_marker) | (!dat$new_ag_marker)),
        dat$standardised_prior_cumulative < 0.001 &

cumulative_incidence_by_ag$standardised_prior_cumulative <-
    cumulative_incidence_by_ag$cumulative_prior_incidence_same_ag /
    cumulative_incidence_by_ag$mean_epi_size

cumulative_incidence_by_ag$standardised_current_season <-
    cumulative_incidence_by_ag$incidence_current_season /
    cumulative_incidence_by_ag$mean_epi_size

dat <- cumulative_incidence_by_ag %>%
    subset(!is.na(incidence_current_season))

dat$city_id <- as.numeric(
    factor(dat$city, 
           levels=unique(dat$city)))

epi_variation <- epi_dat %>% group_by(subtype, city) %>%
    summarise(
        mean_log_inc = mean(
            log(incidence_per_mil),
            na.rm=TRUE),

        sd_log_inc = sd(
            log(incidence_per_mil),
            na.rm=TRUE)
    )

epi_variation$cv_log_inc <- epi_variation$sd_log_inc / epi_variation$mean_log_inc

## calculate epi_z_score
dat <- dat %>% left_join(epi_variation,
                         by = c('subtype', 'city'))

dat$epi_z_score <- (log(dat$incidence_per_mil) - dat$mean_log_inc) / dat$sd_log_inc


no_epi <- dat$epi_alarm=='N' |
    is.na(dat$epi_alarm)
dat$epidemic_flag <- ifelse(
    no_epi,
    0,
    1)


## output only needed columns and rows
columns_wanted = c("city_id",
                   "city",
                   "subtype",
                   "year",
                   "start",
                   "end",
                   "incidence_per_mil",
                   "epidemic_flag",
                   "prior_everything_scaled",
                   "cumulative_prior_incidence_same_ag",
                   "standardised_prior_cumulative",
                   "new_ag_marker",
                   "mean_log_inc",
                   "sd_log_inc",
                   "cv_log_inc",
                   "epi_z_score")
## data integrity checks
cat("\nChecking data integrity...\n\n")
check_data_integrity(dat)

dat <- dat[, columns_wanted]

dat <- dat %>% arrange(city_id, subtype, year, new_ag_marker)

write.csv(dat,
          file = output_path,
          row.names=FALSE)

warnings()
