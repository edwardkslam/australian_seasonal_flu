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
        dat$new_ag_marker)
    if(all(ag_markers_right))
    {
        cat("passed integrity check for prior activity\n\n")
    } else {
        cat("failed integrity check for prior activity\n\n")
        print(dat[which(!ag_markers_right),])
    }
}



##########################
## read command line args,
## load data
##########################
args <- commandArgs(trailingOnly=TRUE)
epi_data_path <- args[1] #  epi_data_path <- "../epi_table.csv"
output_path <- args[2]

## load data
epi_dat <- read_csv(epi_data_path,
                    col_types = cols())

epi_dat$city <- factor(epi_dat$city,
                       levels = unique(epi_dat$city))


## find variants whose cumulative incidence
## cannot be evaluated
no_cumulative_incidence <- epi_dat %>%
    subset(.,is.na(new_ag_marker) | 
             reference_strain=="A/Perth/16/2009-like")

## First find the first and last year in
## which an antigenic variant causes an epidemic in each city

ag_years <- epi_dat %>%
    subset(.,epi_alarm == "Y") %>% 
    subset(.,!(reference_strain %in% no_cumulative_incidence$reference_strain)) %>% 
    group_by(city, subtype, reference_strain) %>%
    summarise(emerge_year=min(year),
              last_recorded_epi=max(year))

## assume that an ag variant could potentially have caused an epidemic right up to the year 
## before the emergence of its replacement new variant
ag_years <- ag_years %>%
  group_by(city, subtype) %>%
  arrange(emerge_year,.by_group=TRUE) %>%
    mutate(last_possible_year=ifelse(
               subtype!="H1sea",
               lead(emerge_year - 1, default = 2015),
               lead(emerge_year - 1, default = 2008)))

## the last possible year in which an epidemic could have been caused is whichever occured last:
## the last recorded epidemic (to account for the case in which old and new variants circulate in same year)
## the year prior to emergence of new variant
ag_years <- ag_years %>%
  rowwise() %>%
    mutate(last_possible_year = max(last_possible_year,
                                    last_recorded_epi))

## correct for mislabelling 
end_rows <- which(ag_years$subtype == "H1sea" &
                  ag_years$last_possible_year >= 2009)

ag_years$last_possible_year[end_rows] <- 2008

## generate the full list of possible years for each variant
cumulative_incidence_by_ag <- ag_years %>% ungroup() %>%
    group_by(city, subtype, reference_strain) %>%
    expand(.,
           year=full_seq(c(emerge_year,
                           last_possible_year), 1))

cumulative_incidence_by_ag <- cumulative_incidence_by_ag %>%
    group_by(city, reference_strain) %>%
    mutate(years_since_emergence = year - min(year))

## populate table with epidemic data
cumulative_incidence_by_ag <- cumulative_incidence_by_ag %>%
    left_join(epi_dat,
              by=c('city',
                   'subtype',
                   'reference_strain',
                   'year'))

no_epi <- cumulative_incidence_by_ag$epi_alarm=='N' |
    is.na(cumulative_incidence_by_ag$epi_alarm)

cumulative_incidence_by_ag$incidence_current_season =
    ifelse(
        no_epi,
        0,
        cumulative_incidence_by_ag$incidence_per_mil)


## add cumulative incidence for each variant
cumulative_incidence_by_ag <- cumulative_incidence_by_ag %>%
    group_by(city, subtype, reference_strain) %>%
    arrange(year, .by_group=TRUE) %>%
    mutate(cumulative_prior_incidence_same_ag = (
        cumsum(incidence_current_season) - incidence_current_season))


## normalize cumulative incidence by city-
## and subtype-specific mean epidemic size
cumulative_incidence_by_ag <- cumulative_incidence_by_ag %>%
    dplyr::group_by(city, subtype) %>%
    mutate(mean_epi_size = mean(incidence_per_mil, na.rm=TRUE)) ## want mean size IF an epi happens

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

dat$subtype_id <- as.numeric(
    factor(dat$subtype, 
           levels=unique(dat$subtype)))

epi_variation <- epi_dat %>% group_by(subtype, city) %>%
    summarise(mean_log_inc = mean(log(incidence_per_mil),na.rm=TRUE),
              mean_inc = mean(incidence_per_mil, na.rm=TRUE),
              sd_log_inc = sd(log(incidence_per_mil),na.rm=TRUE))

epi_variation$cv_log_inc <- epi_variation$sd_log_inc / epi_variation$mean_log_inc

## calculate epi_z_score
dat <- dat %>% left_join(epi_variation,
                         by = c('subtype', 'city'))
dat$epi_z_score <- (log(dat$incidence_per_mil) - dat$mean_log_inc) / dat$sd_log_inc

dat$log_relative_epi_size <- log(dat$incidence_per_mil / dat$mean_inc)
    
no_epi <- dat$epi_alarm=='N' |
    is.na(dat$epi_alarm)
dat$epidemic_flag <- ifelse(
    no_epi,
    0,
    1)

## calculate offset from best possible
## start date (the 13th fortnight
## of year, i.e. the hypothetical
## peak of seasonality
## ranges from 0 to 1
dat$start_date_offset <- abs(dat$start - 13) / 13
    

## output only needed columns and rows
columns_wanted = c("city_id",
                   "city",
                   "subtype",
                   "subtype_id",
                   "year",
                   "start",
                   "start_date_offset",
                   "end",
                   "incidence_per_mil",
                   "epidemic_flag",
                   "prior_everything_scaled",
                   "cumulative_prior_incidence_same_ag",
                   "standardised_prior_cumulative",
                   "new_ag_marker",
                   "mean_epi_ah",
                   "mean_epi_temp",
                   "mean_log_inc",
                   "sd_log_inc",
                   "cv_log_inc",
                   "epi_z_score",
                   "log_relative_epi_size")
## data integrity checks
cat("\nChecking data integrity...\n\n")
check_data_integrity(dat)

dat <- dat[, columns_wanted]

dat <- dat %>% arrange(city_id, subtype, year, new_ag_marker)

write.csv(dat,
          file = output_path,
          row.names=FALSE)

warnings()
