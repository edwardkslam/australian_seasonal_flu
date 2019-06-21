#!/usr/bin/env Rscript

################################
## fit the specified stan model
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
model_src_path <- "multilevel_incidence_model.stan"
datapath <- "../epi_table.csv"
humid_datapath <- "../mean_fortnightly_climate_30years.csv"

## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 1500
nchains = n_cores
adapt_d= 0.8
max_tree = 15
fixed_seed = 232

## load data
dat <- read_csv(datapath)
humid_dat <- read_csv(humid_datapath)

#############################
## format data for model
############################

library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)
library(readr)

#Here we assume and correct for potential mis-identification during antigenic characterisation due to delays in 
#updating vaccine strain nomenclature.

# The following code will reproduce the prior immunity analyses discussed in the main text:
# 1)  The relationship between epidemic incidence and the amount of antigenic variant-specific cumulative incidence 
#     epi_size_cumulative_size_same_variant_plot (Figure 4)
#
# 2)  The relationship between the probability of successful epidemic initiation 
#     and amount of antigenic variant-specific cumulative incidence for each subtype 
#     prob_successful_epi_cumulative_size_same_variant_plot (Figure S16)
#
# 3)  Binary logistic regression assessing the effect of antigenic variant-specific cumulative incidence 
#     on the probability of successful epidemic initiation for each subtype 
#     subtype_logistics_regression (Table S6)

# Loading data ------------------------------------------------------------

epi_table = read_csv("../epi_table.csv")

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


dodgy_prev<-epi_table%>%subset(.,is.na(new_ag_marker) | 
                                 reference_strain=="A/Perth/16/2009-like")

# First find the first and last year in which an ag variant causes an epidemic in each city
first_emerge_table<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::summarise(emerge_year=min(year),
                   last_recorded_epi=max(year))

# assume that an ag variant could potentially have caused an epidemic right up to the year 
# before the emergence of its replacement new variant
first_emerge_table<-first_emerge_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::arrange(emerge_year,.by_group=TRUE)%>%
  dplyr::mutate(last_possible_year=if(subtype!="H1sea"){lead(emerge_year-1,default = 2015)}else{lead(emerge_year-1,default = 2008)})

# the last possible year in which an epidemic could have been caused is whichever occured last:
# the last recorded epidemic (to account for the case in which old and new variants circulate in same year)
# the year prior to emergence of new variant
first_emerge_table<-first_emerge_table%>%
  dplyr::rowwise()%>%
  dplyr::mutate(last_possible_year=max(last_possible_year,last_recorded_epi))

# correcting for mislabelling 
end_rows<-which(first_emerge_table$subtype=="H1sea" & first_emerge_table$last_possible_year>=2009)
first_emerge_table$last_possible_year[end_rows]<-2008

# generate the full list of possible years for each variant
cumulative_incidence_by_ag<-first_emerge_table%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  tidyr::expand(.,year=full_seq(c(emerge_year,last_possible_year),1))

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,reference_strain)%>%
  dplyr::mutate(years_since_emergence=year-min(year))

# filling in which are the years with epidemics and their sizes
cumulative_incidence_by_ag<-left_join(cumulative_incidence_by_ag,epi_table[,c("city","subtype","reference_strain","year","incidence_per_mil")])
cumulative_incidence_by_ag<-rename(cumulative_incidence_by_ag,replace=c("incidence_per_mil"= "incidence_current_season"))
cumulative_incidence_by_ag<-data.frame(cumulative_incidence_by_ag)
cumulative_incidence_by_ag$incidence_current_season[which(is.na(cumulative_incidence_by_ag$incidence_current_season))]<-0

# add cumulative incidence for each variant
cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::arrange(year,.by_group=TRUE)%>%
  dplyr::mutate(cumulative_prior_incidence_same_ag = cumsum(incidence_current_season)- incidence_current_season)

# need to divide cumulative incidence by city-specific mean epidemic size
overall_mean_size<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  dplyr::group_by(city)%>%
  dplyr::summarise(mean_epi_size=mean(incidence_per_mil))

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%left_join(overall_mean_size)

cumulative_incidence_by_ag$standardised_prior_cumulative<-cumulative_incidence_by_ag$cumulative_prior_incidence_same_ag/cumulative_incidence_by_ag$mean_epi_size
cumulative_incidence_by_ag$standardised_current_season<-cumulative_incidence_by_ag$incidence_current_season/cumulative_incidence_by_ag$mean_epi_size

dat<-cumulative_incidence_by_ag%>%left_join(dat, by=c('city','year','subtype','reference_strain'))

dat<-dat[dat$subtype=="H3" & !is.na(dat$start) & !is.na(dat$incidence_per_mil),]

dat <- dat %>% left_join(humid_dat, by=c('city'='city','year'='year','start'='fortnights_since_start_of_year'))
dat<-dat[dat$subtype=="H3" & !is.na(dat$start) & !is.na(dat$incidence_per_mil),]

dat$city_id <- as.numeric(factor(dat$city, 
                                 levels=unique(dat$city)))


## make data into list
data_list <- list(
    n_cities=length(unique(dat$city)),
    n_epidemics=length(dat$incidence_per_mil),
    incidences=log(dat$incidence_per_mil),
    city=dat$city_id,
    antigenic_change=dat$new_ag_marker,
    abs_humidity=dat$mean_AH,
    prior_activity=dat$prior_everything_scaled,
    cumulative_prior_incidence=dat$standardised_prior_cumulative)

hyperparam_list <- list(
    mean_city_reporting_rates_per_hundred=1
    sd_city_reporting_rates_per_hundred=2,
    alpha_average_epi_attack_rate=2,
    beta_average_epi_attack_rate=22,
    sd_sd_incidences=2)

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

warnings()
