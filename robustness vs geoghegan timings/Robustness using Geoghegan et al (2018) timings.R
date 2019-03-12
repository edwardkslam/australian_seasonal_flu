library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)


# loading in data ---------------------------------------------------------

if(Sys.info()['sysname']=="Windows"){
  #   1)  For each season and city which Geoghegan has data for (2007-2015), 
  #       assume that I have misidentified start timing for the LARGEST Influenza Type A epidemic.
  #       and replace it with their timing 
  largest_use_geog<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/largest_use_geog_timing_epi_table.csv")
  
  #   2)  For each season and city which Geoghegan has data for (2007-2015), 
  #       assume that I have misidentified start timing for the EARLIEST ONSET Influenza Type A epidemic.
  #       and replace it with their timing 
  earliest_use_geog<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/earliest_use_geog_timing_epi_table.csv")
  
  #   3)  For each season and city which Geoghegan has data for (2007-2015), 
  #       identify the epidemics within my data set that have poorly defined epidemic time series (ie absence of clear "exponential growth"
  #       assume that I have misidentified start timing for these Influenza Type A epidemics.
  #       and replace it with their timing 
  poorly_defined_used_geog<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/my_dodgy_time_series_use_geog_timing_epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  #   1)  For each season and city which Geoghegan has data for (2007-2015), 
  #       assume that I have misidentified start timing for the LARGEST Influenza Type A epidemic.
  #       and replace it with their timing 
  largest_use_geog<-read.csv("~/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/largest_use_geog_timing_epi_table.csv")
  
  #   2)  For each season and city which Geoghegan has data for (2007-2015), 
  #       assume that I have misidentified start timing for the EARLIEST ONSET Influenza Type A epidemic.
  #       and replace it with their timing 
  earliest_use_geog<-read.csv("~/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/earliest_use_geog_timing_epi_table.csv")
  
  #   3)  For each season and city which Geoghegan has data for (2007-2015), 
  #       identify the epidemics within my data set that have poorly defined epidemic time series (ie absence of clear "exponential growth"
  #       assume that I have misidentified start timing for these Influenza Type A epidemics.
  #       and replace it with their timing 
  poorly_defined_used_geog<-read.csv("~/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/my_dodgy_time_series_use_geog_timing_epi_table.csv")
}


# 1) Ag change and timing (largest_use_geog) ------------------------------


