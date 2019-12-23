#!/usr/bin/env Rscript

#####################################
## create file of epidemic data augmented
## with climate data, for downstream
## analysis
######################################


########################
## functions
#######################

## function to return mean climate values
mean_climate_over_epi <- function(epidemic,
                                  climate_data,
                                  variable = "ah"){
  epidemic <- as.data.frame(epidemic)
  if(epidemic$epi_alarm=="N"){
    if(variable == "ah"){
      return(data.frame(mean_epi_ah=NA))
    }
    
    if(variable == "temp"){
      return(data.frame(mean_epi_temp=NA))
    }
    if(variable == "rainfall"){
      return(data.frame(mean_epi_rainfall=NA))
    }
  }
  
  fortnights <- seq(epidemic$start, epidemic$end, 1)
  temp_clim <- subset(climate_data,
                      city == as.character(epidemic$city))
  temp_clim <- temp_clim %>%
      subset(., (year == epidemic$year &
                 fortnights_since_start_of_year %in% fortnights))
  if(variable == "ah"){
    return(data.frame(mean_epi_ah = mean(temp_clim$mean_AH)))
  }
  if(variable == "temp"){
    return(data.frame(mean_epi_temp=mean(temp_clim$mean_temp)))
  }
  if(variable == "rainfall"){
    return(data.frame(mean_epi_rainfall=mean(temp_clim$mean_rainfall)))
  }
}

#######################
## load needed packages
#######################
script_packages <- c(
    'plyr',       # for adply()
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

##########################
## read command line args,
## load data
##########################
args <- commandArgs(trailingOnly=TRUE)
epi_data_path <- args[1] 
climate_data_path <- args[2]
output_path <- args[3]

## load data
epi_table <- read_csv(epi_data_path,
                    col_types = cols())

climate_data <- read_csv(climate_data_path,
                        col_types = cols())

epi_table$city <- factor(epi_table$city,
                         levels = unique(epi_table$city))

epi_table_with_clim <- epi_table  %>%
    adply(1, mean_climate_over_epi,
          climate_data = climate_data,
          variable = 'ah')

epi_table_with_clim <- epi_table_with_clim %>%
    adply(1, mean_climate_over_epi,
          climate_data = climate_data,
          variable = "temp")

epi_table_with_clim <- epi_table_with_clim %>%
  adply(1, mean_climate_over_epi,
        climate_data = climate_data,
        variable = "rainfall")

write.csv(epi_table_with_clim,
          file = output_path,
          row.names = FALSE)

warnings()
