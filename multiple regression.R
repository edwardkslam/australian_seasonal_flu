library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess whether or not climatic factors act more generally to enhance transmission, rather than as specific triggers
# for epidemic onset
# Figure S3

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  mean_fortnightly_climate_30years<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  mean_fortnightly_climate_30years<-read.csv("~/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table.csv")
}

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


# function to return mean climate values ----------------------------------
mean_climate_over_epi<-function(x,y="ah"){
  x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    if(y=="ah"){
      return(data.frame(mean_epi_ah=NA))
    }
    if(y=="temp"){
      return(data.frame(mean_epi_temp=NA))
    }
  }
  fortnights<-seq(x$start,x$end,1)
  temp_clim<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  temp_clim<-temp_clim%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  if(y=="ah"){
    return(data.frame(mean_epi_ah=mean(temp_clim$mean_AH)))
  }
  if(y=="temp"){
    return(data.frame(mean_epi_temp=mean(temp_clim$mean_temp)))
  }
}


# getting climate for each epidemic ---------------------------------------
epi_table_with_clim<-adply(epi_table%>%subset(.,epi_alarm=="Y"),1,mean_climate_over_epi)
epi_table_with_clim<-adply(epi_table_with_clim,1,mean_climate_over_epi,y="temp")

epi_table_with_clim<-epi_table_with_clim%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(standardised_mean_ah = (mean_epi_ah-mean(mean_epi_ah))/sd(mean_epi_ah),
                standardised_mean_temp = (mean_epi_temp-mean(mean_epi_temp))/sd(mean_epi_temp))

# glm ---------------------------------------------------------------------
temp_data<-epi_table_with_clim %>%
  subset(.,year!=2009 & epi_alarm=="Y" & !is.na(new_ag_marker))%>%
  dplyr::mutate(epi_scaled = incidence_per_mil/mean_epi_size,
                log_prior = log(prior_everything_scaled + 0.001),
                log_relative = log(relative_to_first_n_biggest + 0.001))

temp<-lm(formula= relative_to_first_n_biggest ~ prior_everything_scaled+delay+start+factor(new_ag_marker)+standardised_mean_ah+standardised_mean_temp,
         data=temp_data)

lm(formula= log_relative ~ log_prior+start,data=temp_data)


temp<-lm(formula= epi_scaled ~ prior_everything_scaled+delay+start+factor(new_ag_marker)+standardised_mean_ah+standardised_mean_temp,
         data=temp_data)

epi_table_with_clim %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  dplyr::mutate(log_prior = log(prior_everything_scaled + 0.01),
                log_relative = log(relative_to_first_n_biggest + 0.01))%>%
  lm(formula= log_relative ~ log_prior+start, data=.)

temp<-epi_table_with_clim %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  dplyr::mutate(log_prior = log(prior_everything_scaled + 0.01),
                log_relative = log(relative_to_first_n_biggest + 0.01))%>%
  lm(formula= log_relative ~ log_prior+start+factor(new_ag_marker)+standardised_mean_ah+standardised_mean_temp, data=.)

temp2<-epi_table_with_clim %>%
  subset(.,year!=2009 & epi_alarm=="Y" & delay!=0)%>%
  dplyr::mutate(log_prior = log(prior_everything_scaled + 0.01),
                log_relative = log(relative_to_first_n_biggest + 0.01))%>%
  lm(formula= log_relative ~ log_prior+start+standardised_mean_ah+standardised_mean_temp, data=.)


plot(temp_data$mean_epi_temp,temp_data$relative_to_first_n_biggest)
