library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# loading data ------------------------------------------------------------

epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/epi_table.csv")
mean_fortnightly_climate_30years<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/mean_fortnightly_climate_30years.csv")

epi_table<-epi_table%>%
  dplyr::mutate(scaled_incidence_city = incidence_per_mil/mean_epi_size)

mean_size_subtype_city<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::summarise(mean_epi_size_sc = mean(incidence_per_mil,na.rm=TRUE))

epi_table<-left_join(epi_table,mean_size_subtype_city)
epi_table<-epi_table%>%
  dplyr::mutate(scaled_incidence_subtype_city = incidence_per_mil/mean_epi_size_sc)

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

epi_table_with_clim<-adply(epi_table%>%subset(.,epi_alarm=="Y"),1,mean_climate_over_epi)
epi_table_with_clim<-adply(epi_table_with_clim,1,mean_climate_over_epi,y="temp")

epi_table_with_clim%>%
  subset(.,epi_alarm=="Y")%>%
  lm(scaled_incidence_subtype_city ~ start + log(prior_everything_scaled+0.0001) + as.factor(new_ag_marker) + mean_epi_ah + mean_epi_temp, data=.)%>%
  summary()




# onset -------------------------------------------------------------------
find_preonset_sample<-function(x, n_fortnight = 2){
  x<-as.vector(x)
  temp<-mean_fortnightly_climate_30years%>%
    subset(.,city==as.character(x$city) & year==x$year & fortnights_since_start_of_year%in%c((x$start-n_fortnight):(x$start-1)))
  #print(temp)
  temp<-temp%>%
    dplyr::group_by(city,year)%>%
    dplyr::summarise(start = first(x$start),
                     mean_AH = mean(mean_AH,na.rm=TRUE),
                     mean_temp = mean(mean_temp,na.rm = TRUE),
                     sample_mean_d.AH=mean(d.AH,na.rm=TRUE),
                     #sample_mean_d.SH=mean(d.SH,na.rm=TRUE),
                     sample_mean_d.temp=mean(d.temp,na.rm=TRUE),
                     mean_AH_for_that_fortnight_of_year = mean(mean_AH_for_that_fortnight_of_year,na.rm=TRUE),
                     mean_RH_for_that_fortnight_of_year = mean(mean_RH_for_that_fortnight_of_year,na.rm=TRUE),
                     mean_temp_for_that_fortnight_of_year = mean(mean_temp_for_that_fortnight_of_year,na.rm=TRUE))
  return(data.frame(temp))
}

#function that gets the mean AH' for all fortnights preceding and including onset fortnight in that year.
multi_fortnight_climate_func<-function(x,n_fortnight=1){
  x<-as.data.frame(x)
  fortnights_to_grab<-data.frame(city=as.character(x$city),
                                 year=x$year,
                                 start = seq(n_fortnight+1,x$start,1))
  output<-adply(fortnights_to_grab,1,function(x){find_preonset_sample(x,n_fortnight)})
  output$initiate<-"N"
  output$initiate[nrow(output)]<-"Y"
  return(output)
}

epi_firsts<-epi_table%>%subset(.,first_n_biggest=="Y")

fortnight_before_climate<-adply(epi_firsts,1,function(x){multi_fortnight_climate_func(x,1),})
