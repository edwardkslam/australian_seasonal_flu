library(lubridate)
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
raw_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/raw_data.csv")
raw_table<-raw_table%>%
  dplyr::mutate(fortnights_since_start_of_year = lubridate::yday(specimen_date)%/%14+1)%>%
  dplyr::group_by(city,year,assumed_antigenic_variant,fortnights_since_start_of_year)%>%
  dplyr::summarise(count=n())

epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/epi_table.csv")
mean_fortnightly_climate_30years<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/mean_fortnightly_climate_30years.csv")

epi_table<-epi_table%>%
  dplyr::mutate(scaled_incidence_city = incidence_per_mil/mean_epi_size,
                log_incidence = log(incidence_per_mil))

epi_table<-epi_table%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA))

mean_size_subtype_city<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::summarise(mean_epi_size_sc = mean(incidence_per_mil,na.rm=TRUE))

epi_table<-left_join(epi_table,mean_size_subtype_city)
epi_table<-epi_table%>%
  dplyr::mutate(scaled_incidence_subtype_city = incidence_per_mil/mean_epi_size_sc)

epi_table<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::mutate(incidence_z_score_subtype_city = ifelse(epi_alarm=="Y",
                                                        (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                        NA))

# function to return mean climate values ----------------------------------
mean_climate_over_epi<-function(x,y="ah"){
  #x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    return(data.frame(x,
                      mean_epi_ah=NA,
                      mean_epi_temp=NA))
  }
  fortnights<-seq(x$start,x$end,1)
  temp_clim<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  temp_clim<-temp_clim%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  temp_clim<-temp_clim%>%dplyr::summarise(mean_epi_ah = mean(mean_AH),
                                   mean_epi_temp = mean(mean_temp))
  return(data.frame(x,temp_clim))
}

early_climate<-function(x,y="ah"){
  #mean climate over start to peak fortnight
  x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    return(data.frame(x,
                      early_ah=NA,
                      early_temp=NA))
  }
  
  temp<-raw_table%>%subset(.,city==x$city & year==x$year & assumed_antigenic_variant==as.character(x$reference_strain) &
                             fortnights_since_start_of_year %in%c(x$start+1:x$end))
  peak_fortnight<-temp$fortnights_since_start_of_year[which.max(temp$count)]
  
  fortnights<-seq(x$start,peak_fortnight,1)
  temp_clim<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  temp_clim<-temp_clim%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  temp_clim<-temp_clim%>%dplyr::summarise(early_ah = mean(mean_AH),
                                          early_temp = mean(mean_temp))

}

epi_table_with_clim<-adply(epi_table%>%subset(.,epi_alarm=="Y"),1,mean_climate_over_epi)

epi_table_with_clim<-adply(epi_table_with_clim,1,early_climate)

# prior immunity  for each ag variant -------------------------------------
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

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%left_join(data.frame(city=epi_table$city,reference_strain=epi_table$reference_strain,year=epi_table$year,epi_alarm=factor(epi_table$epi_alarm,levels = c("Y","N"))))
cumulative_incidence_by_ag$epi_alarm[is.na(cumulative_incidence_by_ag$epi_alarm)==TRUE]<-"N"

cumulative_incidence_by_ag$epi_alarm2<-cumulative_incidence_by_ag$epi_alarm%>%
  mapvalues(.,from=c("Y","N"),to=c(1,0))

# append epi_table with cumulative incidence ------------------------------

epi_table_with_clim_2<-left_join(epi_table_with_clim,cumulative_incidence_by_ag[c("city","subtype","reference_strain","year","standardised_prior_cumulative")])



# adding "interaction" between prior everything and first in season -------------------------------------

epi_table_with_clim_2<-epi_table_with_clim_2%>%
  dplyr::mutate(first_in_season = ifelse(delay==0,1,0))

temp_data<-epi_table_with_clim_2%>%
  subset(.,epi_alarm=="Y" & !is.na(standardised_prior_cumulative))

temp_data<-temp_data%>%dplyr::mutate(pes = ifelse(first_in_season==1,0,log(prior_everything_scaled)),
                                     spc = ifelse(new_ag_marker==1,0,standardised_prior_cumulative))

full_model<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ start + as.factor(first_in_season) + pes + 
         as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

full_model%>%
  summary()

model.1<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ start + 
       as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

model.1%>%
  summary()

model.2<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

model.2%>%
  summary()

model.3<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ start + as.factor(first_in_season) + pes + 
       early_ah + early_temp,data=.)

model.3%>%
  summary()

model.4<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + early_temp,data=.)

model.4%>%
  summary()

model.5<-temp_data%>%
  lm(log(scaled_incidence_subtype_city) ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + early_ah,data=.)

model.5%>%
  summary()










# onset -------------------------------------------------------------------
find_preonset_sample<-function(x, n_fortnight = 1){
  x<-as.vector(x)
  temp<-mean_fortnightly_climate_30years%>%
    subset(.,city==as.character(x$city) & year==x$year & fortnights_since_start_of_year%in%c((x$start-n_fortnight):(x$start-1)))
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
  output$initiate<-0
  output$initiate[nrow(output)]<-1
  return(output)
}

reformat_glm_func<-function(x){
  coef_temp<-x%>%
    coef(.)%>%.[2:3]%>%exp(.)%>%melt(.)
  
  temp_SE<-x%>%
    coef(.)%>%.[2:3,2]%>%melt(.)
    
  temp_p_value<-x%>%
    coef(.)%>%.[2:3,4]%>%melt(.)
  
  coef_temp<-cbind(coef_temp,temp_SE,temp_p_value)
  colnames(coef_temp)<-c("OR","SE","p_value")
  
  return(coef_temp)
}


