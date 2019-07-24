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

#Here we assume and correct for potential mis-identification during antigenic characterisation due to delays in 
#updating vaccine strain nomenclature.

# The following code will reproduce the multivariate linear regression models in the main text,
# assessing the joint contributions of climatic and virological factors to epidemic incidence.
#
# The output is as follows:
# 1)  Regression coefficients for full and submodels
#     output_table (Table S15)
#
# 2)  RSE, R-squared, adjusted R-squared for full and submodels
#     rsquared_table (Table S16)


# loading data ------------------------------------------------------------
raw_table<-read.csv("./dat/raw/raw_data.csv")
raw_table<-raw_table%>%
  dplyr::mutate(fortnights_since_start_of_year = lubridate::yday(specimen_date)%/%14+1)%>%
  dplyr::group_by(city,year,assumed_antigenic_variant,fortnights_since_start_of_year)%>%
  dplyr::summarise(count=n())

epi_table<-read.csv("./dat/raw/epi_table.csv")
mean_fortnightly_climate_30years<-read.csv("./dat/raw/mean_fortnightly_climate_30years.csv")

epi_table<-epi_table%>%
  dplyr::mutate(scaled_incidence_city = incidence_per_mil/mean_epi_size,
                log_incidence = log(incidence_per_mil))

epi_table<-epi_table%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA))


epi_table<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::mutate(scaled_incidence_subtype_city = log(incidence_per_mil)-mean(log(incidence_per_mil),na.rm=TRUE))

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
epi_table_with_clim_2$first_in_season<-factor(epi_table_with_clim_2$first_in_season)

temp_data<-epi_table_with_clim_2%>%
  subset(.,epi_alarm=="Y" & !is.na(standardised_prior_cumulative))

temp_data<-temp_data%>%dplyr::mutate(pes = ifelse(first_in_season==1,0,log(prior_everything_scaled)),
                                     spc = ifelse(new_ag_marker==1,0,standardised_prior_cumulative))


# lm (mean climatic) ------------------------------------------------------

full_model<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc+ mean_epi_ah + mean_epi_temp,data=.)

summary_fm<-full_model%>%
  summary()

model.1<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc+ mean_epi_ah + mean_epi_temp,data=.)

summary_1<-model.1%>%
  summary()

model.2<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + 
       as.factor(new_ag_marker)+spc+ mean_epi_ah + mean_epi_temp,data=.)

summary_2<-model.2%>%
  summary()

model.3<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       mean_epi_ah + mean_epi_temp,data=.)

summary_3<-model.3%>%
  summary()

model.4<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + mean_epi_temp,data=.)

summary_4<-model.4%>%
  summary()

model.5<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + mean_epi_ah,data=.)

summary_5<-model.5%>%
  summary()

coef_fm<-summary_fm$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_fm$coefficients%>%as.data.frame()),model = "Full Model")
coef_1<-summary_1$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_1$coefficients%>%as.data.frame()),model = "Model 1")
coef_2<-summary_2$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_2$coefficients%>%as.data.frame()),model = "Model 2")
coef_3<-summary_3$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_3$coefficients%>%as.data.frame()),model = "Model 3")
coef_4<-summary_4$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_4$coefficients%>%as.data.frame()),model = "Model 4")
coef_5<-summary_5$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_5$coefficients%>%as.data.frame()),model = "Model 5")

output_table<-rbind(coef_fm,coef_1,coef_2,coef_3,coef_4,coef_5)
col_order<-c("model","predictor","Estimate","Std. Error","t value","Pr(>|t|)")
output_table<-output_table[,col_order]
output_table$predictor<-factor(output_table$predictor)

output_table$predictor<-mapvalues(output_table$predictor,
                                  c("as.factor(first_in_season)1","pes","as.factor(new_ag_marker)1","spc"),
                                  c("first_in_season","prior_activity_other_subtypes|first_in_season==F","new_ag_maker","cumulative_activity_same_variant|new_ag_maker==F"))


rsquared_fm<-data.frame(model = "Full Model",
                        RSE = summary_fm$sigma,
                        Rsquared = summary_fm$r.squared,
                        adjusted_Rsquared = summary_fm$adj.r.squared)
rsquared_1<-data.frame(model = "Model 1",
                       RSE = summary_1$sigma,
                       Rsquared = summary_1$r.squared,
                       adjusted_Rsquared = summary_1$adj.r.squared)
rsquared_2<-data.frame(model = "Model 2",
                       RSE = summary_2$sigma,
                       Rsquared = summary_2$r.squared,
                       adjusted_Rsquared = summary_2$adj.r.squared)
rsquared_3<-data.frame(model = "Model 3",
                       RSE = summary_3$sigma,
                       Rsquared = summary_3$r.squared,
                       adjusted_Rsquared = summary_3$adj.r.squared)
rsquared_4<-data.frame(model = "Model 4",
                       RSE = summary_4$sigma,
                       Rsquared = summary_4$r.squared,
                       adjusted_Rsquared = summary_4$adj.r.squared)
rsquared_5<-data.frame(model = "Model 5",
                       RSE = summary_5$sigma,
                       Rsquared = summary_5$r.squared,
                       adjusted_Rsquared = summary_5$adj.r.squared)

rsquared_table<-rbind(rsquared_fm,rsquared_1,rsquared_2,rsquared_3,rsquared_4,rsquared_5)

# save tables -------------------------------------------------------------

write.csv(output_table%>%dplyr::mutate_if(is.numeric,signif,digits=3),
          "./tables/table_S15.csv",row.names = FALSE)

write.csv(rsquared_table%>%dplyr::mutate_if(is.numeric,signif,digits=3),
          "./tables/table_S16.csv",row.names = FALSE)


stop("main analyses complete")

#continue past if you want to repeat but with mean climatic values from onset to peak

# lm (early climatic) -----------------------------------------------------

full_model<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
         as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

summary_fm<-full_model%>%
  summary()

model.1<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

summary_1<-model.1%>%
  summary()

model.2<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + 
       as.factor(new_ag_marker)+spc+ early_ah + early_temp,data=.)

summary_2<-model.2%>%
  summary()

model.3<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       early_ah + early_temp,data=.)

summary_3<-model.3%>%
  summary()

model.4<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + early_temp,data=.)

summary_4<-model.4%>%
  summary()

model.5<-temp_data%>%
  lm(scaled_incidence_subtype_city ~ start + as.factor(first_in_season) + pes + 
       as.factor(new_ag_marker)+spc + early_ah,data=.)

summary_5<-model.5%>%
  summary()

coef_fm<-summary_fm$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_fm$coefficients%>%as.data.frame()),model = "Full Model")
coef_1<-summary_1$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_1$coefficients%>%as.data.frame()),model = "Model 1")
coef_2<-summary_2$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_2$coefficients%>%as.data.frame()),model = "Model 2")
coef_3<-summary_3$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_3$coefficients%>%as.data.frame()),model = "Model 3")
coef_4<-summary_4$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_4$coefficients%>%as.data.frame()),model = "Model 4")
coef_5<-summary_5$coefficients%>%as.data.frame()%>%dplyr::mutate(predictor = row.names(summary_5$coefficients%>%as.data.frame()),model = "Model 5")

output_table<-rbind(coef_fm,coef_1,coef_2,coef_3,coef_4,coef_5)
col_order<-c("model","predictor","Estimate","Std. Error","t value","Pr(>|t|)")
output_table<-output_table[,col_order]
output_table$predictor<-factor(output_table$predictor)

output_table$predictor<-mapvalues(output_table$predictor,
                                  c("as.factor(first_in_season)1","pes","as.factor(new_ag_marker)1","spc"),
                                  c("first_in_season","prior_activity_other_subtypes|first_in_season==F","new_ag_maker","cumulative_activity_same_variant|new_ag_maker==F"))


rsquared_fm<-data.frame(model = "Full Model",
                        RSE = summary_fm$sigma,
                        Rsquared = summary_fm$r.squared,
                        adjusted_Rsquared = summary_fm$adj.r.squared)
rsquared_1<-data.frame(model = "Model 1",
                       RSE = summary_1$sigma,
                       Rsquared = summary_1$r.squared,
                       adjusted_Rsquared = summary_1$adj.r.squared)
rsquared_2<-data.frame(model = "Model 2",
                       RSE = summary_2$sigma,
                       Rsquared = summary_2$r.squared,
                       adjusted_Rsquared = summary_2$adj.r.squared)
rsquared_3<-data.frame(model = "Model 3",
                       RSE = summary_3$sigma,
                       Rsquared = summary_3$r.squared,
                       adjusted_Rsquared = summary_3$adj.r.squared)
rsquared_4<-data.frame(model = "Model 4",
                       RSE = summary_4$sigma,
                       Rsquared = summary_4$r.squared,
                       adjusted_Rsquared = summary_4$adj.r.squared)
rsquared_5<-data.frame(model = "Model 5",
                       RSE = summary_5$sigma,
                       Rsquared = summary_5$r.squared,
                       adjusted_Rsquared = summary_5$adj.r.squared)

rsquared_table<-rbind(rsquared_fm,rsquared_1,rsquared_2,rsquared_3,rsquared_4,rsquared_5)


# save tables -------------------------------------------------------------

write.csv(output_table%>%dplyr::mutate_if(is.numeric,signif,digits=3),
          "./tables/glm_earlyClim_coef.csv",row.names = FALSE)

write.csv(rsquared_table%>%dplyr::mutate_if(is.numeric,signif,digits=3),
          "./tables/glm_earlyClim_rsquared.csv",row.names = FALSE)
