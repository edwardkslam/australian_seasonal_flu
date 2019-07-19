library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess the robustness of our analyses on climatic fluctuations preceding epidemic onset potential inaccuracies in the estimation of epidemic onset timing.
# We incorporate alternative estimates given by Geoghegan et al. (2018), based on a set of different assumptions numbered 1 to 4.
# We considered whether or not these preceding climatic fluctuations were
#   i)  anomalous when compared with against "typical wintertime" fluctuations using the bootstrap sampling method presented by Shaman (2010). 
#         For each set of assumptions, the output is a data frame (Table S2, S3, S4)
#   ii) anomalous when compared against historical averages that typical for that time of the year. 
#         For each set of assumptions, the output is two plots which are combined to form the figures in the manuscript: T_plot and AH-plot (Figure S4, S5, S6, S7).


########ASSUMPTIONS########

#   1)  For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the DOMINANT influenza A subtype 
#       was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       The timing value used for this set of analyses is stored in geog_epi_table$largest_geog_start
#       Outputs
#         i) largest_geog_results_pooled (Table S3)
#         ii) largest_geog_results(Table S4)
#         iii) stacked_climate_plots1 (Figure S4)
#   2)  For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the EARLIEST influenza A subtype 
#       was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       The timing value used for this set of analyses is stored in geog_epi_table$earliest_geog_start
#       Outputs
#         i) earliest_geog_results_pooled (Table S5)
#         ii) earliest_geog_results (Table S6)
#         ii) stacked_climate_plots2 (Figure S5)
#   3)  For seasons between 2007 and 2015 in which the number of cases for the dominant influenza A subtype were small or 
#       it was difficult to discern the period of epidemic from background activity, 
#       we assumed that our timing estimate was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       Outputs
#         i) poor_geog_results_pooled (Table S7)
#         ii) poor_geog_results (Table S8)
#         ii) stacked_climate_plots3 (Figure S6)
#   4)  Utilising only the timing estimates by Geoghegan et al. (2018), we assessed if more generally, the onset of influenza A epidemic activity 
#       in the seasons from 2007 to 2015 was preceded by periods of anomalous climatic conditions.
#         i) just_geog_results_pooled (Table S9)
#         ii) just_geog_results (Table S10)
#         ii) stacked_climate_plots3 (Figure S7)

# loading in data ---------------------------------------------------------


mean_fortnightly_climate_30years<-read.csv("./dat/raw/mean_fortnightly_climate_30years.csv")
geog_epi_table<-read.csv("./dat/raw/geoghegan/epi_table_with_geoghegan_estimates.csv")
just_geog_estimates<-read.csv("./dat/raw/geoghegan/Geoghegan_2018_estimated_A_onsets.csv")%>%
  subset(.,year%in%c(2007:2015))

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")
geog_epi_table$city<-factor(geog_epi_table$city,levels = cities)

#subsetting to retain only earliest (based on various assumptions) and largest epidemic of each season and city
largest_use_geog_earliest<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & year!=2009)%>%
  dplyr::group_by(city,year)%>%dplyr::summarise(start=min(largest_geog_start))

earliest_use_geog_earliest<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & year!=2009)%>%
  dplyr::group_by(city,year)%>%dplyr::summarise(start=min(earliest_geog_start))

poor_use_geog_earliest<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & year!=2009)%>%
  dplyr::group_by(city,year)%>%dplyr::summarise(start=min(poor_timeseries_geog_start))


# Comparison against general wintertime conditions Shaman (2010) ----------

# n week block sampling function
# find the mean climatic value for n_fortnight before a certain timepoint

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


# 1) i) largest use Geoghegan ; compare against "typical wintertime conditions" ------------------------------------
largest_use_geog_sample1<-adply(largest_use_geog_earliest,1,function(x){find_preonset_sample(x,1)})
largest_use_geog_sample2<-adply(largest_use_geog_earliest,1,function(x){find_preonset_sample(x,2)})
largest_use_geog_sample3<-adply(largest_use_geog_earliest,1,function(x){find_preonset_sample(x,3)})

#quick calculator of d.AH -> d.RH
#largest_use_geog_sample1%>%dplyr::group_by(city)%>%dplyr::summarise(mean_d.RH = mean_relative_humidity_calc(mean(mean_AH),mean(mean_temp))-mean_relative_humidity_calc(mean(mean_AH)+mean(sample_mean_d.AH),mean(mean_temp)))

# data frame listing the years and fortnights during "winter", from which bootstrap sample will be drawn from 
year_fortnight_1ftn<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_1ftn<-year_fortnight_1ftn%>%
  subset(.,year!=2009)

year_fortnight_2ftn<-expand.grid(year=c(1985:2015),start=c(9:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_2ftn<-year_fortnight_2ftn%>%
  subset(.,year!=2009)

year_fortnight_3ftn<-expand.grid(year=c(1985:2015),start=c(10:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_3ftn<-year_fortnight_3ftn%>%
  subset(.,year!=2009)

samples_per_bootstrap_sample<-largest_use_geog_earliest$year%>%unique%>%length()
bootstrap_n<-1000000
bootstrap_sample_list_1ftn<-list()
bootstrap_sample_list_2ftn<-list()
bootstrap_sample_list_3ftn<-list()

largest_bootstrap_ttest_results_1ftn<-list()
largest_bootstrap_ttest_results_2ftn<-list()
largest_bootstrap_ttest_results_3ftn<-list()

for(i in 1: length(cities)){
  min_possible_year<-mean_fortnightly_climate_30years%>%
    subset(.,city==cities[i])%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(n=n())
  min_possible_year<-min_possible_year$year[min(which(min_possible_year$n==26))]
  
  # mean of 2-week block samples
  sample_means_1ftn<-adply(year_fortnight_1ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,1)})

  temp<-sample_n(sample_means_1ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_1ftn[[i]]<-temp
  
  empirical_values_1ftn<-largest_use_geog_sample1%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  largest_bootstrap_ttest_results_1ftn[[i]]<-data.frame(city = cities[i],
                                                mean_emperical_sample_mean_d.AH = empirical_values_1ftn$empirical_mean_d.AH,
                                                d.AH_pvalue = sum(empirical_values_1ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                
                                                mean_emperical_sample_mean_d.temp = empirical_values_1ftn$empirical_mean_d.temp,
                                                d.temp_pvalue = sum(empirical_values_1ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                num_preceding_fortnight = 1
                                                )
  
  
  
  
  # mean of 4-week block samples
  sample_means_2ftn<-adply(year_fortnight_2ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,2)})
  
  temp<-sample_n(sample_means_2ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_2ftn[[i]]<-temp
  
  empirical_values_2ftn<-largest_use_geog_sample2%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  largest_bootstrap_ttest_results_2ftn[[i]]<-data.frame(city = cities[i],
                                                        mean_emperical_sample_mean_d.AH = empirical_values_2ftn$empirical_mean_d.AH,
                                                        d.AH_pvalue = sum(empirical_values_2ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                        
                                                        mean_emperical_sample_mean_d.temp = empirical_values_2ftn$empirical_mean_d.temp,
                                                        d.temp_pvalue = sum(empirical_values_2ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                        num_preceding_fortnight = 2
  )
  # mean of 6-week block samples
  sample_means_3ftn<-adply(year_fortnight_3ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,3)})
  
  temp<-sample_n(sample_means_3ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_3ftn[[i]]<-temp
  
  empirical_values_3ftn<-largest_use_geog_sample3%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  largest_bootstrap_ttest_results_3ftn[[i]]<-data.frame(city = cities[i],
                                                        mean_emperical_sample_mean_d.AH = empirical_values_3ftn$empirical_mean_d.AH,
                                                        d.AH_pvalue = sum(empirical_values_3ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                        
                                                        mean_emperical_sample_mean_d.temp = empirical_values_3ftn$empirical_mean_d.temp,
                                                        d.temp_pvalue = sum(empirical_values_3ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                        num_preceding_fortnight = 3
  )
}


largest_geog_results<-rbind(ldply(largest_bootstrap_ttest_results_1ftn),
                            ldply(largest_bootstrap_ttest_results_2ftn),
                            ldply(largest_bootstrap_ttest_results_3ftn))

largest_geog_results<-largest_geog_results%>%
  dplyr::group_by(num_preceding_fortnight)%>%
  dplyr::mutate(adjusted_d.AH_pvalue = p.adjust(d.AH_pvalue,method = "holm"),
                adjusted_d.temp_pvalue = p.adjust(d.temp_pvalue,method = "holm"))

largest_geog_results<-largest_geog_results%>%
  mutate(.,
         ManuscriptTable_Temp = paste(round(mean_emperical_sample_mean_d.temp,3),"  (",round(d.temp_pvalue,3),")",sep=""),
         ManuscriptTable_AH = paste(round(mean_emperical_sample_mean_d.AH,3),"  (",round(d.AH_pvalue,3),")",sep=""))


#pooled
pooled_bootstrap_sample_1ftn<-list()
pooled_bootstrap_sample_2ftn<-list()
pooled_bootstrap_sample_3ftn<-list()

#subsample 1/5 of each city specific bootstrap distribution and then combine it together so that overall bootstrap n = 100,000
for(i in 1:length(cities)){
  pooled_bootstrap_sample_1ftn[[i]]<-sample_frac(bootstrap_sample_list_1ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_2ftn[[i]]<-sample_frac(bootstrap_sample_list_2ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_3ftn[[i]]<-sample_frac(bootstrap_sample_list_3ftn[[i]],1/length(cities))
}

pooled_bootstrap_sample_1ftn<-ldply(pooled_bootstrap_sample_1ftn)
pooled_bootstrap_sample_2ftn<-ldply(pooled_bootstrap_sample_2ftn)
pooled_bootstrap_sample_3ftn<-ldply(pooled_bootstrap_sample_3ftn)

empirical_sample_pooled_1ftn<-largest_use_geog_sample1%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_2ftn<-largest_use_geog_sample2%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_3ftn<-largest_use_geog_sample3%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

pooled_results_1ftn<-data.frame(num_preceding_fortnight = 1,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_1ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_1ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

pooled_results_2ftn<-data.frame(num_preceding_fortnight = 2,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_2ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_2ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)


pooled_results_3ftn<-data.frame(num_preceding_fortnight = 3,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_3ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_3ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

largest_geog_results_pooled<-rbind(pooled_results_1ftn,pooled_results_2ftn,pooled_results_3ftn)
largest_geog_results_pooled<-largest_geog_results_pooled%>%mutate(.,
                                                    ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                                    ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


print(largest_geog_results_pooled)






# 2) i) earliest use Geoghegan ; compare against "typical wintertime conditions" ------------------------------------
earliest_use_geog_sample1<-adply(earliest_use_geog_earliest,1,function(x){find_preonset_sample(x,1)})
earliest_use_geog_sample2<-adply(earliest_use_geog_earliest,1,function(x){find_preonset_sample(x,2)})
earliest_use_geog_sample3<-adply(earliest_use_geog_earliest,1,function(x){find_preonset_sample(x,3)})

# data frame listing the years and fortnights during "winter", from which bootstrap sample will be drawn from 
year_fortnight_1ftn<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_1ftn<-year_fortnight_1ftn%>%
  subset(.,year!=2009)

year_fortnight_2ftn<-expand.grid(year=c(1985:2015),start=c(9:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_2ftn<-year_fortnight_2ftn%>%
  subset(.,year!=2009)

year_fortnight_3ftn<-expand.grid(year=c(1985:2015),start=c(10:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_3ftn<-year_fortnight_3ftn%>%
  subset(.,year!=2009)

samples_per_bootstrap_sample<-earliest_use_geog_earliest$year%>%unique%>%length()
bootstrap_n<-1000000
earliest_bootstrap_ttest_results_1ftn<-list()
earliest_bootstrap_ttest_results_2ftn<-list()
earliest_bootstrap_ttest_results_3ftn<-list()
for(i in 1: length(cities)){
  min_possible_year<-mean_fortnightly_climate_30years%>%
    subset(.,city==cities[i])%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(n=n())
  min_possible_year<-min_possible_year$year[min(which(min_possible_year$n==26))]
  
  # mean of 2-week block samples
  sample_means_1ftn<-adply(year_fortnight_1ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,1)})
  
  temp<-sample_n(sample_means_1ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_1ftn[[i]]<-temp
  
  empirical_values_1ftn<-earliest_use_geog_sample1%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  earliest_bootstrap_ttest_results_1ftn[[i]]<-data.frame(city = cities[i],
                                                         mean_emperical_sample_mean_d.AH = empirical_values_1ftn$empirical_mean_d.AH,
                                                         d.AH_pvalue = sum(empirical_values_1ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                         
                                                         mean_emperical_sample_mean_d.temp = empirical_values_1ftn$empirical_mean_d.temp,
                                                         d.temp_pvalue = sum(empirical_values_1ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                         num_preceding_fortnight = 1
  )
  
  
  
  
  # mean of 4-week block samples
  sample_means_2ftn<-adply(year_fortnight_2ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,2)})
  
  temp<-sample_n(sample_means_2ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_2ftn[[i]]<-temp
  
  empirical_values_2ftn<-earliest_use_geog_sample2%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  earliest_bootstrap_ttest_results_2ftn[[i]]<-data.frame(city = cities[i],
                                                         mean_emperical_sample_mean_d.AH = empirical_values_2ftn$empirical_mean_d.AH,
                                                         d.AH_pvalue = sum(empirical_values_2ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                         
                                                         mean_emperical_sample_mean_d.temp = empirical_values_2ftn$empirical_mean_d.temp,
                                                         d.temp_pvalue = sum(empirical_values_2ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                         num_preceding_fortnight = 2
  )
  # mean of 6-week block samples
  sample_means_3ftn<-adply(year_fortnight_3ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,3)})
  
  temp<-sample_n(sample_means_3ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_3ftn[[i]]<-temp
  
  empirical_values_3ftn<-earliest_use_geog_sample3%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  earliest_bootstrap_ttest_results_3ftn[[i]]<-data.frame(city = cities[i],
                                                         mean_emperical_sample_mean_d.AH = empirical_values_3ftn$empirical_mean_d.AH,
                                                         d.AH_pvalue = sum(empirical_values_3ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                         
                                                         mean_emperical_sample_mean_d.temp = empirical_values_3ftn$empirical_mean_d.temp,
                                                         d.temp_pvalue = sum(empirical_values_3ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                         num_preceding_fortnight = 3
  )
}


earliest_geog_results<-rbind(ldply(earliest_bootstrap_ttest_results_1ftn),
                             ldply(earliest_bootstrap_ttest_results_2ftn),
                             ldply(earliest_bootstrap_ttest_results_3ftn))

earliest_geog_results<-earliest_geog_results%>%
  dplyr::group_by(num_preceding_fortnight)%>%
  dplyr::mutate(adjusted_d.AH_pvalue = p.adjust(d.AH_pvalue,method = "holm"),
                adjusted_d.temp_pvalue = p.adjust(d.temp_pvalue,method = "holm"))

earliest_geog_results<-earliest_geog_results%>%
  mutate(.,
         ManuscriptTable_Temp = paste(round(mean_emperical_sample_mean_d.temp,3),"  (",round(d.temp_pvalue,3),")",sep=""),
         ManuscriptTable_AH = paste(round(mean_emperical_sample_mean_d.AH,3),"  (",round(d.AH_pvalue,3),")",sep=""))


#pooled
pooled_bootstrap_sample_1ftn<-list()
pooled_bootstrap_sample_2ftn<-list()
pooled_bootstrap_sample_3ftn<-list()

#subsample 1/5 of each city specific bootstrap distribution and then combine it together so that overall bootstrap n = 100,000
for(i in 1:length(cities)){
  pooled_bootstrap_sample_1ftn[[i]]<-sample_frac(bootstrap_sample_list_1ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_2ftn[[i]]<-sample_frac(bootstrap_sample_list_2ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_3ftn[[i]]<-sample_frac(bootstrap_sample_list_3ftn[[i]],1/length(cities))
}

pooled_bootstrap_sample_1ftn<-ldply(pooled_bootstrap_sample_1ftn)
pooled_bootstrap_sample_2ftn<-ldply(pooled_bootstrap_sample_2ftn)
pooled_bootstrap_sample_3ftn<-ldply(pooled_bootstrap_sample_3ftn)

empirical_sample_pooled_1ftn<-earliest_use_geog_sample1%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_2ftn<-earliest_use_geog_sample2%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_3ftn<-earliest_use_geog_sample3%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

pooled_results_1ftn<-data.frame(num_preceding_fortnight = 1,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_1ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_1ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

pooled_results_2ftn<-data.frame(num_preceding_fortnight = 2,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_2ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_2ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)


pooled_results_3ftn<-data.frame(num_preceding_fortnight = 3,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_3ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_3ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

earliest_geog_results_pooled<-rbind(pooled_results_1ftn,pooled_results_2ftn,pooled_results_3ftn)
earliest_geog_results_pooled<-earliest_geog_results_pooled%>%mutate(.,
                                                                  ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                                                  ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


print(earliest_geog_results_pooled)


# 3) i) replace poorly defined timings with Geoghegan ; compare against "typical wintertime conditions"------------------------------------
poor_use_geog_sample1<-adply(poor_use_geog_earliest,1,function(x){find_preonset_sample(x,1)})
poor_use_geog_sample2<-adply(poor_use_geog_earliest,1,function(x){find_preonset_sample(x,2)})
poor_use_geog_sample3<-adply(poor_use_geog_earliest,1,function(x){find_preonset_sample(x,3)})

# data frame listing the years and fortnights during "winter", from which bootstrap sample will be drawn from 
year_fortnight_1ftn<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_1ftn<-year_fortnight_1ftn%>%
  subset(.,year!=2009)

year_fortnight_2ftn<-expand.grid(year=c(1985:2015),start=c(9:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_2ftn<-year_fortnight_2ftn%>%
  subset(.,year!=2009)

year_fortnight_3ftn<-expand.grid(year=c(1985:2015),start=c(10:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_3ftn<-year_fortnight_3ftn%>%
  subset(.,year!=2009)

samples_per_bootstrap_sample<-poor_use_geog_earliest$year%>%unique%>%length()
bootstrap_n<-1000000
poor_bootstrap_ttest_results_1ftn<-list()
poor_bootstrap_ttest_results_2ftn<-list()
poor_bootstrap_ttest_results_3ftn<-list()


for(i in 1: length(cities)){
  min_possible_year<-mean_fortnightly_climate_30years%>%
    subset(.,city==cities[i])%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(n=n())
  min_possible_year<-min_possible_year$year[min(which(min_possible_year$n==26))]
  
  # mean of 2-week block samples
  sample_means_1ftn<-adply(year_fortnight_1ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,1)})
  
  temp<-sample_n(sample_means_1ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_1ftn[[i]]<-temp
  
  empirical_values_1ftn<-poor_use_geog_sample1%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  poor_bootstrap_ttest_results_1ftn[[i]]<-data.frame(city = cities[i],
                                                     mean_emperical_sample_mean_d.AH = empirical_values_1ftn$empirical_mean_d.AH,
                                                     d.AH_pvalue = sum(empirical_values_1ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                     
                                                     mean_emperical_sample_mean_d.temp = empirical_values_1ftn$empirical_mean_d.temp,
                                                     d.temp_pvalue = sum(empirical_values_1ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                     num_preceding_fortnight = 1
  )
  
  
  
  
  # mean of 4-week block samples
  sample_means_2ftn<-adply(year_fortnight_2ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,2)})
  
  temp<-sample_n(sample_means_2ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_2ftn[[i]]<-temp
  
  empirical_values_2ftn<-poor_use_geog_sample2%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  poor_bootstrap_ttest_results_2ftn[[i]]<-data.frame(city = cities[i],
                                                     mean_emperical_sample_mean_d.AH = empirical_values_2ftn$empirical_mean_d.AH,
                                                     d.AH_pvalue = sum(empirical_values_2ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                     
                                                     mean_emperical_sample_mean_d.temp = empirical_values_2ftn$empirical_mean_d.temp,
                                                     d.temp_pvalue = sum(empirical_values_2ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                     num_preceding_fortnight = 2
  )
  # mean of 6-week block samples
  sample_means_3ftn<-adply(year_fortnight_3ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,3)})
  
  temp<-sample_n(sample_means_3ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_3ftn[[i]]<-temp
  
  empirical_values_3ftn<-poor_use_geog_sample3%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  poor_bootstrap_ttest_results_3ftn[[i]]<-data.frame(city = cities[i],
                                                     mean_emperical_sample_mean_d.AH = empirical_values_3ftn$empirical_mean_d.AH,
                                                     d.AH_pvalue = sum(empirical_values_3ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                     
                                                     mean_emperical_sample_mean_d.temp = empirical_values_3ftn$empirical_mean_d.temp,
                                                     d.temp_pvalue = sum(empirical_values_3ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                     num_preceding_fortnight = 3
  )
}

poor_geog_results<-rbind(ldply(poor_bootstrap_ttest_results_1ftn),
                         ldply(poor_bootstrap_ttest_results_2ftn),
                         ldply(poor_bootstrap_ttest_results_3ftn))

poor_geog_results<-poor_geog_results%>%
  dplyr::group_by(num_preceding_fortnight)%>%
  dplyr::mutate(adjusted_d.AH_pvalue = p.adjust(d.AH_pvalue,method = "holm"),
                adjusted_d.temp_pvalue = p.adjust(d.temp_pvalue,method = "holm"))

poor_geog_results<-poor_geog_results%>%
  mutate(.,
         ManuscriptTable_Temp = paste(round(mean_emperical_sample_mean_d.temp,3),"  (",round(d.temp_pvalue,3),")",sep=""),
         ManuscriptTable_AH = paste(round(mean_emperical_sample_mean_d.AH,3),"  (",round(d.AH_pvalue,3),")",sep=""))

#pooled
pooled_bootstrap_sample_1ftn<-list()
pooled_bootstrap_sample_2ftn<-list()
pooled_bootstrap_sample_3ftn<-list()

#subsample 1/5 of each city specific bootstrap distribution and then combine it together so that overall bootstrap n = 100,000
for(i in 1:length(cities)){
  pooled_bootstrap_sample_1ftn[[i]]<-sample_frac(bootstrap_sample_list_1ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_2ftn[[i]]<-sample_frac(bootstrap_sample_list_2ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_3ftn[[i]]<-sample_frac(bootstrap_sample_list_3ftn[[i]],1/length(cities))
}

pooled_bootstrap_sample_1ftn<-ldply(pooled_bootstrap_sample_1ftn)
pooled_bootstrap_sample_2ftn<-ldply(pooled_bootstrap_sample_2ftn)
pooled_bootstrap_sample_3ftn<-ldply(pooled_bootstrap_sample_3ftn)

empirical_sample_pooled_1ftn<-poor_use_geog_sample1%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_2ftn<-poor_use_geog_sample2%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_3ftn<-poor_use_geog_sample3%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

pooled_results_1ftn<-data.frame(num_preceding_fortnight = 1,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_1ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_1ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

pooled_results_2ftn<-data.frame(num_preceding_fortnight = 2,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_2ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_2ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)


pooled_results_3ftn<-data.frame(num_preceding_fortnight = 3,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_3ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_3ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

poor_geog_results_pooled<-rbind(pooled_results_1ftn,pooled_results_2ftn,pooled_results_3ftn)
poor_geog_results_pooled<-poor_geog_results_pooled%>%mutate(.,
                                                                  ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                                                  ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


print(poor_geog_results_pooled)



# 4) i) Just Geoghegan ; compare against "typical wintertime conditions" ------------------------------------
just_geog_sample1<-adply(just_geog_estimates%>%subset(.,year!=2009),1,function(x){find_preonset_sample(x,1)})
just_geog_sample2<-adply(just_geog_estimates%>%subset(.,year!=2009),1,function(x){find_preonset_sample(x,2)})
just_geog_sample3<-adply(just_geog_estimates%>%subset(.,year!=2009),1,function(x){find_preonset_sample(x,3)})

# data frame listing the years and fortnights during "winter", from which bootstrap sample will be drawn from 
year_fortnight_1ftn<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_1ftn<-year_fortnight_1ftn%>%
  subset(.,year!=2009)

year_fortnight_2ftn<-expand.grid(year=c(1985:2015),start=c(9:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_2ftn<-year_fortnight_2ftn%>%
  subset(.,year!=2009)

year_fortnight_3ftn<-expand.grid(year=c(1985:2015),start=c(10:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_3ftn<-year_fortnight_3ftn%>%
  subset(.,year!=2009)

samples_per_bootstrap_sample<-just_geog_estimates%>%subset(.,year!=2009)%>%.$year%>%unique%>%length()
bootstrap_n<-1000000
just_geog_bootstrap_ttest_results_1ftn<-list()
just_geog_bootstrap_ttest_results_2ftn<-list()
just_geog_bootstrap_ttest_results_3ftn<-list()

for(i in 1: length(cities)){
  min_possible_year<-mean_fortnightly_climate_30years%>%
    subset(.,city==cities[i])%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(n=n())
  min_possible_year<-min_possible_year$year[min(which(min_possible_year$n==26))]
  
  # mean of 2-week block samples
  sample_means_1ftn<-adply(year_fortnight_1ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,1)})
  
  temp<-sample_n(sample_means_1ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_1ftn[[i]]<-temp
  
  empirical_values_1ftn<-just_geog_sample1%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  just_geog_bootstrap_ttest_results_1ftn[[i]]<-data.frame(city = cities[i],
                                                          mean_emperical_sample_mean_d.AH = empirical_values_1ftn$empirical_mean_d.AH,
                                                          d.AH_pvalue = sum(empirical_values_1ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                          
                                                          mean_emperical_sample_mean_d.temp = empirical_values_1ftn$empirical_mean_d.temp,
                                                          d.temp_pvalue = sum(empirical_values_1ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                          num_preceding_fortnight = 1
  )
  
  
  
  
  # mean of 4-week block samples
  sample_means_2ftn<-adply(year_fortnight_2ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,2)})
  
  temp<-sample_n(sample_means_2ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_2ftn[[i]]<-temp
  
  empirical_values_2ftn<-just_geog_sample2%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  just_geog_bootstrap_ttest_results_2ftn[[i]]<-data.frame(city = cities[i],
                                                          mean_emperical_sample_mean_d.AH = empirical_values_2ftn$empirical_mean_d.AH,
                                                          d.AH_pvalue = sum(empirical_values_2ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                          
                                                          mean_emperical_sample_mean_d.temp = empirical_values_2ftn$empirical_mean_d.temp,
                                                          d.temp_pvalue = sum(empirical_values_2ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                          num_preceding_fortnight = 2
  )
  # mean of 6-week block samples
  sample_means_3ftn<-adply(year_fortnight_3ftn%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,3)})
  
  temp<-sample_n(sample_means_3ftn,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_3ftn[[i]]<-temp
  
  empirical_values_3ftn<-just_geog_sample3%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  just_geog_bootstrap_ttest_results_3ftn[[i]]<-data.frame(city = cities[i],
                                                          mean_emperical_sample_mean_d.AH = empirical_values_3ftn$empirical_mean_d.AH,
                                                          d.AH_pvalue = sum(empirical_values_3ftn$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                                          
                                                          mean_emperical_sample_mean_d.temp = empirical_values_3ftn$empirical_mean_d.temp,
                                                          d.temp_pvalue = sum(empirical_values_3ftn$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n,
                                                          num_preceding_fortnight = 3
  )
}


just_geog_results<-rbind(ldply(just_geog_bootstrap_ttest_results_1ftn),
                         ldply(just_geog_bootstrap_ttest_results_2ftn),
                         ldply(just_geog_bootstrap_ttest_results_3ftn))

just_geog_results<-just_geog_results%>%
  dplyr::group_by(num_preceding_fortnight)%>%
  dplyr::mutate(adjusted_d.AH_pvalue = p.adjust(d.AH_pvalue,method = "holm"),
                adjusted_d.temp_pvalue = p.adjust(d.temp_pvalue,method = "holm"))

just_geog_results<-just_geog_results%>%
  mutate(.,
         ManuscriptTable_Temp = paste(round(mean_emperical_sample_mean_d.temp,3),"  (",round(d.temp_pvalue,3),")",sep=""),
         ManuscriptTable_AH = paste(round(mean_emperical_sample_mean_d.AH,3),"  (",round(d.AH_pvalue,3),")",sep=""))

#pooled
pooled_bootstrap_sample_1ftn<-list()
pooled_bootstrap_sample_2ftn<-list()
pooled_bootstrap_sample_3ftn<-list()

#subsample 1/5 of each city specific bootstrap distribution and then combine it together so that overall bootstrap n = 100,000
for(i in 1:length(cities)){
  pooled_bootstrap_sample_1ftn[[i]]<-sample_frac(bootstrap_sample_list_1ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_2ftn[[i]]<-sample_frac(bootstrap_sample_list_2ftn[[i]],1/length(cities))
  pooled_bootstrap_sample_3ftn[[i]]<-sample_frac(bootstrap_sample_list_3ftn[[i]],1/length(cities))
}

pooled_bootstrap_sample_1ftn<-ldply(pooled_bootstrap_sample_1ftn)
pooled_bootstrap_sample_2ftn<-ldply(pooled_bootstrap_sample_2ftn)
pooled_bootstrap_sample_3ftn<-ldply(pooled_bootstrap_sample_3ftn)

empirical_sample_pooled_1ftn<-just_geog_sample1%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_2ftn<-just_geog_sample2%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_3ftn<-just_geog_sample3%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

pooled_results_1ftn<-data.frame(num_preceding_fortnight = 1,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_1ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_1ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_1ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_1ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

pooled_results_2ftn<-data.frame(num_preceding_fortnight = 2,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_2ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_2ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_2ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_2ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)


pooled_results_3ftn<-data.frame(num_preceding_fortnight = 3,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_3ftn$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.AH>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_3ftn$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_3ftn$empirical_mean_d.temp>=pooled_bootstrap_sample_3ftn$bootstrap_sample_mean_d.temp)/bootstrap_n)

just_geog_geog_results_pooled<-rbind(pooled_results_1ftn,pooled_results_2ftn,pooled_results_3ftn)
just_geog_geog_results_pooled<-just_geog_geog_results_pooled%>%mutate(.,
                                                                  ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                                                  ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


print(just_geog_geog_results_pooled)


# saving bootstrap result tables ------------------------------------------

write.csv(largest_geog_results_pooled,"./tables/table_S3.csv",row.names = FALSE)
write.csv(largest_geog_results,"./tables/table_S4.csv",row.names = FALSE)

write.csv(earliest_geog_results_pooled,"./tables/table_S5.csv",row.names = FALSE)
write.csv(earliest_geog_results,"./tables/table_S6.csv",row.names = FALSE)

write.csv(poor_geog_results_pooled,"./tables/table_S7.csv",row.names = FALSE)
write.csv(poor_geog_results,"./tables/table_S8.csv",row.names = FALSE)

write.csv(just_geog_geog_results_pooled,"./tables/table_S9.csv",row.names = FALSE)
write.csv(just_geog_results,"./tables/table_S10.csv",row.names = FALSE)




########################################################################################

# Comparison against historic climatic values of that particular time of year --------------------------------------------------------------

# Defining functions 

# Takes an entry from epidemic table and spits out the climatic values 5 fortnights either side of onset (by default)
# mean_AH is the year specific value
# mean_AH_for_that_fortnight_of_year etc is the 31 year mean value for that particular fortnight of the year
# d.AH is mean_AH - mean_AH_for_that_fortnight_of_year
# f2f_mean_AH is the difference in mean AH between current fortnight and previous one
# delay is copied from epidemic table and is to denote delay in onset relative to first epidemic of that season in that city
# strain_year denotes current year and the strain

climate_centered<-function(x,y=c(-5,5),point="start"){
  x<-as.data.frame(x)
  if(point=="start"){
    reference<-x$start
  }
  #for use in replicates
  if(point=="earliest_geog"){
    reference<-x$earliest_geog_start
  }
  if(point =="largest_geog"){
    reference<-x$largest_geog_start
  }
  if(point=="poor_time_series_geog"){
    reference<-x$poor_timeseries_geog_start
  }
  if(point=="just_geog"){
    reference<-x$start
  }
  
  temp_clim<-mean_fortnightly_climate_30years%>%
    subset(.,city== as.character(x$city) & year==x$year)
  temp_clim<-temp_clim[order(temp_clim$year,temp_clim$fortnights_since_start_of_year),]
  
  temp_clim<-temp_clim%>%
    dplyr::group_by(year)%>%
    dplyr::mutate(f2f_mean_temp = mean_temp-lag(mean_temp,1),
                  f2f_mean_AH = mean_AH-lag(mean_AH,1))
  
  temp_clim<-temp_clim%>%subset(.,fortnights_since_start_of_year %in% (c(y[1]:y[2])+x$start))
  
  temp_clim<-temp_clim%>%
    dplyr::mutate(relative_fortnight = seq_along(fortnights_since_start_of_year) - 0.5 -length(fortnights_since_start_of_year)/2)
  
  temp_clim$delay <- x$delay
  temp_clim$strain_year<-x$strain_year
  
  return(data.frame(temp_clim))
  
}

# once the mean AH and temp for the fortnights before/after epidemic onset
# have been calculated, function for calculating equivalent RH %
# 
mean_relative_humidity_calc<-function(mean_ah,mean_temp){
  #https://carnotcycle.wordpress.com/2012/08/04/how-to-convert-relative-humidity-to-absolute-humidity/
  #http://bioma.jrc.ec.europa.eu/components/componentstools/evapotranspiration/help/Actual_vapor_pressure.html
  #http://www.reahvac.com/tools/humidity-formulas/
  denom <-6.112 * exp(17.67*mean_temp / (mean_temp+243.5)) *  2.1674
  ans<-mean_ah * (mean_temp+273.15)/denom
  return(ans)
}



# 1) ii) largest use Geoghegan ; compare against "historic climatic values of that particular time of year" ------------------------------------------------------------------
centered_df1<-adply(geog_epi_table%>%subset(.,firstNbiggest_largest_geog=="Y" & year!=2009),1,climate_centered,point="largest_geog",.expand=FALSE,.id = NULL)

mean_stats1<-centered_df1%>%
  dplyr::group_by(relative_fortnight)%>%
  dplyr::summarise(mean_temp=mean(mean_temp),
                   mean_d.temp=mean(d.temp),
                   sd_d.temp=sd(d.temp),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH),
                   mean_d.AH=mean(d.AH),
                   sd_d.AH=sd(d.AH),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_fortnight_of_year = mean(mean_AH_for_that_fortnight_of_year),
                   mean_temp_for_that_fortnight_of_year = mean(mean_temp_for_that_fortnight_of_year))

mean_stats1<-mean_stats1%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_fortnight_of_year = mean_relative_humidity_calc(mean_AH_for_that_fortnight_of_year,mean_temp_for_that_fortnight_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_fortnight_of_year)


AT_plot1<-mean_stats1%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df1,
             aes(x=relative_fortnight,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,5,0.5), limits = c(-3.5,5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Temperature ( ",degree,"C)")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

AH_plot1<-mean_stats1%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df1,
             aes(x=relative_fortnight,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,3.5,0.5), limits = c(-3.5,3.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Absolute Humidity "," (g/",m^{3},")",sep="")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig_S4<-grid.arrange(AT_plot1,AH_plot1,ncol=1)


# 2) ii) earliest use Geoghegan ; compare against "historic climatic values of that particular time of year" -----------------------------------------------------------------
centered_df2<-adply(geog_epi_table%>%subset(.,firstNbiggest_earliest_geog=="Y" & year!=2009),1,climate_centered,point="earliest_geog",.expand=FALSE,.id = NULL)

mean_stats2<-centered_df2%>%
  dplyr::group_by(relative_fortnight)%>%
  dplyr::summarise(mean_temp=mean(mean_temp),
                   mean_d.temp=mean(d.temp),
                   sd_d.temp=sd(d.temp),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH),
                   mean_d.AH=mean(d.AH),
                   sd_d.AH=sd(d.AH),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_fortnight_of_year = mean(mean_AH_for_that_fortnight_of_year),
                   mean_temp_for_that_fortnight_of_year = mean(mean_temp_for_that_fortnight_of_year))

mean_stats2<-mean_stats2%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_fortnight_of_year = mean_relative_humidity_calc(mean_AH_for_that_fortnight_of_year,mean_temp_for_that_fortnight_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_fortnight_of_year)

AT_plot2<-mean_stats2%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df2,
             aes(x=relative_fortnight,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,5,0.5), limits = c(-3.5,5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Temperature ( ",degree,"C)")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

AH_plot2<-mean_stats2%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df2,
             aes(x=relative_fortnight,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,3.5,0.5), limits = c(-3.5,3.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Absolute Humidity "," (g/",m^{3},")",sep="")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig_S5<-grid.arrange(AT_plot2,AH_plot2,ncol=1)



# 3) ii) replace poorly defined timings with Geoghegan ; compare against "historic climatic values of that particular time of year"  -----------------------------------------------------------------
centered_df3<-adply(geog_epi_table%>%subset(.,firstNbiggest_poor_geog=="Y" & year!=2009),1,climate_centered,point="poor_time_series_geog",.expand=FALSE,.id = NULL)


mean_stats3<-centered_df3%>%
  dplyr::group_by(relative_fortnight)%>%
  dplyr::summarise(mean_temp=mean(mean_temp),
                   mean_d.temp=mean(d.temp),
                   sd_d.temp=sd(d.temp),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH),
                   mean_d.AH=mean(d.AH),
                   sd_d.AH=sd(d.AH),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_fortnight_of_year = mean(mean_AH_for_that_fortnight_of_year),
                   mean_temp_for_that_fortnight_of_year = mean(mean_temp_for_that_fortnight_of_year))

mean_stats3<-mean_stats3%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_fortnight_of_year = mean_relative_humidity_calc(mean_AH_for_that_fortnight_of_year,mean_temp_for_that_fortnight_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_fortnight_of_year)

AT_plot3<-mean_stats3%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df3,
             aes(x=relative_fortnight,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,5,0.5), limits = c(-3.5,5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Temperature ( ",degree,"C)")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

AH_plot3<-mean_stats3%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df3,
             aes(x=relative_fortnight,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,3.5,0.5), limits = c(-3.5,3.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Absolute Humidity "," (g/",m^{3},")",sep="")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig_S6<-grid.arrange(AT_plot3,AH_plot3,ncol=1)


# 4) ii) Just Geoghegan ; compare against "historic climatic values of that particular time of year"  -------------------------------------------------
centered_df4<-adply(just_geog_estimates%>%subset(.,year!=2009),1,climate_centered,point="start",.expand=FALSE,.id=NULL)

mean_stats4<-centered_df4%>%
  dplyr::group_by(relative_fortnight)%>%
  dplyr::summarise(mean_temp=mean(mean_temp),
                   mean_d.temp=mean(d.temp),
                   sd_d.temp=sd(d.temp),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH),
                   mean_d.AH=mean(d.AH),
                   sd_d.AH=sd(d.AH),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_fortnight_of_year = mean(mean_AH_for_that_fortnight_of_year),
                   mean_temp_for_that_fortnight_of_year = mean(mean_temp_for_that_fortnight_of_year))

mean_stats4<-mean_stats4%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_fortnight_of_year = mean_relative_humidity_calc(mean_AH_for_that_fortnight_of_year,mean_temp_for_that_fortnight_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_fortnight_of_year)

AT_plot4<-mean_stats4%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df4,
             aes(x=relative_fortnight,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,5,0.5), limits = c(-3.5,5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Temperature ( ",degree,"C)")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

AH_plot4<-mean_stats4%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df4,
             aes(x=relative_fortnight,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-3.5,3.5,0.5), limits = c(-3.5,3.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Two week intervals relative to onset")+
  ylab(expression(paste("Anomalous Absolute Humidity "," (g/",m^{3},")",sep="")))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig_S7<-grid.arrange(AT_plot4,AH_plot4,ncol=1)


# save plots --------------------------------------------------------------

base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/figures/reviewer comments/"

ggsave(plot = fig_S4,filename = "./figures/supp/figure_S4.png", 
       width=12, height=11,limitsize=FALSE)
ggsave(plot = fig_S5,filename = "./figures/supp/figure_S5.png", 
       width=12, height=11,limitsize=FALSE)
ggsave(plot = fig_S6,filename = "./figures/supp/figure_S6.png", 
       width=12, height=11,limitsize=FALSE)
ggsave(plot = fig_S7,filename = "./figures/supp/figure_S7.png", 
       width=12, height=11,limitsize=FALSE)

