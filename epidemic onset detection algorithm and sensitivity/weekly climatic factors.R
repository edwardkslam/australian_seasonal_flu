library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)


# The following code will reproduce the climatic factor analyses discussed in the main text:
# 1)  Comparison of observed climatic fluctuations prior to epidemic onset against typical wintertime
#     fluctuations using the bootstrap sampling method presented by Shaman (2010)
#       - final_results_pooled (Table S1)
#       - final_results (Table S2)
# 2)  Comparison of observed climatic fluctuations prior to epidemic onset against historical averages
#     typical for that time of the year
#       - The output is two plots which are combined to form the figures in the manuscript: T_plot and AH-plot 
#         (aggregated across all cities Figure 2; by city Figure S2)


# Loading data ------------------------------------------------------------

mean_weekly_climate_30years<-read.csv("./dat/raw/mean_weekly_climate_30years.csv")
epi_table<-read.csv("./dat/raw/sensitivity testing/epi_table_weekly.csv")



cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)

# essentially, this contains the epidemic onset timing for the earliest epidemic of each season and city.
first_epi<-epi_table%>%
  subset(.,year!=2009)%>%
  subset(.,first_n_biggest=="Y")


# Comparison against general wintertime conditions Shaman (2010) ----------

# n week block sampling function
# find the mean climatic value for n_week before a certain timepoint

find_preonset_sample<-function(x, n_week = 2){
  x<-as.vector(x)
  temp<-mean_weekly_climate_30years%>%
    subset(.,city==as.character(x$city) & year==x$year & weeks_since_start_of_year%in%c((x$start-n_week):(x$start-1)))
  #print(temp)
  temp<-temp%>%
    dplyr::group_by(city,year)%>%
    dplyr::summarise(start = first(x$start),
                     mean_AH = mean(mean_AH,na.rm=TRUE),
                     mean_temp = mean(mean_temp,na.rm = TRUE),
                     sample_mean_d.AH=mean(d.AH,na.rm=TRUE),
                     #sample_mean_d.SH=mean(d.SH,na.rm=TRUE),
                     sample_mean_d.temp=mean(d.temp,na.rm=TRUE),
                     mean_AH_for_that_week_of_year = mean(mean_AH_for_that_week_of_year,na.rm=TRUE),
                     mean_RH_for_that_week_of_year = mean(mean_RH_for_that_week_of_year,na.rm=TRUE),
                     mean_temp_for_that_week_of_year = mean(mean_temp_for_that_week_of_year,na.rm=TRUE))
  return(data.frame(temp))
}



# mean climatic conditions in the 2-, 4- and 6- weeks prior to empirically observed epidemic onsets
preonset_sample_1wk<-adply(first_epi,1,function(x){find_preonset_sample(x,1)})
preonset_sample_2wk<-adply(first_epi,1,function(x){find_preonset_sample(x,2)})
preonset_sample_3wk<-adply(first_epi,1,function(x){find_preonset_sample(x,3)})


# data frame listing the years and weeks during "winter", from which bootstrap sample will be drawn from
year_week_1wk<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_week_1wk<-year_week_1wk%>%
  subset(.,year!=2009)

year_week_2wk<-expand.grid(year=c(1985:2015),start=c(9:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_week_2wk<-year_week_2wk%>%
  subset(.,year!=2009)

year_week_3wk<-expand.grid(year=c(1985:2015),start=c(10:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_week_3wk<-year_week_3wk%>%
  subset(.,year!=2009)


# bootstrap by city -------------------------------------------------------
# lists for storing bootstrap samples and summary statistics
samples_per_bootstrap_sample<-first_epi$year%>%unique%>%length()
bootstrap_n<-1000000
bootstrap_sample_list_1wk<-list()
bootstrap_sample_list_2wk<-list()
bootstrap_sample_list_3wk<-list()

bootstrap_results_1wk<-list()
bootstrap_results_2wk<-list()
bootstrap_results_3wk<-list()



for(i in 1: length(cities)){
  print(cities[i])
  min_possible_year<-mean_weekly_climate_30years%>%
    subset(.,city==cities[i])%>%
    dplyr::group_by(year)%>%
    dplyr::summarise(n=n())
  min_possible_year<-min_possible_year$year[min(which(min_possible_year$n==52))]
  
  # mean of 2-week block samples
  sample_means_1wk<-adply(year_week_1wk%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,1)})
  
  
  temp<-sample_n(sample_means_1wk,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_1wk[[i]]<-temp
  
  empirical_values_1wk<-preonset_sample_1wk%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  bootstrap_results_1wk[[i]]<-data.frame(city = cities[i],
                                          mean_emperical_sample_mean_d.AH = empirical_values_1wk$empirical_mean_d.AH,
                                          d.AH_pvalue = sum(empirical_values_1wk$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                          
                                          mean_emperical_sample_mean_d.temp = empirical_values_1wk$empirical_mean_d.temp,
                                          d.temp_pvalue = sum(empirical_values_1wk$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n
  )
  
  # mean of 4-week block samples
  sample_means_2wk<-adply(year_week_2wk%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,2)})
  
  temp<-sample_n(sample_means_2wk,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_2wk[[i]]<-temp
  
  empirical_values_2wk<-preonset_sample_2wk%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  
  bootstrap_results_2wk[[i]]<-data.frame(city = cities[i],
                                          mean_emperical_sample_mean_d.AH = empirical_values_2wk$empirical_mean_d.AH,
                                          d.AH_pvalue = sum(empirical_values_2wk$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                          
                                          mean_emperical_sample_mean_d.temp = empirical_values_2wk$empirical_mean_d.temp,
                                          d.temp_pvalue = sum(empirical_values_2wk$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n
  )
  
  # mean of 6-week block samples
  sample_means_3wk<-adply(year_week_3wk%>%
                             subset(.,year>=min_possible_year)%>%
                             dplyr::mutate(city=cities[i]),1,
                           function(x){find_preonset_sample(x,3)})
  
  temp<-sample_n(sample_means_3wk,bootstrap_n*samples_per_bootstrap_sample,replace=TRUE)%>%
    #mean by every 15 rows to produce the 100kx 15
    dplyr::group_by(G=trunc(samples_per_bootstrap_sample:(n()+samples_per_bootstrap_sample-1)/samples_per_bootstrap_sample))%>%
    dplyr::summarise(bootstrap_sample_mean_d.AH = mean(sample_mean_d.AH,na.rm=TRUE),
                     bootstrap_sample_mean_d.temp = mean(sample_mean_d.temp,na.rm=TRUE))
  
  bootstrap_sample_list_3wk[[i]]<-temp
  
  empirical_values_3wk<-preonset_sample_3wk%>%subset(.,city==cities[i])%>%
    dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                     empirical_mean_d.temp = mean(sample_mean_d.temp))
  
  bootstrap_results_3wk[[i]]<-data.frame(city = cities[i],
                                          mean_emperical_sample_mean_d.AH = empirical_values_3wk$empirical_mean_d.AH,
                                          d.AH_pvalue = sum(empirical_values_3wk$empirical_mean_d.AH>=temp$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                          
                                          mean_emperical_sample_mean_d.temp = empirical_values_3wk$empirical_mean_d.temp,
                                          d.temp_pvalue = sum(empirical_values_3wk$empirical_mean_d.temp>=temp$bootstrap_sample_mean_d.temp)/bootstrap_n
  )
  
}

# converting the lists to data frames and adding column denoting whether it is 2-/4-/6- week block samples
bootstrap_results_1wk<-ldply(bootstrap_results_1wk)
bootstrap_results_1wk$num_preceding_week <- 1

bootstrap_results_2wk<-ldply(bootstrap_results_2wk)
bootstrap_results_2wk$num_preceding_week <- 2

bootstrap_results_3wk<-ldply(bootstrap_results_3wk)
bootstrap_results_3wk$num_preceding_week <- 3


# Bootstrap results
final_results<-rbind(bootstrap_results_1wk,bootstrap_results_2wk,bootstrap_results_3wk)
print(final_results)

final_results<-final_results%>%mutate(.,
                                      ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                      ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


# bootstrap pooling across all cities -------------------------------------
pooled_bootstrap_sample_1wk<-list()
pooled_bootstrap_sample_2wk<-list()
pooled_bootstrap_sample_3wk<-list()

#subsample 1/5 of each city specific bootstrap distribution and then combine it together so that overall bootstrap n = 100,000
for(i in 1:length(cities)){
  pooled_bootstrap_sample_1wk[[i]]<-sample_frac(bootstrap_sample_list_1wk[[i]],1/length(cities))
  pooled_bootstrap_sample_2wk[[i]]<-sample_frac(bootstrap_sample_list_2wk[[i]],1/length(cities))
  pooled_bootstrap_sample_3wk[[i]]<-sample_frac(bootstrap_sample_list_3wk[[i]],1/length(cities))
}

pooled_bootstrap_sample_1wk<-ldply(pooled_bootstrap_sample_1wk)
pooled_bootstrap_sample_2wk<-ldply(pooled_bootstrap_sample_2wk)
pooled_bootstrap_sample_3wk<-ldply(pooled_bootstrap_sample_3wk)

empirical_sample_pooled_1wk<-preonset_sample_1wk%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_2wk<-preonset_sample_2wk%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

empirical_sample_pooled_3wk<-preonset_sample_3wk%>%
  dplyr::summarise(empirical_mean_d.AH = mean(sample_mean_d.AH),
                   empirical_mean_d.temp = mean(sample_mean_d.temp))

pooled_results_1wk<-data.frame(num_preceding_week = 1,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_1wk$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_1wk$empirical_mean_d.AH>=pooled_bootstrap_sample_1wk$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_1wk$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_1wk$empirical_mean_d.temp>=pooled_bootstrap_sample_1wk$bootstrap_sample_mean_d.temp)/bootstrap_n)

pooled_results_2wk<-data.frame(num_preceding_week = 2,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_2wk$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_2wk$empirical_mean_d.AH>=pooled_bootstrap_sample_2wk$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_2wk$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_2wk$empirical_mean_d.temp>=pooled_bootstrap_sample_2wk$bootstrap_sample_mean_d.temp)/bootstrap_n)


pooled_results_3wk<-data.frame(num_preceding_week = 3,
                                mean_emperical_sample_mean_d.AH = empirical_sample_pooled_3wk$empirical_mean_d.AH,
                                d.AH_pvalue = sum(empirical_sample_pooled_3wk$empirical_mean_d.AH>=pooled_bootstrap_sample_3wk$bootstrap_sample_mean_d.AH)/bootstrap_n,
                                
                                mean_emperical_sample_mean_d.temp = empirical_sample_pooled_3wk$empirical_mean_d.temp,
                                d.temp_pvalue = sum(empirical_sample_pooled_3wk$empirical_mean_d.temp>=pooled_bootstrap_sample_3wk$bootstrap_sample_mean_d.temp)/bootstrap_n)


final_results_pooled<-rbind(pooled_results_1wk,pooled_results_2wk,pooled_results_3wk)
final_results_pooled<-final_results_pooled%>%mutate(.,
                                                    ManuscriptTable_Temp = paste(signif(mean_emperical_sample_mean_d.temp,3),"  (",signif(d.temp_pvalue,3),")",sep=""),
                                                    ManuscriptTable_AH = paste(signif(mean_emperical_sample_mean_d.AH,3),"  (",signif(d.AH_pvalue,3),")",sep=""))


print(final_results_pooled)



# saving bootstrap result tables ------------------------------------------
write.csv(final_results_pooled,"./nat comms reviewers comments/climate sensitivity/table_S1_weekly.csv",row.names = FALSE)
write.csv(final_results,"./nat comms reviewers comments/climate sensitivity/table_S2_weekly.csv",row.names = FALSE)



################################

# Comparison against historic climatic values of that particular time of year --------

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
  
  temp_clim<-mean_weekly_climate_30years%>%
    subset(.,city== as.character(x$city) & year==x$year)
  temp_clim<-temp_clim[order(temp_clim$year,temp_clim$weeks_since_start_of_year),]
  
  temp_clim<-temp_clim%>%
    dplyr::group_by(year)%>%
    dplyr::mutate(f2f_mean_temp = mean_temp-lag(mean_temp,1),
                  f2f_mean_AH = mean_AH-lag(mean_AH,1))
  
  temp_clim<-temp_clim%>%subset(.,weeks_since_start_of_year %in% (c(y[1]:y[2])+x$start))
  
  temp_clim<-temp_clim%>%
    dplyr::mutate(relative_week= seq_along(weeks_since_start_of_year) - 0.5 -length(weeks_since_start_of_year)/2)
  
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


# data manipulation
# climatic values in the five fortnights before and after each epidemic onset
centered_df<-adply(first_epi,1,climate_centered,.expand=FALSE,.id = NULL)

#mean climatic values for the five weeks before and after epidemic onset
mean_stats_aggregated<-centered_df%>%
  dplyr::group_by(relative_week)%>%
  dplyr::summarise(mean_temp=mean(mean_temp,na.rm=TRUE),
                   mean_d.temp=mean(d.temp,na.rm=TRUE),
                   sd_d.temp=sd(d.temp,na.rm=TRUE),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH,na.rm=TRUE),
                   mean_d.AH=mean(d.AH,na.rm=TRUE),
                   sd_d.AH=sd(d.AH,na.rm=TRUE),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_week_of_year = mean(mean_AH_for_that_week_of_year,na.rm=TRUE),
                   mean_temp_for_that_week_of_year = mean(mean_temp_for_that_week_of_year,na.rm=TRUE))

mean_stats_aggregated<-mean_stats_aggregated%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_week_of_year = mean_relative_humidity_calc(mean_AH_for_that_week_of_year,mean_temp_for_that_week_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_week_of_year)


# Figure 2: plots showing T' and AH' in the 5x 2-week periods before and after epidemic onset-------------


T_plot<-mean_stats_aggregated%>%
  ggplot(data=.,aes(x=relative_week,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_week,
                    ymin=mean_d.temp-sd_d.temp,
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
             aes(x=relative_week,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-5,7.5,0.5), limits = c(-5,7.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  theme_bw()+
  xlab("Weeks relative to onset")+
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

AH_plot<-mean_stats_aggregated%>%
  ggplot(data=.,aes(x=relative_week,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_week,
                    ymin=mean_d.AH-sd_d.AH,
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
             aes(x=relative_week,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-4,4.5,0.5), limits = c(-4.2,4.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Weeks relative to onset")+
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


# city by city climate ----------------------------------------------------
mean_stats_city<-centered_df%>%
  dplyr::group_by(city,relative_week)%>%
  dplyr::summarise(mean_temp=mean(mean_temp,na.rm=TRUE),
                   mean_d.temp=mean(d.temp,na.rm=TRUE),
                   sd_d.temp=sd(d.temp,na.rm=TRUE),
                   wt_d.temp=wilcox.test(d.temp,alternative = c("less"))$p.value,
                   
                   mean_AH=mean(mean_AH,na.rm=TRUE),
                   mean_d.AH=mean(d.AH,na.rm=TRUE),
                   sd_d.AH=sd(d.AH,na.rm=TRUE),
                   wt_d.AH=wilcox.test(d.AH,alternative = c("less"))$p.value,
                   
                   mean_AH_for_that_week_of_year = mean(mean_AH_for_that_week_of_year,na.rm=TRUE),
                   mean_temp_for_that_week_of_year = mean(mean_temp_for_that_week_of_year,na.rm=TRUE))

mean_stats_city<-mean_stats_city%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH_for_that_week_of_year = mean_relative_humidity_calc(mean_AH_for_that_week_of_year,mean_temp_for_that_week_of_year),
                mean_d.RH = mean_relative_humidity_calc(mean_AH ,mean_temp) - mean_RH_for_that_week_of_year)

AT_plot2<-mean_stats_city%>%
  ggplot(data=.,aes(x=relative_week,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_week,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
             aes(x=relative_week,
                 y=d.temp),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.temp),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-5,7.5,0.5), limits = c(-5,7.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Weeks relative to onset")+
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
        panel.grid.minor = element_blank())+
  facet_grid(~ as.factor(city),labeller = label_wrap_gen(width=10))

AH_plot2<-mean_stats_city%>%
  ggplot(data=.,aes(x=relative_week,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_week,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
             aes(x=relative_week,
                 y=d.AH),
             position = position_jitter(w = 0.1, h = 0),
             size=3,
             alpha = 0.2)+
  geom_point(aes(colour=signif_neg.d.AH),
             size=4)+
  scale_color_manual(name="Statistically significant decrease temp (<0.05)",
                     values=c("TRUE"="#D55E00",
                              "FALSE"="#56B4E9")) +
  
  scale_y_continuous(breaks=seq(-4,4.5,0.5), limits = c(-4.2,4.5))+
  scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5.5,5.5))+
  
  theme_bw()+
  xlab("Weeks relative to onset")+
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
        panel.grid.minor = element_blank())+
  facet_grid(~ as.factor(city),labeller = label_wrap_gen(width=10))



# save plot ---------------------------------------------------------------
fig2<-grid.arrange(T_plot,AH_plot,ncol=1)
ggsave(plot = fig2,"./nat comms reviewers comments/climate sensitivity/figure_2_weekly.png",
       width=12, height=11,limitsize=FALSE)

fig_S2<-grid.arrange(AT_plot2,AH_plot2,ncol=1)
ggsave(plot = fig_S2,"./nat comms reviewers comments/climate sensitivity/figure_S2_weekly.png",
       width=12, height=11,limitsize=FALSE)
