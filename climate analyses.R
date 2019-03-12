library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  mean_fortnightly_climate_30years<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  mean_fortnightly_climate_30years<-read.csv("~/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table.csv")
}

cities<-c("BRISBANE","PERTH","ADELAIDE","SYDNEY","MELBOURNE")

first_epi<-epi_table%>%
  subset(.,year!=2009)%>%
  subset(.,first_n_biggest=="Y")


# Shaman (2010) bootstrap method ------------------------------------------
# n week block sampling function
#find the mean climatic value for n_fortnight before a certain timepoint
find_preonset_sample<-function(x, n_fortnight = 2){
  x<-as.vector(x)
  temp<-mean_fortnightly_climate_30years%>%
    subset(.,city==as.character(x$city) & year==x$year & fortnights_since_start_of_year%in%c((x$start-n_fortnight):(x$start-1)))
  #print(temp)
  temp<-temp%>%
    dplyr::group_by(city,year)%>%
    dplyr::summarise(start = first(x$start),
                     sample_mean_d.AH=mean(d.AH,na.rm=TRUE),
                     #sample_mean_d.SH=mean(d.SH,na.rm=TRUE),
                     sample_mean_d.temp=mean(d.temp,na.rm=TRUE))
  return(data.frame(temp))
}

preonset_sample_1ftn<-adply(first_epi,1,function(x){find_preonset_sample(x,1)})
preonset_sample_2ftn<-adply(first_epi,1,function(x){find_preonset_sample(x)})
preonset_sample_3ftn<-adply(first_epi,1,function(x){find_preonset_sample(x,3)})

bootstrap_n<-1000000
bootstrap_sample_list_1ftn<-list()
bootstrap_sample_list_2ftn<-list()
bootstrap_sample_list_3ftn<-list()

bootstrap_ttest_results_1ftn<-list()
bootstrap_ttest_results_2ftn<-list()
bootstrap_ttest_results_3ftn<-list()

year_fortnight_1ftn<-expand.grid(year=c(1985:2015),start=c(8:17))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_1ftn<-year_fortnight_1ftn%>%
  subset(.,year!=2009)

year_fortnight_2ftn<-expand.grid(year=c(1985:2015),start=c(9:18))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_2ftn<-year_fortnight_2ftn%>%
  subset(.,year!=2009)

year_fortnight_3ftn<-expand.grid(year=c(1985:2015),start=c(10:19))    # start of 7 = 01 April ; end of 16 = 31 August
year_fortnight_3ftn<-year_fortnight_3ftn%>%
  subset(.,year!=2009)


# Bootstrap sampling
stopCluster(cl)
start_cl<-Sys.time()
no_cores<-detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)

for(i in 1: length(cities)){
  print(paste("Bootstrapping // Current City = ",cities[i]," // Time elapsed = ", difftime(Sys.time(),start_cl,units="mins"), " mins",sep=""))
  samples_to_select_1ftn<-sample_n(year_fortnight_1ftn,bootstrap_n,replace=TRUE)
  samples_to_select_1ftn$city<-cities[i]
  
  bootstrap_sample_list_1ftn[[i]]<-adply(samples_to_select_1ftn,1,function(x){find_preonset_sample(x,1)},
                                         .parallel = TRUE,
                                         .paropts = list(.export=c("find_preonset_sample","mean_fortnightly_climate_30years"),
                                                         .packages=(.packages()))
  )
  
  bootstrap_ttest_results_1ftn[[i]]<-data.frame(city = cities[i],
                                                mean_synthetic_sample_mean_d.AH = bootstrap_sample_list_1ftn[[i]]$sample_mean_d.AH %>% mean(),
                                                mean_emperical_sample_mean_d.AH = preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH %>%mean(),
                                                d.AH_stat = t.test(bootstrap_sample_list_1ftn[[i]]$sample_mean_d.AH,preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$statistic,
                                                d.AH_pvalue = t.test(bootstrap_sample_list_1ftn[[i]]$sample_mean_d.AH,preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$p.value,
                                                mean_synthetic_sample_mean_d.temp = bootstrap_sample_list_1ftn[[i]]$sample_mean_d.temp %>% mean(),
                                                mean_emperical_sample_mean_d.temp = preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp %>%mean(),
                                                d.temp_stat = t.test(bootstrap_sample_list_1ftn[[i]]$sample_mean_d.temp,preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$statistic,
                                                d.temp_pvalue = t.test(bootstrap_sample_list_1ftn[[i]]$sample_mean_d.temp,preonset_sample_1ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$p.value
  )
  
  samples_to_select_2ftn<-sample_n(year_fortnight_2ftn,bootstrap_n,replace=TRUE)
  samples_to_select_2ftn$city<-cities[i]
  
  bootstrap_sample_list_2ftn[[i]]<-adply(samples_to_select_2ftn,1,function(x){find_preonset_sample(x,2)},
                                         .parallel = TRUE,
                                         .paropts = list(.export=c("find_preonset_sample","mean_fortnightly_climate_30years"),
                                                         .packages=(.packages()))
  )
  
  bootstrap_ttest_results_2ftn[[i]]<-data.frame(city = cities[i],
                                                mean_synthetic_sample_mean_d.AH = bootstrap_sample_list_2ftn[[i]]$sample_mean_d.AH %>% mean(),
                                                mean_emperical_sample_mean_d.AH = preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH %>%mean(),
                                                d.AH_stat = t.test(bootstrap_sample_list_2ftn[[i]]$sample_mean_d.AH,preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$statistic,
                                                d.AH_pvalue = t.test(bootstrap_sample_list_2ftn[[i]]$sample_mean_d.AH,preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$p.value,
                                                mean_synthetic_sample_mean_d.temp = bootstrap_sample_list_2ftn[[i]]$sample_mean_d.temp %>% mean(),
                                                mean_emperical_sample_mean_d.temp = preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp %>%mean(),
                                                d.temp_stat = t.test(bootstrap_sample_list_2ftn[[i]]$sample_mean_d.temp,preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$statistic,
                                                d.temp_pvalue = t.test(bootstrap_sample_list_2ftn[[i]]$sample_mean_d.temp,preonset_sample_2ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$p.value
  )
  
  samples_to_select_3ftn<-sample_n(year_fortnight_3ftn,bootstrap_n,replace=TRUE)
  samples_to_select_3ftn$city<-cities[i]
  
  bootstrap_sample_list_3ftn[[i]]<-adply(samples_to_select_3ftn,1,function(x){find_preonset_sample(x,3)},
                                         .parallel = TRUE,
                                         .paropts = list(.export=c("find_preonset_sample","mean_fortnightly_climate_30years"),
                                                         .packages=(.packages()))
  )
  
  bootstrap_ttest_results_3ftn[[i]]<-data.frame(city = cities[i],
                                                mean_synthetic_sample_mean_d.AH = bootstrap_sample_list_3ftn[[i]]$sample_mean_d.AH %>% mean(),
                                                mean_emperical_sample_mean_d.AH = preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH %>%mean(),
                                                d.AH_stat = t.test(bootstrap_sample_list_3ftn[[i]]$sample_mean_d.AH,preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$statistic,
                                                d.AH_pvalue = t.test(bootstrap_sample_list_3ftn[[i]]$sample_mean_d.AH,preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.AH)%>%.$p.value,
                                                mean_synthetic_sample_mean_d.temp = bootstrap_sample_list_3ftn[[i]]$sample_mean_d.temp %>% mean(),
                                                mean_emperical_sample_mean_d.temp = preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp %>%mean(),
                                                d.temp_stat = t.test(bootstrap_sample_list_3ftn[[i]]$sample_mean_d.temp,preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$statistic,
                                                d.temp_pvalue = t.test(bootstrap_sample_list_3ftn[[i]]$sample_mean_d.temp,preonset_sample_3ftn%>%subset(.,city==cities[i])%>%.$sample_mean_d.temp)%>%.$p.value
  )
}

bootstrap_ttest_results_1ftn<-ldply(bootstrap_ttest_results_1ftn)
bootstrap_ttest_results_1ftn$num_preceding_fortnight <- 1

bootstrap_ttest_results_2ftn<-ldply(bootstrap_ttest_results_2ftn)
bootstrap_ttest_results_2ftn$num_preceding_fortnight <- 2

bootstrap_ttest_results_3ftn<-ldply(bootstrap_ttest_results_3ftn)
bootstrap_ttest_results_3ftn$num_preceding_fortnight <- 3


# Bootstrap results
final_results<-rbind(bootstrap_ttest_results_1ftn,bootstrap_ttest_results_2ftn,bootstrap_ttest_results_3ftn)
print(final_results)
stopCluster(cl)
end_cl<-Sys.time()
print(difftime(end_cl,start_cl,units="mins"))

final_results<-final_results%>%mutate(.,
                                      ManuscriptTable_Temp = paste(round(mean_emperical_sample_mean_d.temp,3),"  (",round(d.temp_pvalue,3),")",sep=""),
                                      ManuscriptTable_AH = paste(round(mean_emperical_sample_mean_d.AH,3),"  (",round(d.AH_pvalue,3),")",sep=""))

#Bootstrap output
if(Sys.info()['sysname']=="Windows"){
  #write.csv(final_results,"C:/Users/el382/Dropbox/PhD/shaman_bootstrap/bootstrap_results.csv",row.names = FALSE)
  
}

if(Sys.info()['sysname']=="Darwin"){
  write.csv(final_results,"~/Dropbox/PhD/shaman_bootstrap/bootstrap_results.csv",row.names = FALSE)
  write.csv(ldply(bootstrap_sample_list_1ftn),"~/Dropbox/PhD/shaman_bootstrap/bootstrap_samples_1ftn.csv",row.names = FALSE)
  write.csv(ldply(bootstrap_sample_list_2ftn),"~/Dropbox/PhD/shaman_bootstrap/bootstrap_samples_2ftn.csv",row.names = FALSE)
  write.csv(ldply(bootstrap_sample_list_3ftn),"~/Dropbox/PhD/shaman_bootstrap/bootstrap_samples_3ftn.csv",row.names = FALSE)
}



#############################################################################################################
#############################################################################################################



# Comparison against historic climatic values of that particular time of year --------

# Defining functions 

#Takes an entry from epidemic table and spits out the climatic values 5 fortnights either side of onset (by default)
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
  if(point=="peak"){
    reference<-x$max_epi_fortnight
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
mean_relative_humidity_calc<-function(mean_ah,mean_temp){
  #https://carnotcycle.wordpress.com/2012/08/04/how-to-convert-relative-humidity-to-absolute-humidity/
  #http://bioma.jrc.ec.europa.eu/components/componentstools/evapotranspiration/help/Actual_vapor_pressure.html
  #http://www.reahvac.com/tools/humidity-formulas/
  denom <-6.112 * exp(17.67*mean_temp / (mean_temp+243.5)) *  2.1674
  ans<-mean_ah /denom *100
  return(ans)
}

# data manipulation
# five fortnights before and after each epidemic onset
centered_df<-adply(first_epi,1,climate_centered,.expand=FALSE,.id = NULL)

#mean climatic values for the five fortnights before and after epidemic onset
mean_stats<-centered_df%>%
  dplyr::group_by(city,relative_fortnight)%>%
  dplyr::summarise(mean_temp=mean(d.temp),
                   mean_d.temp=mean(d.temp),
                   sd_d.temp=sd(d.temp),
                   wt_d.temp=wilcox.test(d.temp)$p.value,
                   
                   mean_AH=mean(mean_AH),
                   mean_d.AH=mean(d.AH),
                   sd_d.AH=sd(d.AH),
                   wt_d.AH=wilcox.test(d.AH)$p.value)

mean_stats<-mean_stats%>%
  dplyr::rowwise()%>%
  dplyr::mutate(signif_neg.d.temp = (wt_d.temp < 0.05  & mean_d.temp < 0),
                signif_neg.d.AH = (wt_d.AH < 0.05  & mean_d.AH < 0),
                mean_RH = mean_relative_humidity_calc(mean_AH,mean_temp),
                mean_d.RH = mean_relative_humidity_calc(mean_AH + mean_d.AH,mean_temp) - mean_RH)



# plots

AT_plot<-mean_stats%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.temp))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.temp-sd_d.temp, 
                    ymax=mean_d.temp+sd_d.temp),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
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
  xlab("Fortnights Relative to Onset")+
  ylab(expression(paste("Anomalous Temperature ( ",degree,"C)")))+
  facet_grid(~as.factor(city),labeller = label_wrap_gen(width=10))+
  theme(strip.text = element_text(size=15),
        axis.title=element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

AH_plot<-mean_stats%>%
  ggplot(data=.,aes(x=relative_fortnight,y=mean_d.AH))+
  
  geom_hline(aes(yintercept = 0),size=0.4,color="black",linetype="solid") + 
  geom_vline(aes(xintercept = 0),size=0.4,color="black",linetype="solid") +
  geom_errorbar(aes(x=relative_fortnight,
                    ymin=mean_d.AH-sd_d.AH, 
                    ymax=mean_d.AH+sd_d.AH),
                size = 0.01)+
  geom_line()+
  geom_point(data=centered_df,
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
  xlab("Fortnights Relative to Onset")+
  ylab(expression(paste("Anomalous Absolute Humidity "," (g/",m^{3},")",sep="")))+
  facet_grid(~ as.factor(city),labeller = label_wrap_gen(width=10))+
  theme(strip.text = element_text(size=15),
        axis.title=element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

stacked_climate_plots<-grid.arrange(AT_plot,AH_plot,ncol=1)
ggsave(plot = stacked_climate_plots,filename = paste("C:/Users/el382/Dropbox/PhD/shaman_bootstrap/start_absolute_humidity_temp_firsts_stacked",".pdf",sep=""), 
       width=12, height=11,limitsize=FALSE)

ggsave(plot = stacked_climate_plots,filename = paste("C:/Users/el382/Dropbox/PhD/shaman_bootstrap/start_absolute_humidity_temp_firsts_stacked",".png",sep=""), 
       width=12, height=11,limitsize=FALSE)


