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

ag_change_size_1<-largest_use_geog%>%
  subset(.,!(year==2009 & subtype!="H1pdm09"))%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=new_start))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=new_start,colour=subtype),
              position=position_jitter(width=0.1,height=0.1),
              alpha=0.6,
              size=5)+
  
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26.2))+
  ylab("Start Fortnight") +
  labs(x=NULL)+
  theme_bw()+
  facet_grid(subtype~city, labeller = label_wrap_gen(width=10))+ 
  theme(strip.text = element_text(size=15),
        axis.title=element_text(size=13),
        axis.text = element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


present_in_all_1<-plyr::count(largest_use_geog%>%
                              subset(.,epi_alarm=="Y")%>%.$strain_year)%>%.[.$freq==5,]%>%.$x

synchrony_ag_change_1<-largest_use_geog%>%subset(.,epi_alarm=="Y")%>%plyr::count(.,c("strain_year","new_ag_marker"))%>%.[.$freq==5,]%>%.$strain_year
paired_df_ag_1<-largest_use_geog%>%.[.$strain_year %in% synchrony_ag_change_1,]

synchrony_measure_1<-paired_df_ag_1%>% 
  aggregate(start~strain_year + subtype + new_ag_marker,data=.,function(x) 1/var(x))

colnames(synchrony_measure_1)[ncol(synchrony_measure_1)]<-"synchrony"

synchrony_plot_overall<-synchrony_measure_1%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=log(synchrony)))+
  geom_boxplot(aes(group=new_ag_marker),outlier.size=0)+ 
  geom_jitter(aes(color=subtype),position=position_jitter(width=0.05,height=0.005),alpha=0.6,size=5)+
  #scale_y_continuous(breaks =seq(0,5,1), limits = c(-0.1,5.1))+
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=8)+
  #stat_summary(fun.y = mean, geom = "point",colour = "black", size=2) +
  #stat_summary(fun.data = meanFunction, geom ="text", color = "black", size = 3, vjust = 1.3,hjust=-0.5)+
  #ggtitle("Synchrony for ag variants that cause epidemics in all 4 cities")+
  xlab(NULL)+
  ylab("ln(Synchrony) (1/Var)")+
  theme_bw()+
  theme(axis.title=element_text(size=13),
        axis.text = element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())