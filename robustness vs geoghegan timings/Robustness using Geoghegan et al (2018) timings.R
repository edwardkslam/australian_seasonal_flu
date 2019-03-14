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
  geog_epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/epi_table_with_geoghegan_estimates.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  geog_epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/epi_table_with_geoghegan_estimates.csv")
}


########ASSUMPTIONS########

#   1)  For each season and city which Geoghegan has data for (2007-2015), 
#       assume that I have misidentified start timing for the LARGEST Influenza Type A epidemic.
#       and replace it with their timing 
#       geog_epi_table$largest_geog_start

#   2)  For each season and city which Geoghegan has data for (2007-2015), 
#       assume that I have misidentified start timing for the EARLIEST ONSET Influenza Type A epidemic.
#       and replace it with their timing 
#       geog_epi_table$earliest_geog_start

#   3)  For each season and city which Geoghegan has data for (2007-2015), 
#       identify the epidemics within my data set that have poorly defined epidemic time series (ie absence of clear "exponential growth"
#       assume that I have misidentified start timing for these Influenza Type A epidemics.
#       and replace it with their timing 
#       geog_epi_table$poor_timeseries_geog_start


# setting directories -----------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  plot_dir<-"~/Age structure project/New folder/geohegan/"
}

# 1) Ag change and timing (largest_geog_start) ------------------------------

ag_change_start_1<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=largest_geog_start))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=largest_geog_start,colour=subtype),
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


present_in_all_1<-plyr::count(geog_epi_table%>%
                              subset(.,epi_alarm=="Y")%>%.$strain_year)%>%.[.$freq==5,]%>%.$x

synchrony_ag_change_1<-geog_epi_table%>%subset(.,epi_alarm=="Y")%>%plyr::count(.,c("strain_year","new_ag_marker"))%>%.[.$freq==5,]%>%.$strain_year
paired_df_ag_1<-geog_epi_table%>%.[.$strain_year %in% synchrony_ag_change_1,]

synchrony_measure_1<-paired_df_ag_1%>% 
  aggregate(largest_geog_start~strain_year + subtype + new_ag_marker,data=.,function(x) 1/var(x))

colnames(synchrony_measure_1)[ncol(synchrony_measure_1)]<-"synchrony"

synchrony_plot_overall_1<-synchrony_measure_1%>%
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

ggsave(plot = ag_change_start_1,filename = paste(plot_dir,"largest_geog_antigenic_change_start_4x4",".pdf",sep=""), 
       width=13, height=13,limitsize=FALSE)

ggsave(plot = synchrony_plot_overall_1,filename = paste(plot_dir,"largest_geog_ag_change_synchrony_overall",".pdf",sep=""), 
       width=12, height=5,limitsize=FALSE)

# 2) Ag change and timing (earliest_geog_start) ------------------------------
ag_change_start_2<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=earliest_geog_start))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=earliest_geog_start,colour=subtype),
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


present_in_all_2<-plyr::count(geog_epi_table%>%
                                subset(.,epi_alarm=="Y")%>%.$strain_year)%>%.[.$freq==5,]%>%.$x

synchrony_ag_change_2<-geog_epi_table%>%subset(.,epi_alarm=="Y")%>%plyr::count(.,c("strain_year","new_ag_marker"))%>%.[.$freq==5,]%>%.$strain_year
paired_df_ag_2<-geog_epi_table%>%.[.$strain_year %in% synchrony_ag_change_2,]

synchrony_measure_2<-paired_df_ag_2%>% 
  aggregate(earliest_geog_start~strain_year + subtype + new_ag_marker,data=.,function(x) 1/var(x))

colnames(synchrony_measure_2)[ncol(synchrony_measure_2)]<-"synchrony"

synchrony_plot_overall_2<-synchrony_measure_2%>%
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


ggsave(plot = ag_change_start_2,filename = paste(plot_dir,"earliest_geog_antigenic_change_start_4x4",".pdf",sep=""), 
       width=13, height=13,limitsize=FALSE)

ggsave(plot = synchrony_plot_overall_2,filename = paste(plot_dir,"earliest_geog_ag_change_synchrony_overall",".pdf",sep=""), 
       width=12, height=5,limitsize=FALSE)

# 3) Ag change and timing (poor_timeseries_geog_start) ------------------------------
ag_change_start_3<-geog_epi_table%>%
  subset(.,!(year==2009 & subtype!="H1pdm09"))%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=poor_timeseries_geog_start))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=poor_timeseries_geog_start,colour=subtype),
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


present_in_all_3<-plyr::count(geog_epi_table%>%
                                subset(.,epi_alarm=="Y")%>%.$strain_year)%>%.[.$freq==5,]%>%.$x

synchrony_ag_change_3<-geog_epi_table%>%subset(.,epi_alarm=="Y")%>%plyr::count(.,c("strain_year","new_ag_marker"))%>%.[.$freq==5,]%>%.$strain_year
paired_df_ag_3<-geog_epi_table%>%.[.$strain_year %in% synchrony_ag_change_3,]

synchrony_measure_3<-paired_df_ag_3%>% 
  aggregate(poor_timeseries_geog_start~strain_year + subtype + new_ag_marker,data=.,function(x) 1/var(x))

colnames(synchrony_measure_3)[ncol(synchrony_measure_3)]<-"synchrony"

synchrony_plot_overall_3<-synchrony_measure_3%>%
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


ggsave(plot = ag_change_start_3,filename = paste(plot_dir,"poorly_defined_geog_antigenic_change_start_4x4",".pdf",sep=""), 
       width=13, height=13,limitsize=FALSE)

ggsave(plot = synchrony_plot_overall_3,filename = paste(plot_dir,"poorly_defined_geog_ag_change_synchrony_overall",".pdf",sep=""), 
       width=12, height=5,limitsize=FALSE)