library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess the robustness of our analyses of the effects of antigenic change on i) epidemic onset timing; ii) temporal synchrony of epidemics 
# across all five cities; to potential inaccuracies in the estimation of epidemic onset timing.
# We incorporate alternative estimates given by Geoghegan et al. (2018), based on a series of assumptions 1 to 3.


########ASSUMPTIONS########

#   1)  For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the DOMINANT influenza A subtype 
#       was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       The timing value used for this set of analyses is stored in geog_epi_table$largest_geog_start
#       Outputs
#       i) ag_change_start_1
#       ii) synchrony_plot_overall_1 

#   2)  For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the EARLIEST influenza A subtype 
#       was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       The timing value used for this set of analyses is stored in geog_epi_table$earliest_geog_start
#       Outputs
#       i) ag_change_start_1
#       ii) synchrony_plot_overall_1

#   3)  For seasons between 2007 and 2015 in which the number of cases for the dominant influenza A subtype were small or 
#       it was difficult to discern the period of epidemic from background activity, 
#       we assumed that our timing estimate was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       Outputs
#       i) ag_change_start_1
#       ii) synchrony_plot_overall_1 


# loading in data ---------------------------------------------------------

geog_epi_table<-read.csv("./dat/raw/geoghegan/epi_table_with_geoghegan_estimates.csv")

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

geog_epi_table$city<-factor(geog_epi_table$city,levels = cities)

# setting directories -----------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  plot_dir<-"~/Age structure project/New folder/geohegan/"
}

# 1) i) Ag change and timing (largest_geog_start) ------------------------------
ag_change_start_1<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09" & year!=2009)%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=largest_geog_start))+
  geom_boxplot(outlier.size=0)+
  geom_quasirandom(aes(colour=city),
                   dodge.width=.7,
                   cex=5,
                   alpha=0.6)+
  
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26.2))+
  ylab("Start Fortnight") +
  labs(x=NULL)+
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(~subtype,scales = "free_y", labeller = label_wrap_gen(width=10))


# 1) ii) Ag change and temporal synchrony (largest_geog_start) -----------------------------------------------------------------------
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
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=8)+
  xlab(NULL)+
  ylab("ln(Synchrony) (1/Var)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# 2) Ag change and timing (earliest_geog_start) ------------------------------
ag_change_start_2<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09"& year!=2009)%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=earliest_geog_start))+
  geom_boxplot(outlier.size=0)+
  geom_quasirandom(aes(colour=city),
                   dodge.width=.7,
                   cex=5,
                   alpha=0.6)+
  
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26.2))+
  ylab("Start Fortnight") +
  labs(x=NULL)+
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(~subtype,scales = "free_y", labeller = label_wrap_gen(width=10))


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
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=8)+
  xlab(NULL)+
  ylab("ln(Synchrony) (1/Var)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# 3) Ag change and timing (poor_timeseries_geog_start) ------------------------------
ag_change_start_3<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09" & year!=2009)%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=poor_timeseries_geog_start))+
  geom_boxplot(outlier.size=0)+
  geom_quasirandom(aes(colour=city),
                   dodge.width=.7,
                   cex=5,
                   alpha=0.6)+
  
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26.2))+
  ylab("Start Fortnight") +
  labs(x=NULL)+
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(~subtype,scales = "free_y", labeller = label_wrap_gen(width=10))


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
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=8)+
  xlab(NULL)+
  ylab("ln(Synchrony) (1/Var)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


