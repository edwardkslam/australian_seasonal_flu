library(magrittr)
library(ggplot2)
library(ggpubr)
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
#       i) ag_change_start_1 (Figure S10)
#       ii) synchrony_plot_overall_1 (Figure S11)

#   2)  For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the EARLIEST influenza A subtype 
#       was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       The timing value used for this set of analyses is stored in geog_epi_table$earliest_geog_start
#       Outputs
#       i) ag_change_start_1 (Figure S12)
#       ii) synchrony_plot_overall_1 (Figure S13)

#   3)  For seasons between 2007 and 2015 in which the number of cases for the dominant influenza A subtype were small or 
#       it was difficult to discern the period of epidemic from background activity, 
#       we assumed that our timing estimate was incorrect and replaced it with estimates from Geoghegan et al. (2018).
#       Outputs
#       i) ag_change_start_1 (Figure S14)
#       ii) synchrony_plot_overall_1 (Figure S15)


# loading in data ---------------------------------------------------------

if(Sys.info()['sysname']=="Windows"){
  geog_epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/epi_table_with_geoghegan_estimates.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  geog_epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/robustness vs geoghegan timings/epi_table_with_geoghegan_estimates.csv")
}


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
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   legend.position="none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(subtype~city,scales = "free_y", labeller = label_wrap_gen(width=10))


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

ggsave(plot = ag_change_start_1,filename = paste(plot_dir,"largest_geog_antigenic_change_start_4x4",".pdf",sep=""),
       width=13, height=13,limitsize=FALSE)

ggsave(plot = synchrony_plot_overall_1,filename = paste(plot_dir,"largest_geog_ag_change_synchrony_overall",".pdf",sep=""),
       width=12, height=5,limitsize=FALSE)

# 2) Ag change and timing (earliest_geog_start) ------------------------------
ag_change_start_2<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09"& year!=2009)%>%
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
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   legend.position="none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(subtype~city,scales = "free_y", labeller = label_wrap_gen(width=10))


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


ggsave(plot = ag_change_start_2,filename = paste(plot_dir,"earliest_geog_antigenic_change_start_4x4",".pdf",sep=""),
       width=13, height=13,limitsize=FALSE)

ggsave(plot = synchrony_plot_overall_2,filename = paste(plot_dir,"earliest_geog_ag_change_synchrony_overall",".pdf",sep=""),
       width=12, height=5,limitsize=FALSE)

# 3) Ag change and timing (poor_timeseries_geog_start) ------------------------------
ag_change_start_3<-geog_epi_table%>%
  subset(.,epi_alarm=="Y" & subtype!="H1pdm09" & year!=2009)%>%
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
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(size=20),
                   axis.title=element_text(size=20),
                   axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
                   axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
                   axis.ticks.length = unit(0.4,"cm"),
                   panel.border = element_rect(colour = "black"),
                   legend.position="none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())+
  facet_grid(subtype~city,scales = "free_y", labeller = label_wrap_gen(width=10))


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



# save plots --------------------------------------------------------------

base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/figures/supp/"

ggsave(plot = ag_change_start_1,filename = paste(base_dir2,"figure_S10.png",sep=""), 
       width=13, height=13,limitsize=FALSE)
ggsave(plot = synchrony_plot_overall_1,filename = paste(base_dir2,"figure_S11.png",sep=""), 
       width=12, height=5,limitsize=FALSE)

ggsave(plot = ag_change_start_2,filename = paste(base_dir2,"figure_S12.png",sep=""), 
       width=13, height=13,limitsize=FALSE)
ggsave(plot = synchrony_plot_overall_2,filename = paste(base_dir2,"figure_S13.png",sep=""), 
       width=12, height=5,limitsize=FALSE)

ggsave(plot = ag_change_start_3,filename = paste(base_dir2,"figure_S14.png",sep=""), 
       width=13, height=13,limitsize=FALSE)
ggsave(plot = synchrony_plot_overall_3,filename = paste(base_dir2,"figure_S15.png",sep=""), 
       width=12, height=5,limitsize=FALSE)