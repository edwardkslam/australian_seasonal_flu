library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

#Here we assume and correct for potential mis-identification during antigenic characterisation due to delays in 
#updating vaccine strain nomenclature.

# The following code will reproduce the analyses assessing the effect of antigneic change discussed in the main text:
# 1)  Comparing epidemic sizes between seasons with and without antigenic change 
#     ag_change_incidence_plot (Figure 3)
#
# 2)  Comparing epidemic onset timing between seasons with and without antigenic change
#     ag_change_start_plot (Figure S8)
#
# 3)  Comparing temporal synchrony of epidemics across cities between seasons with and without antigenic change 
#     ag_change_synchrony_plot (Figure S9)


# Loading in data ---------------------------------------------------------
epi_table<-read.csv("./dat/raw/epi_table.csv")
epi_table<-epi_table%>%
  dplyr::mutate(log_incidence = log(incidence_per_mil))

# -------------------------------------------------------------------------


epi_table<-epi_table%>%
  subset(year!=2009)%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA),
                scaled_incidence_city = log_incidence-mean(log_incidence,na.rm=TRUE))

epi_table<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::mutate(incidence_z_score_subtype_city = ifelse(epi_alarm=="Y",
                                                        (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                        NA),
                scaled_incidence_subtype_city = log_incidence - mean(log_incidence,na.rm=TRUE))

# Figure 3: Comparing epidemic sizes between seasons with and without antigenic change  -------------------------------------
# Grouped by city and subtype
ag_change_incidence_plot<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=scaled_incidence_city))+
  geom_boxplot(outlier.size=0)+ 
  
  geom_quasirandom(aes(colour=city),dodge.width=.7,cex=5,alpha=0.6)+
  
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  #scale_y_continuous(breaks=)+
  ylab("Lab confirmed incidence") +
  labs(x=NULL)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=30),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype,scales = "free_y", labeller = label_wrap_gen(width=10))


# Figure S8 Comparing epidemic onset timing between seasons with and without antigenic change ---------------------------------------------------------------
# Grouped by city and subtype

ag_change_start_plot<-epi_table%>%
  subset(epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=start))+
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


# Figure S9  Comparing temporal synchrony of epidemics across cities between seasons with and without antigenic change ---------------------------------------------------------------

# First identify the seasons in which an antigenic variant 
# causes above baseline levels of epidemic activity across all 5 cities
# remove 2009
# also remove any years, in which we are uncertain if it is the first epidemic by an antigenic variant
# ie A/H3/Moscow/10/99 in 2000, which is marked by new_ag_marker = NA

#synchrony is defined as the inverse of the variance in epidemic onset timing across the five cities for an antigenic variant in a season.

present_in_all<-epi_table%>%
  subset(.,epi_alarm=="Y" & year!=2009 & !is.na(new_ag_marker))%>%
  plyr::count(.,c("strain_year","new_ag_marker"))%>%.[.$freq==5,]%>%
  .$strain_year

synchrony_table<-epi_table%>%
  subset(.,strain_year%in%present_in_all)%>%
  dplyr::group_by(strain_year,subtype,new_ag_marker)%>%
  dplyr::summarise(synchrony=1/var(start))

ag_change_synchrony_plot<-synchrony_table%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=log(synchrony)))+
  geom_boxplot(aes(group=new_ag_marker),outlier.size=0)+ 
  geom_quasirandom(aes(colour=subtype),
                   dodge.width=.3,
                   cex=5,
                   alpha=0.6)+
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
# 
# ggsave(plot = ag_change_incidence_plot,"./figures/main/figure_3.png",
#        width=13, height=8,limitsize=FALSE)
# 
# ggsave(plot = ag_change_start_plot,"./figures/supp/figure_S8.png",
#        width=13, height=8,limitsize=FALSE)
# 
# ggsave(plot = ag_change_synchrony_plot,"./figures/supp/figure_S9.png",
#        width=12, height=5,limitsize=FALSE)




stop()


# additional plots --------------------------------------------------------
#alternative normalisation and also aggregating across subtypes

#pooled
ag_change_incidence_plot_2<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=scaled_incidence_subtype_city))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=scaled_incidence_subtype_city,colour=subtype),
              position=position_jitter(width=0.06,height=0),
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
  #scale_y_continuous(breaks=)+
  ylab("Lab confirmed incidence") +
  labs(x=NULL)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=30),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#normalise z score
ag_change_z_score_plot<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=z_score_incidence_city))+
  geom_boxplot(outlier.size=0)+ 
  geom_quasirandom(aes(colour=city),dodge.width=.7,cex=5,alpha=0.6)+
  
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  #scale_y_continuous(breaks=)+
  ylab("Lab confirmed incidence") +
  labs(x=NULL)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=30),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype,scales = "free_y", labeller = label_wrap_gen(width=10))

#z score pooled
ag_change_z_score_plot_3<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=incidence_z_score_subtype_city))+
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
  #scale_y_continuous(breaks=)+
  ylab("Lab confirmed incidence") +
  labs(x=NULL)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=30),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#z score pooled coloured by subtype
ag_change_z_score_plot_3<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=incidence_z_score_subtype_city))+
  geom_boxplot(outlier.size=0)+ 
  geom_quasirandom(aes(colour=subtype),
                   dodge.width=.7,
                   cex=5,
                   alpha=0.6)+
  
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  stat_compare_means(method = "wilcox.test", label = "p.format",label.x.npc="middle",size=5)+
  scale_x_discrete(labels=c("0"="No Ag \n Change",
                            "1"="Ag \n Change"))+
  #scale_y_continuous(breaks=)+
  ylab("Lab confirmed incidence") +
  labs(x=NULL)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title=element_text(size=30),
        axis.text.x =element_text(size=15,margin=margin(t=5,r=0,b=0,l=0)),
        axis.text.y =element_text(size=15,margin=margin(t=0,r=5,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


