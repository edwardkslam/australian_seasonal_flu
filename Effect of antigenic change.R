library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
# library(ggmap)
library(tidyr)

epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")


# Figure 3: Comparing epidemic sizes between seasons with and without antigenic change  -------------------------------------
# Grouped by city and subtype
ag_change_incidence_plot<-epi_table%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=incidence_per_mil))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=incidence_per_mil,colour=subtype),
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
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  labs(x=NULL)+
  theme_bw()+
  facet_grid(subtype~city,scales = "free_y", labeller = label_wrap_gen(width=10))+ 
  theme(strip.text = element_text(size=15),
        axis.title=element_text(size=13),
        axis.text = element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Figure S8 Comparing epidemic onset timing between seasons with and without antigenic change ---------------------------------------------------------------
# Grouped by city and subtype

ag_change_start_plot<-epi_table%>%
  subset(epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,!is.na(new_ag_marker))%>%
  ggplot(.,aes(x=as.factor(new_ag_marker),y=start))+
  geom_boxplot(outlier.size=0)+ 
  geom_jitter(aes(x=as.factor(new_ag_marker),y=start,colour=subtype),
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

synchrony_table%>%
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
  xlab(NULL)+
  ylab("ln(Synchrony) (1/Var)")+
  theme_bw()+
  theme(axis.title=element_text(size=13),
        axis.text = element_text(size=13),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


