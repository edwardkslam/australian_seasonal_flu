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



epi_table%>%subset(.,epi_alarm=="Y" & first_n_biggest=="Y")%>%dplyr::group_by(year>2009)%>%dplyr::summarise(mean_dur = mean(end-start))%>%as.data.frame()

epi_table<-epi_table%>%
  dplyr::mutate(log_incidence = log(incidence_per_mil))

pre_2009<-epi_table%>%
  subset(year<2009)%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA),
                scaled_incidence_city = log_incidence-mean(log_incidence,na.rm=TRUE))



pre_2009_ag<-pre_2009%>%
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


post_2009<-epi_table%>%
  subset(year>2009)%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA),
                scaled_incidence_city = log_incidence-mean(log_incidence,na.rm=TRUE))


post_2009_ag<-post_2009%>%
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



rbind(pre_2009,post_2009)%>%
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

