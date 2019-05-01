library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess the effect of prior epidemic activity (within the same season andcity) by antigenic variants of other subtypes on
# the size of an epidemic.

# The following code will reproduce the prior epidemic activity analyses discussed in the main text:
# 1)  The relationship between prior epidemic activity by other subtypes and relative size of an epidemic.
#     Prior Activity is measured relative to the city-specific mean epidemic size.
#     The size of a specific epidemic is measured relative to that of the earliest onset epidemic of the season.
#     log_prior_activity_vs_relative_size_plot (Figure 5)
#
# 2)  The relationship between delay in onset timing and relative size of an epidemic.
#     Delay is measured as the difference in onset timing between a specific epidemic and the earliest onset epidemic 
#     within a particular season and city
#     The size of a specific epidemic is measured relative to that of the earliest onset epidemic of the season.
#     log_delayed_delay_vs_relative_size_plot (Figure 5)


# load in summary data ----------------------------------------------------
epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")


# 1)  Relationship between Prior Activity and Relative Size of an Epidemic --------
log_prior_activity_vs_relative_size_plot<-epi_table %>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,delay!=0)%>%
  ggplot(.,aes(x=log(prior_everything_scaled + 0.01),y= log(relative_to_first_n_biggest  +0.01)))+
  geom_hline(aes(yintercept = 0),size=0.3,color="black",linetype="solid") +
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",size=8,label.x.npc = "centre")+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  xlab("ln(Prior Activity of other Antigenic Variants)")+
  ylab("ln(Size relative to First of Season)")+
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# 2)  Delay and Relative Size of an Epidemic ----------------------------------

log_delayed_delay_vs_relative_size_plot<-epi_table %>%
  #subset(.,city!="TOWNSVILLE")%>%
  subset(.,year!=2009)%>%
  subset(.,delay!=0)%>%
  ggplot(.,aes(x=delay,y= log(relative_to_first_n_biggest  +0.01)))+
  geom_hline(aes(yintercept = 0),size=0.3,color="black",linetype="solid") +
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.05,height=0.05),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",size=8,label.x.npc = "centre")+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_x_continuous(breaks=seq(0,9,1), limits = c(-0.1,9.1))+
  xlab("Delay")+
  ylab("ln(Size relative to First of Season)")+
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
