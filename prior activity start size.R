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

epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/shaman_bootstrap/epi_table.csv")
# plots -------------------------------------------------------------------

size_just_earliest<-epi_table %>%
  subset(.,year!=2009)%>%
  subset(.,first_n_biggest=="Y")%>%
  #subset(.,first_n_biggest=="Y")%>%
  ggplot(.,aes(x=start,y= incidence_per_mil))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.1,height=0.01),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",label.x = 13,size=8)+
  #scale_y_continuous(breaks =seq(0,1,0.25), limits = c(0,1))+
  #scale_x_continuous(breaks = seq(0,200,100), limits = c(0,200))+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_x_continuous(breaks=seq(1,26,2), limits = c(-0.9,26.1))+
  scale_y_continuous(breaks=seq(0,250,50), limits = c(0,250))+
  ggtitle("Size vs start time for the first and largest of a season")+
  xlab("Start")+
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.text = element_text(size=25),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~city, labeller = label_wrap_gen(width=10))

ggsave(plot = size_just_earliest,filename = "C:/Users/el382/Desktop/size_of_firsts.pdf", 
       width=18, height=6,units = "in",limitsize=FALSE)


size_all<-epi_table %>%
  subset(.,year!=2009)%>%
  #subset(.,first_n_biggest=="Y")%>%
  ggplot(.,aes(x=start,y= incidence_per_mil))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.1,height=0.01),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",label.x = 13,size=8)+
  #scale_y_continuous(breaks =seq(0,1,0.25), limits = c(0,1))+
  #scale_x_continuous(breaks = seq(0,200,100), limits = c(0,200))+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_x_continuous(breaks=seq(1,26,2), limits = c(-0.9,26.1))+
  scale_y_continuous(breaks=seq(0,250,50), limits = c(0,250))+
  ggtitle("Size vs start time for ALL")+
  xlab("Start")+
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.text = element_text(size=25),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~city, labeller = label_wrap_gen(width=10))

ggsave(plot = size_all,filename = "C:/Users/el382/Desktop/size_of_all.pdf", 
       width=18, height=6,units = "in",limitsize=FALSE)


size_delayed<-epi_table %>%
  subset(.,year!=2009)%>%
  subset(.,delay!=0)%>%
  ggplot(.,aes(x=start,y= incidence_per_mil))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.1,height=0.01),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",label.x = 13,size=8)+
  #scale_y_continuous(breaks =seq(0,1,0.25), limits = c(0,1))+
  #scale_x_continuous(breaks = seq(0,200,100), limits = c(0,200))+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_x_continuous(breaks=seq(1,26,2), limits = c(-0.9,26.1))+
  scale_y_continuous(breaks=seq(0,250,50), limits = c(0,250))+
  ggtitle("Size vs start time for delayed")+
  xlab("Start")+
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.text = element_text(size=25),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~city, labeller = label_wrap_gen(width=10))

ggsave(plot = size_delayed,filename = "C:/Users/el382/Desktop/size_delayed.pdf", 
       width=18, height=6,units = "in",limitsize=FALSE)


# multilevel models -------------------------------------------------------


m0<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lm(formula= (incidence_per_mil) ~ prior_everything_scaled, data=.)

m1<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lmer(formula= (incidence_per_mil) ~ prior_everything_scaled + (1|city), data=.)

m2<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lmer(formula= (incidence_per_mil) ~ prior_everything_scaled + start + (1|city), data=.)

anova(m2,m1,m0)


m0.1<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lm(formula= (incidence_per_mil) ~ start, data=.)

m1.1<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lmer(formula= (incidence_per_mil) ~ start + (1|city), data=.)

m2.1<-epi_table %>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  lmer(formula= (incidence_per_mil) ~ prior_everything_scaled + start + (1|city), data=.)

anova(m2.1,m1.1,m0.1)
