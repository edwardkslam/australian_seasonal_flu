library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess the effect of prior epidemic activity (within the same season and city) by antigenic variants of other subtypes on
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


# load data ---------------------------------------------------------------
epi_table<-read.csv("./dat/raw/epi_table.csv")
epi_table<-epi_table%>%
  dplyr::mutate(log_incidence = log(incidence_per_mil))

# 1)  Relationship between Prior Activity and Relative Size of an Epidemic --------
log_prior_activity_vs_relative_size_plot<-epi_table %>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,year!=2009)%>%
  subset(.,delay!=0)%>%
  ggplot(.,aes(x=log(prior_everything_scaled),y= log(relative_to_first_n_biggest)))+
  geom_hline(aes(yintercept = 0),size=0.3,color="black",linetype="solid") +
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=5)+
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
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=20),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        plot.margin=margin(1.5,0.5,0,0,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# 2)  Delay and Relative Size of an Epidemic ----------------------------------

log_delayed_delay_vs_relative_size_plot<-epi_table %>%
  subset(.,year!=2009)%>%
  subset(.,delay!=0)%>%
  ggplot(.,aes(x=delay,y= log(relative_to_first_n_biggest)))+
  geom_hline(aes(yintercept = 0),size=0.3,color="black",linetype="solid") +
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.05,height=0.05),alpha=0.6,size=5)+
  stat_cor(method = "pearson",size=8,label.x.npc = "centre")+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_x_continuous(breaks=seq(0,9,1), limits = c(-0.1,9.1))+
  xlab("Delay (two week increments)")+
  ylab("ln(Size relative to First of Season)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=20),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        plot.margin=margin(1.5,0,0,0.5,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# save plot ---------------------------------------------------------------
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(log_prior_activity_vs_relative_size_plot)

fig5<-ggarrange(log_prior_activity_vs_relative_size_plot+ theme(legend.position="none"),
                log_delayed_delay_vs_relative_size_plot + theme(legend.position="none"),
                mylegend,
                labels = c("a","b",NA),
                font.label = list(size=22),
                ncol = 3, nrow = 1,widths = c(5,5,1))

ggsave(plot = fig5,filename = "./figures/main/figure_5.png", 
       width=20, height=8,limitsize=FALSE)

ggsave(plot = fig5,filename = "./figures/main/figure_5.eps", 
       width=20, height=8,limitsize=FALSE,device = cairo_ps)
