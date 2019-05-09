library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess compare whether the onset timings between onset timings of epidemics by subytpe between cities
# Figure S1

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table.csv")
}

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


# plot --------------------------------------------------------------------

subtype_list<-c("H3","H1sea","H1pdm09","B/Vic","B/Yam")
subtype_comparisons<-combn(subtype_list,2,simplify=FALSE)

subtype_start_plot<-epi_table %>%
  subset(.,year!=2009)%>%
  ggplot(.,aes(x=subtype,y= start))+
  geom_jitter(aes(color=subtype),
              position=position_jitter(width=0.1,height=0.01),alpha=0.6,size=3.5)+
  ggpubr::stat_compare_means(method = "wilcox.test",comparisons = subtype_comparisons,textsize =10)+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,38))+
  xlab("Subtype")+
  ylab("Start Fortnight")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        panel.border = element_rect(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=25),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~city)


# save plot ---------------------------------------------------------------
base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/figures/supp/"
ggsave(plot = subtype_start_plot,filename = paste(base_dir2,"figure_S1.png",sep=""), 
       width=15, height=8,limitsize=FALSE)

