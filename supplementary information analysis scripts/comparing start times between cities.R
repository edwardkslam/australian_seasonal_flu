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
# Figure S2

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table.csv")
}

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


# Figure S2 ---------------------------------------------------------------

city_comparisons<-combn(cities[order(cities)],2,simplify=FALSE)

present_in_all<-plyr::count(epi_table %>% subset(.,epi_alarm=="Y")%>%.$strain_year)%>%.[.$freq==5,]%>%.$x

paired_df<-epi_table%>%.[.$strain_year %in% present_in_all,]

paired_plot<-paired_df %>%subset(.,year!=2009) %>%
  ggplot(.,aes(x=city,y= start))+
  geom_jitter(aes(color=subtype),
              position=position_jitter(width=0.1,height=0.01),alpha=0.6,size=3.5)+
  ggpubr::stat_compare_means(method = "wilcox.test",comparisons = city_comparisons)+
  scale_x_discrete(breaks=cities,
                   labels=substr(cities,1,3))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,35))+
  
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  xlab("City")+
  ylab("Start Fortnight")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=15),
        axis.text.y =element_text(size=15),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype)

paired_plot$layers[[2]]$aes_params$textsize<-5.5


# save plot ---------------------------------------------------------------
base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/figures/supp/"
ggsave(plot = paired_plot,filename = paste(base_dir2,"figure_S2.png",sep=""), 
       width=15, height=8,limitsize=FALSE)
