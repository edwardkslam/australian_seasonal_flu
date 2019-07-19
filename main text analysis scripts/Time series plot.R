library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)


# This produces time_series_plot (Figure 1).
# Note, submitted samples from 2009 were omitted since the pandemic overwhelmed laboratory testing
# so we were not able to accurately estimate epidemic size, timing etc.


# load in data ------------------------------------------------------------


raw_table<-read.csv("./dat/raw/raw_data.csv")

time_series_plot<-raw_table%>%
  subset(.,year!=2009)%>%
  dplyr::mutate(fortnight_date = cut.Date(as.Date(specimen_date),breaks = "2 week")%>%as.Date(.))%>%
  plyr::count(.,c("city","subtype","fortnight_date"))%>%
  ggplot(., aes(x = fortnight_date,y = freq, group = subtype))+
  geom_bar(aes(fill = subtype),alpha=0.7, position = "identity" ,stat="identity") + 
  annotate("rect", xmin=as.Date("2009-01-01"), xmax=as.Date("2009-12-31"),
           ymin=c(0) , ymax=c(140), alpha=0.2, fill="red")+
  annotate("text", x = as.Date("2009-07-01"), y = 70, label = c("Pandemic"),
           color = "black", size=7 , angle=90, fontface="bold")+
  scale_fill_manual(name = "Subtype",
                    values=c("B/Yam"="#CC79A7",
                             "B/Vic"="#009E73",
                             "H1sea"="#56B4E9",
                             "H1pdm09"="#999999",
                             "H3"="#E69F00"))+
  xlab("Year")+
  ylab("No. Lab Confirmed Cases")+
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("year"),
               limits =  as.Date(c('2000-01-01','2016-01-01')),expand=c(0,0))+
  scale_y_continuous(breaks =seq(0,140,20), limits = c(0,140))+
  theme_bw()+
  theme(axis.title=element_text(size=25),
        strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=15),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17), 
        panel.grid.major.x = element_line(colour = "grey",size = 0.01),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(city~.,labeller = label_wrap_gen(width=10))
  


# saving plots ------------------------------------------------------------

ggsave(plot = time_series_plot,"./figures/main/figure_1.png",
       width=18, height=10,limitsize=FALSE)
