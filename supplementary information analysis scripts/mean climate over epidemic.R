library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

# Here we assess whether or not climatic factors act more generally to enhance transmission, rather than as specific triggers
# for epidemic onset
# Figure S3

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  mean_fortnightly_climate_30years<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  mean_fortnightly_climate_30years<-read.csv("~/Dropbox/PhD/code for manuscript/mean_fortnightly_climate_30years.csv")
  epi_table<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table.csv")
}

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


# function to return mean climate values ----------------------------------
mean_climate_over_epi<-function(x,y="ah"){
  x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    if(y=="ah"){
      return(data.frame(mean_epi_ah=NA))
    }
    if(y=="temp"){
      return(data.frame(mean_epi_temp=NA))
    }
  }
  fortnights<-seq(x$start,x$end,1)
  temp_clim<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  temp_clim<-temp_clim%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  if(y=="ah"){
    return(data.frame(mean_epi_ah=mean(temp_clim$mean_AH)))
  }
  if(y=="temp"){
    return(data.frame(mean_epi_temp=mean(temp_clim$mean_temp)))
  }
}


# getting climate for each epidemic ---------------------------------------
epi_table_with_clim<-adply(epi_table%>%subset(.,epi_alarm=="Y"),1,mean_climate_over_epi)
epi_table_with_clim<-adply(epi_table_with_clim,1,mean_climate_over_epi,y="temp")

# Mean AH over epidemic period --------------------------------------------
mean_epi_ah_plot<-epi_table_with_clim%>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  ggplot(.,aes(x=mean_epi_ah,y= incidence_per_mil))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",size=8)+
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_y_continuous(breaks=seq(0,200,50), limits = c(0,205))+
  xlab(expression(paste("Mean Absolute Humidity over Epidemic Period ","(g/",m^{3},")",sep="")))+
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(.~city,scales = "free_x",labeller = label_wrap_gen(width=10))


# mean temp over epidemic period -------------------------------------------------

mean_epi_temp_plot<-epi_table_with_clim%>%
  subset(.,year!=2009)%>%
  ggplot(.,aes(x=mean_epi_temp,y= incidence_per_mil))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",size=8)+
  scale_colour_manual(name = "Subtype",
                      values=c("B/Yam"="#CC79A7",
                               "B/Vic"="#009E73",
                               "H1sea"="#56B4E9",
                               "H1pdm09"="#999999",
                               "H3"="#E69F00"))+
  scale_y_continuous(breaks=seq(0,200,50), limits = c(0,205))+
  xlab(expression(paste("Mean Temperature over Epidemic Period (",degree,"C)",sep="")))+
  ylab(expression(paste("Lab confirmed incidence (",10^{-6},")"))) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(.~city,scales = "free_x",labeller = label_wrap_gen(width=10))


# save plot ---------------------------------------------------------------

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(mean_epi_temp_plot)

yaxis_common<- textGrob(expression(paste("Lab confirmed incidence (",10^{-6},")")),
                        gp = gpar(fontsize = 25),
                        rot = 90, vjust = 1)

figS3<-grid.arrange(yaxis_common,
                    arrangeGrob(mean_epi_temp_plot+ theme(legend.position="none",axis.title.y = element_blank()),
                               mean_epi_ah_plot + theme(legend.position="none",axis.title.y = element_blank()),nrow=2),
                    mylegend, ncol=3,widths=c(1,14,2))

base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/figures/supp/"
ggsave(plot = figS3,filename = paste(base_dir2,"figure_S3.png",sep=""), 
       width=20, height=10,limitsize=FALSE)
