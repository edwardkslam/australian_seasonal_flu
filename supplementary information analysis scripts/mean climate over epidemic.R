library(lubridate)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)
library(zoo)

# Here we assess whether or not climatic factors act more generally to enhance transmission, rather than as specific triggers
# for epidemic onset
# Figure S3

# Loading data ------------------------------------------------------------

mean_fortnightly_climate_30years<-read.csv("./dat/raw/mean_fortnightly_climate_30years.csv")
epi_table<-read.csv("./dat/raw/epi_table.csv")
  
raw_table<-read.csv("./dat/raw/raw_data.csv")
raw_table<-raw_table%>%
  dplyr::mutate(fortnights_since_start_of_year = lubridate::yday(specimen_date)%/%14+1)%>%
  dplyr::group_by(city,year,assumed_antigenic_variant,fortnights_since_start_of_year)%>%
  dplyr::summarise(count=n())

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


# function to return mean climate values ----------------------------------
mean_climate_over_epi<-function(x){
  #x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    return(data.frame(x,
                      mean_epi_ah=NA,
                      mean_epi_temp=NA))
  }
  fortnights<-seq(x$start,x$end,1)
  temp_clim1<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  temp_clim<-temp_clim1%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  
  #roll_mean_AH<-temp_clim1$mean_AH%>%rollmean(.,k=length(fortnights))
  #roll_mean_temp<-temp_clim1$mean_temp%>%rollmean(.,k=length(fortnights))
  
  temp_clim<-temp_clim%>%dplyr::summarise(mean_epi_ah = mean(mean_AH),
                                          mean_epi_temp = mean(mean_temp)#,
                                          #z_score_mean_epi_ah = (mean_epi_ah-mean(roll_mean_AH,na.rm=TRUE))/sd(roll_mean_AH,na.rm=TRUE),
                                          #z_score_mean_epi_temp = (mean_epi_temp-mean(roll_mean_temp,na.rm=TRUE))/sd(roll_mean_temp,na.rm=TRUE)
                                          )
  return(data.frame(x,temp_clim))
}

early_climate<-function(x){
  #mean climate over start to peak fortnight
  x<-as.data.frame(x)
  if(x$epi_alarm=="N"){
    return(data.frame(x,
                      early_ah=NA,
                      early_temp=NA))
  }
  
  temp<-raw_table%>%subset(.,city==x$city & year==x$year & assumed_antigenic_variant==as.character(x$reference_strain) &
                             fortnights_since_start_of_year %in%c(x$start+1:x$end))
  peak_fortnight<-temp$fortnights_since_start_of_year[which.max(temp$count)]
  
  fortnights<-seq(x$start,peak_fortnight,1)
  temp_clim1<-subset(mean_fortnightly_climate_30years,city ==as.character(x$city))
  
  #roll_mean_AH<-temp_clim1$mean_AH%>%rollmean(.,k=length(fortnights))
  #roll_mean_temp<-temp_clim1$mean_temp%>%rollmean(.,k=length(fortnights))
  
  temp_clim<-temp_clim1%>%subset(.,year==x$year & fortnights_since_start_of_year%in% fortnights)
  temp_clim<-temp_clim%>%dplyr::summarise(early_ah = mean(mean_AH),
                                          early_temp = mean(mean_temp)#,
                                          #z_score_early_ah = (early_ah-mean(roll_mean_AH,na.rm=TRUE))/sd(roll_mean_AH,na.rm=TRUE),
                                          #z_score_early_temp = (early_temp-mean(roll_mean_temp,na.rm=TRUE))/sd(roll_mean_temp,na.rm=TRUE)
                                          )
  
}




# getting climate for each epidemic ---------------------------------------
epi_table_with_clim<-adply(epi_table%>%subset(.,epi_alarm=="Y" & year!=2009),1,mean_climate_over_epi)

epi_table_with_clim<-adply(epi_table_with_clim,1,early_climate)


# calculating z score for epidemic size and climate  ----------------------
#in order to make it comparable between cities
epi_table_with_clim<-epi_table_with_clim%>%
  dplyr::mutate(scaled_incidence_city = incidence_per_mil/mean_epi_size,
                log_incidence = log(incidence_per_mil))

mean_size_subtype_city<-epi_table_with_clim%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::summarise(mean_epi_size_sc = mean(log_incidence,na.rm=TRUE))

epi_table_with_clim<-left_join(epi_table_with_clim,mean_size_subtype_city)
epi_table_with_clim<-epi_table_with_clim%>%
  dplyr::mutate(scaled_incidence_subtype_city = log_incidence-mean_epi_size_sc)

epi_table_with_clim<-epi_table_with_clim%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::mutate(incidence_z_score_subtype_city = ifelse(epi_alarm=="Y",
                                                        (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                        NA))
epi_table_with_clim<-epi_table_with_clim%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_mean_epi_ah = ifelse(epi_alarm=="Y",
                                         (mean_epi_ah-mean(mean_epi_ah,na.rm=TRUE))/sd(mean_epi_ah,na.rm=TRUE)))

# Mean AH over epidemic period --------------------------------------------
mean_epi_ah_plot<-epi_table_with_clim%>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  ggplot(.,aes(x=mean_epi_ah,y= scaled_incidence_subtype_city))+
  geom_jitter(aes(group = city, color=city),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  
  scale_x_continuous(breaks=seq(6,16,2),limits = c(5,16))+
  scale_y_continuous(breaks=seq(-4,2,1),limits = c(-4,2))+
  xlab(expression(paste("Mean Absolute Humidity over Epidemic Period "," (g/",m^{3},")",sep="")))+
  ylab("Lab confirmed incidence") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=18),
        axis.text.x =element_text(size=16,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=16,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# mean temp over epidemic period -------------------------------------------------

mean_epi_temp_plot<-epi_table_with_clim%>%
  subset(.,year!=2009  & epi_alarm=="Y")%>%
  ggplot(.,aes(x=mean_epi_temp,y= scaled_incidence_subtype_city))+
  geom_jitter(aes(group = subtype, color=city),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  scale_x_continuous(breaks=seq(8,24,4),limits = c(8,24))+
  scale_y_continuous(breaks=seq(-4,2,1),limits = c(-4,2))+
  xlab(expression(paste("Mean Temperature over Epidemic Period (",degree,"C)",sep="")))+
  ylab("Lab confirmed incidence") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=18),
        axis.text.x =element_text(size=16,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=16,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# mean AH over early epidemic ---------------------------------------------
early_epi_ah_plot<-epi_table_with_clim%>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  ggplot(.,aes(x=early_ah,y= scaled_incidence_subtype_city))+
  geom_jitter(aes(group = city, color=city),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  
  scale_x_continuous(breaks=seq(6,16,2),limits = c(5,16))+
  scale_y_continuous(breaks=seq(-4,2,1),limits = c(-4,2))+
  xlab(expression(paste("Mean Absolute Humidity over Early Epidemic "," (g/",m^{3},")",sep="")))+
  ylab("Lab confirmed incidence") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=18),
        axis.text.x =element_text(size=16,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=16,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# early temp over epidemic period -------------------------------------------------

early_epi_temp_plot<-epi_table_with_clim%>%
  subset(.,year!=2009 & epi_alarm=="Y")%>%
  ggplot(.,aes(x=early_temp,y= scaled_incidence_subtype_city))+
  geom_jitter(aes(group = subtype, color=city),
              position=position_jitter(width=0.1,height=0.05),alpha=0.6,size=3.5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  
  scale_x_continuous(breaks=seq(8,24,4),limits = c(8,24))+
  scale_y_continuous(breaks=seq(-4,2,1),limits = c(-4,2))+
  xlab(expression(paste("Mean Temperature over Early Epidemic (",degree,"C)",sep="")))+
  ylab("Lab confirmed incidence") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=18),
        axis.text.x =element_text(size=16,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=16,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
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

mylegend<-g_legend(mean_epi_temp_plot)

yaxis_common<- textGrob(expression(paste("Lab confirmed incidence (",10^{-6},")")),
                        gp = gpar(fontsize = 25),
                        rot = 90, vjust = 1)

figS3<-grid.arrange(yaxis_common,
                    arrangeGrob(mean_epi_temp_plot+ theme(legend.position="none",axis.title.y = element_blank()),
                               mean_epi_ah_plot + theme(legend.position="none",axis.title.y = element_blank()),
                               early_epi_temp_plot + theme(legend.position="none",axis.title.y = element_blank()),
                               early_epi_ah_plot + theme(legend.position="none",axis.title.y = element_blank()),nrow=2),
                    mylegend, ncol=3,widths=c(1,14,2))

ggsave(plot = figS3,"./figures/supp/figure_S3.png",
       width=18, height=12,limitsize=FALSE)
