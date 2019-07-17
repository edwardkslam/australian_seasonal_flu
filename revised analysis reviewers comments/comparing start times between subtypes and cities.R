library(car)
library(knitr)
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

epi_table<-epi_table%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(normalised_start = (start-mean(start,na.rm=TRUE))/sd(start,na.rm=TRUE))

# plot --------------------------------------------------------------------


subtype_start_plot<-epi_table %>%
  subset(.,year!=2009)%>%
  ggplot(.,aes(x=subtype,y= start))+
  geom_boxplot(outlier.size=0)+ 
  geom_quasirandom(aes(colour=city),dodge.width=.7,cex=5,alpha=0.6)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  xlab("Subtype")+
  ylab("Start Fortnight")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        panel.border = element_rect(colour = "black"),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x =element_text(size=20),
        axis.title.y=element_text(size=25),
        axis.text.y =element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=17),
        #legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# anova -------------------------------------------------------------------

lm_model<-epi_table %>%
  subset(.,year!=2009)%>%lm(start ~ subtype+city,.)

sstable <- Anova(lm_model, type = 3)
sstable<-round(sstable,3)
TukeyHSD(anova_model)


# save plot ---------------------------------------------------------------
base_dir2<-"C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/figures/reviewer comments/"
ggsave(plot = subtype_start_plot,filename = paste(base_dir2,"figure_S1.png",sep=""), 
       width=15, height=8,limitsize=FALSE)

write.csv(sstable,paste(base_dir2,"anova_table.csv",sep=""))
