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


# Here we assess compare whether the onset timings vary between subtypes or between cities

# subtype_start_plot (Figure S1)
# city_start_plot (Figure S2)

# Loading data ------------------------------------------------------------

epi_table<-read.csv("./dat/raw/epi_table.csv")


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
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=20),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

city_start_plot<-epi_table%>%subset(.,year!=2009) %>%
  ggplot(.,aes(x=city,y= start))+
  geom_boxplot(outlier.size=0)+ 
  geom_quasirandom(aes(colour=subtype),dodge.width=.7,cex=5,alpha=0.6)+
  scale_x_discrete(breaks=cities,
                   labels=substr(cities,1,3))+
  scale_y_continuous(breaks =seq(1,26,2), limits = c(0,26))+
  
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
        axis.text.x =element_text(size=20),
        axis.text.y =element_text(size=20),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# anova -------------------------------------------------------------------

lm_model<-epi_table %>%
  subset(.,year!=2009)%>%lm(start ~ subtype+city,.)

sstable <- Anova(lm_model, type = 3)
sstable<-round(sstable,3)
TukeyHSD(anova_model)


# save plot ---------------------------------------------------------------

ggsave(plot = subtype_start_plot,"./figures/supp/figure_S1.png",
       width=15, height=8,limitsize=FALSE)

ggsave(plot = city_start_plot,"./figures/supp/figure_S2.png",
       width=15, height=8,limitsize=FALSE)

write.csv(sstable,"./tables/anova_table.csv",row.names = FALSE)
