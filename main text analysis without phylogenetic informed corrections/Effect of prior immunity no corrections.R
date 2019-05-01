library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

#Here we DO NOT correct for potential mis-identification during antigenic characterisation due to delays in 
#updating vaccine strain nomenclature.

# The following code will reproduce the prior immunity analyses discussed in the main text:
# 1)  The relationship between epidemic incidence and the amount of antigenic variant-specific cumulative incidence 
#     epi_size_cumulative_size_same_variant_plot (Figure S20)
#
# 2)  The relationship between the probability of successful epidemic initiation 
#     and amount of antigenic variant-specific cumulative incidence for each subtype 
#     prob_successful_epi_cumulative_size_same_variant_plot (Figure S21)
#
# 3)  Binary logistics regression assessing the effect of antigenic variant-specific cumulative incidence 
#     on the probability of successful epidemic initiation for each subtype 
#     subtype_logistics_regression (Table S9)

# Loading data ------------------------------------------------------------
if(Sys.info()['sysname']=="Windows"){
  epi_table_no_corrections<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/epi_table_no_corrections.csv")
}

if(Sys.info()['sysname']=="Darwin"){
  epi_table_no_corrections<-read.csv("~/Dropbox/PhD/code for manuscript/epi_table_no_corrections.csv")
}

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table_no_corrections$city<-factor(epi_table_no_corrections$city,levels = cities)


# 1) Effect of ag variant specific cumulative incidence on epidemic incidence --------

# new_ag_marker = NA means that the variant emerged prior to start of study period, 
# ie we are unable to calculate the cumulative incidence
# also, due to the pandemic H1N1 emergence in 2009, epidemic activity could not be reliable quantified for A/H3/Perth/16/2009
# these need to be removed
dodgy_prev<-epi_table_no_corrections%>%subset(.,is.na(new_ag_marker) | 
                                 reference_strain=="A/Perth/16/2009-like")

# First find the first and last year in which an ag variant causes an epidemic in each city
first_emerge_table<-epi_table_no_corrections%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::summarise(emerge_year=min(year),
                   last_recorded_epi=max(year))

# assume that an ag variant could potentially have caused an epidemic right up to the year 
# before the emergence of its replacement new variant
first_emerge_table<-first_emerge_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::arrange(emerge_year,.by_group=TRUE)%>%
  dplyr::mutate(last_possible_year=if(subtype!="H1sea"){lead(emerge_year-1,default = 2015)}else{lead(emerge_year-1,default = 2008)})

# the last possible year in which an epidemic could have been caused is whichever occured last:
# the last recorded epidemic (to account for the case in which old and new variants circulate in same year)
# the year prior to emergence of new variant
first_emerge_table<-first_emerge_table%>%
  dplyr::rowwise()%>%
  dplyr::mutate(last_possible_year=max(last_possible_year,last_recorded_epi))

# correcting for mislabelling 
end_rows<-which(first_emerge_table$subtype=="H1sea" & first_emerge_table$last_possible_year>=2009)
first_emerge_table$last_possible_year[end_rows]<-2008

# generate the full list of possible years for each variant
cumulative_incidence_by_ag<-first_emerge_table%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  tidyr::expand(.,year=full_seq(c(emerge_year,last_possible_year),1))

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,reference_strain)%>%
  dplyr::mutate(years_since_emergence=year-min(year))

# filling in which are the years with epidemics and their sizes
cumulative_incidence_by_ag<-left_join(cumulative_incidence_by_ag,epi_table_no_corrections[,c("city","subtype","reference_strain","year","incidence_per_mil")])
cumulative_incidence_by_ag<-rename(cumulative_incidence_by_ag,replace=c("incidence_per_mil"= "incidence_current_season"))
cumulative_incidence_by_ag<-data.frame(cumulative_incidence_by_ag)
cumulative_incidence_by_ag$incidence_current_season[which(is.na(cumulative_incidence_by_ag$incidence_current_season))]<-0

# add cumulative incidence for each variant
cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::arrange(year,.by_group=TRUE)%>%
  dplyr::mutate(cumulative_prior_incidence_same_ag = cumsum(incidence_current_season)- incidence_current_season)

# need to divide cumulative incidence by city-specific mean epidemic size
overall_mean_size<-epi_table_no_corrections%>%
  subset(.,epi_alarm=="Y")%>%
  dplyr::group_by(city)%>%
  dplyr::summarise(mean_epi_size=mean(incidence_per_mil))

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%left_join(overall_mean_size)

cumulative_incidence_by_ag$standardised_prior_cumulative<-cumulative_incidence_by_ag$cumulative_prior_incidence_same_ag/cumulative_incidence_by_ag$mean_epi_size
cumulative_incidence_by_ag$standardised_current_season<-cumulative_incidence_by_ag$incidence_current_season/cumulative_incidence_by_ag$mean_epi_size

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%left_join(data.frame(city=epi_table_no_corrections$city,reference_strain=epi_table_no_corrections$reference_strain,year=epi_table_no_corrections$year,epi_alarm=factor(epi_table_no_corrections$epi_alarm,levels = c("Y","N"))))
cumulative_incidence_by_ag$epi_alarm[is.na(cumulative_incidence_by_ag$epi_alarm)==TRUE]<-"N"

cumulative_incidence_by_ag$epi_alarm2<-cumulative_incidence_by_ag$epi_alarm%>%
  mapvalues(.,from=c("Y","N"),to=c(1,0))


epi_size_cumulative_size_same_variant_plot<-cumulative_incidence_by_ag %>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x=standardised_prior_cumulative,y= standardised_current_season))+
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=3.5)+
  geom_smooth(method='lm',formula=y~x, se = FALSE)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  xlab("Antigenic Variant-Specific Cumulative Incidence")+
  ylab("Epidemic Incidence")+
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=25),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype)


# 2) Logistics regression plot: effect of cumulative incidence on p(successful epidemic) --------
prob_successful_epi_cumulative_size_same_variant_plot<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x=standardised_prior_cumulative,y= as.numeric(as.character(epi_alarm2))))+
  
  geom_jitter(aes(group = subtype, color=subtype),
              position=position_jitter(width=0.05,height=0.05),alpha=0.6,size=3.5)+
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0,1))+
  scale_color_manual(name = "Subtype",
                     values=c("B/Yam"="#CC79A7",
                              "B/Vic"="#009E73",
                              "H1sea"="#56B4E9",
                              "H1pdm09"="#999999",
                              "H3"="#E69F00"))+
  xlab("Antigenic Variant-Specific Cumulative Incidence")+
  ylab("Epidemic Present/Absent")+
  theme_bw()+
  theme(axis.title=element_text(size=16),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+facet_grid(.~subtype)+ 
  theme(strip.text = element_text(size=25))


# 3) Logistics regression and Odds Ratios: effect of cumulative incidence on p(successful epidemic) ------------------------------------
subtype_logistics_regression<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(subtype)%>%
  dplyr::summarise(OR= glm(as.numeric(as.character(epi_alarm2))~standardised_prior_cumulative,family=binomial)%>%coef(.)%>%.[2]%>%exp(.),
                   SE= glm(as.numeric(as.character(epi_alarm2))~standardised_prior_cumulative,family=binomial)%>%summary(.)%>%coef(.)%>%.[2,2],
                   p_value= glm(as.numeric(as.character(epi_alarm2))~standardised_prior_cumulative,family=binomial)%>%summary(.)%>%coef(.)%>%.[2,4])
subtype_logistics_regression$holm_adjusted_p<-p.adjust(subtype_logistics_regression$p_value,method = "holm")
