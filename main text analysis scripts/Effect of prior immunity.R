library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)
library(readr)

#Here we assume and correct for potential mis-identification during antigenic characterisation due to delays in 
#updating vaccine strain nomenclature.

# The following code will reproduce the prior immunity analyses discussed in the main text:
# 1)  The relationship between epidemic incidence and the amount of antigenic variant-specific cumulative incidence 
#     epi_size_cumulative_size_same_variant_plot (Figure 4)
#
# 2)  The relationship between the probability of successful epidemic initiation 
#     and amount of antigenic variant-specific cumulative incidence for each subtype 
#     prob_successful_epi_cumulative_size_same_variant_plot (Figure S8)
#
# 3)  Binary logistic regression assessing the effect of antigenic variant-specific cumulative incidence 
#     on the probability of successful epidemic initiation for each subtype 
#     subtype_logistics_regression (Table S5)

# Loading data ------------------------------------------------------------

epi_table<-read.csv("./dat/raw/epi_table.csv")

cities<-c("ADELAIDE","BRISBANE","MELBOURNE","PERTH","SYDNEY")

epi_table$city<-factor(epi_table$city,levels = cities)


#normalise incidence  log(incidence )
epi_table<-epi_table%>%
  dplyr::mutate(log_incidence = log(incidence_per_mil))

epi_table<-epi_table%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_city = ifelse(epi_alarm=="Y",
                                                (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                NA),
                scaled_incidence_city = log_incidence-mean(log_incidence,na.rm=TRUE))

epi_table<-epi_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::mutate(incidence_z_score_subtype_city = ifelse(epi_alarm=="Y",
                                                        (log_incidence-mean(log_incidence,na.rm=TRUE))/sd(log_incidence,na.rm = TRUE),
                                                        NA),
                scaled_incidence_subtype_city = log_incidence - mean(log_incidence,na.rm=TRUE))


# 1) Effect of ag variant specific cumulative incidence on epidemic incidence --------

# new_ag_marker = NA means that the variant emerged prior to start of study period, 
# ie we are unable to calculate the cumulative incidence
# also, due to the pandemic H1N1 emergence in 2009, epidemic activity could not be reliable quantified for A/H3/Perth/16/2009
# B/Wisconsin/1/2010-like since it never caused an epidemic at any time.
#these need to be removed
#
dodgy_prev<-epi_table%>%subset(.,is.na(new_ag_marker) | 
                                 reference_strain%in%c("A/Perth/16/2009-like","B/Wisconsin/1/2010-like"))

# First find the first and last year in which an ag variant was detected in each city
first_emerge_table<-epi_table%>%
  #subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::summarise(first_detected=min(year),
                   last_recorded_epi=max(year))

# assume that an ag variant could potentially have caused an epidemic right up to the year 
# before the emergence of its replacement new variant
first_emerge_table<-first_emerge_table%>%
  dplyr::group_by(city,subtype)%>%
  dplyr::arrange(first_detected,.by_group=TRUE)%>%
  dplyr::mutate(last_possible_year=if(subtype!="H1sea"){lead(first_detected-1,default = 2015)}else{lead(first_detected-1,default = 2008)})

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
  tidyr::expand(.,year=full_seq(c(first_detected,last_possible_year),1))

cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,reference_strain)%>%
  dplyr::mutate(years_since_emergence=year-min(year))

# filling in which are the years with epidemics and their sizes
cumulative_incidence_by_ag<-left_join(cumulative_incidence_by_ag,epi_table[,c("city","subtype","reference_strain","year","epi_alarm","incidence_per_mil")])
cumulative_incidence_by_ag<-rename(cumulative_incidence_by_ag,replace=c("incidence_per_mil"= "incidence_current_season"))

cumulative_incidence_by_ag<-data.frame(cumulative_incidence_by_ag)
cumulative_incidence_by_ag$incidence_current_season[which(is.na(cumulative_incidence_by_ag$incidence_current_season))]<-0

# add cumulative incidence for each variant
cumulative_incidence_by_ag<-cumulative_incidence_by_ag%>%
  dplyr::group_by(city,subtype,reference_strain)%>%
  dplyr::arrange(year,.by_group=TRUE)%>%
  dplyr::mutate(cumulative_prior_incidence_same_ag = cumsum(incidence_current_season)- incidence_current_season)


size_plot_table<-cumulative_incidence_by_ag%>%
  subset(.,epi_alarm=="Y")

#cumulative incidence measured in terms of number of mean epidemics
mean_epi_size1<-epi_table%>%
  subset(.,epi_alarm=="Y" &year!=2009)%>%
  dplyr::group_by(city)%>%
  dplyr::summarise(mean_epi_size_c=mean(incidence_per_mil))

size_plot_table<-left_join(size_plot_table,mean_epi_size1)

mean_epi_size2<-epi_table%>%
  subset(.,epi_alarm=="Y" &year!=2009)%>%
  dplyr::group_by(subtype,city)%>%
  dplyr::summarise(mean_epi_size_sc=mean(incidence_per_mil))
size_plot_table<-left_join(size_plot_table,mean_epi_size2)

size_plot_table<-size_plot_table%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(z_score_incidence_current_season_c = ifelse(epi_alarm=="Y",
                                                            (log(incidence_current_season)-mean(log(incidence_current_season),na.rm=TRUE))/sd(log(incidence_current_season),na.rm = TRUE),
                                                            NA),
                scaled_incidence_current_season_c= log(incidence_current_season)-mean(log(incidence_current_season),na.rm=TRUE),
                scaled_cumulative_incidence_c = cumulative_prior_incidence_same_ag/mean_epi_size_c)

size_plot_table<-size_plot_table%>%
  dplyr::group_by(subtype,city)%>%
  dplyr::mutate(z_score_incidence_current_season_sc = ifelse(epi_alarm=="Y",
                                                            (log(incidence_current_season)-mean(log(incidence_current_season),na.rm=TRUE))/sd(log(incidence_current_season),na.rm = TRUE),
                                                            NA),
                scaled_incidence_current_season_sc= log(incidence_current_season)-mean(log(incidence_current_season),na.rm=TRUE),
                scaled_cumulative_incidence_sc = cumulative_prior_incidence_same_ag/mean_epi_size_sc)


cumulative_incidence_by_ag<-left_join(cumulative_incidence_by_ag,mean_epi_size1)%>%
  left_join(.,mean_epi_size2)%>%
  dplyr::mutate(scaled_cumulative_incidence_c = cumulative_prior_incidence_same_ag/mean_epi_size_c,
                scaled_cumulative_incidence_sc = cumulative_prior_incidence_same_ag/mean_epi_size_sc)
  
cumulative_incidence_by_ag$epi_alarm[is.na(cumulative_incidence_by_ag$epi_alarm)==TRUE]<-"N"

cumulative_incidence_by_ag$epi_alarm2<-cumulative_incidence_by_ag$epi_alarm%>%
  mapvalues(.,from=c("Y","N"),to=c(1,0))


# plot 1 ------------------------------------------------------------------
epi_size_cumulative_size_same_variant_plot<-size_plot_table %>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x= scaled_cumulative_incidence_c ,y=scaled_incidence_current_season_c ))+
  geom_jitter(aes(group = subtype, color=city),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  xlab("Antigenic variant-specific cumulative incidence")+
  ylab("Epidemic incidence")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype)

pooled_epi_size_cumulative_size_same_variant_plot<-size_plot_table %>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x=scaled_cumulative_incidence_sc,y= scaled_incidence_current_season_sc))+
  geom_jitter(aes(group = subtype, color=city),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  xlab("Antigenic variant-specific cumulative incidence")+
  ylab("Epidemic incidence")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 2) Logistics regression plot: effect of cumulative incidence on p(successful epidemic) --------
prob_successful_epi_cumulative_size_same_variant_plot<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x=scaled_cumulative_incidence_c,y= as.numeric(as.character(epi_alarm2))))+

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
  xlab("Antigenic variant-specific cumulative incidence")+
  ylab("Epidemic Present/Absent")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_grid(~subtype)


# 3) Logistic regression and Odds Ratios: effect of cumulative incidence on p(successful epidemic) ------------------------------------
subtype_logistic_regression<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(subtype)%>%
  dplyr::summarise(term = "Cumulative incidence",
                   OR= glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial)%>%coef(.)%>%.[2]%>%exp(.),
                   OR_adjusted_SE = sqrt(OR^2 * diag(vcov(glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial))))%>%.[2],
                   p_value= glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial)%>%summary(.)%>%coef(.)%>%.[2,4])
subtype_logistic_regression$holm_adjusted_p<-p.adjust(subtype_logistic_regression$p_value,method = "holm")


subtype_logistic_regression2<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  dplyr::group_by(subtype)%>%
  dplyr::summarise(term = "(intercept)",
                   OR= glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial)%>%coef(.)%>%.[1]%>%exp(.),
                   OR_adjusted_SE = sqrt(OR^2 * diag(vcov(glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial))))%>%.[1],
                   p_value= glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial)%>%summary(.)%>%coef(.)%>%.[1,4])
subtype_logistic_regression2$holm_adjusted_p<-p.adjust(subtype_logistic_regression2$p_value,method = "holm")

subtype_logistic_output<-rbind(subtype_logistic_regression,subtype_logistic_regression2)
subtype_logistic_output<-dplyr::arrange(subtype_logistic_output,subtype,term)

# 4) AGGREGATED ACROSS SUBTYPES: Effect of ag variant specific cumulative incidence on epidemic incidence --------

aggregated_cumulative_size_plot<-size_plot_table %>%
  subset(.,subtype!="H1pdm09")%>%
  subset(.,epi_alarm=="Y")%>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x= scaled_cumulative_incidence_sc ,y=scaled_incidence_current_season_c ))+
  geom_jitter(aes(group = subtype, color=city),
              position=position_jitter(width=0.01,height=0.01),alpha=0.6,size=5)+
  stat_cor(method = "pearson",size=8)+
  scale_color_manual(name = "City",
                     values=c("ADELAIDE"="#CC79A7",
                              "BRISBANE"="#009E73",
                              "MELBOURNE"="#56B4E9",
                              "PERTH"="#999999",
                              "SYDNEY"="#E69F00"))+
  xlab("Antigenic variant-specific cumulative incidence")+
  ylab("Epidemic incidence")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 5) Aggregated Logistics regression plot: effect of cumulative incidence on p(successful epidemic) --------
aggregated_prob_successful_epi_plot<-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  ggplot(.,aes(x=scaled_cumulative_incidence_sc,y= as.numeric(as.character(epi_alarm2))))+
  
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
  xlab("Antigenic variant-specific cumulative incidence")+
  ylab("Epidemic Present/Absent")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=25),
        axis.title=element_text(size=25),
        axis.text.x =element_text(size=20,margin=margin(t=7,r=0,b=0,l=0)),
        axis.text.y =element_text(size=20,margin=margin(t=0,r=7,b=0,l=0)),
        axis.ticks.length = unit(0.4,"cm"),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



aggregated_logistic_regression <-cumulative_incidence_by_ag %>%
  subset(.,!(reference_strain%in%dodgy_prev$reference_strain))%>%
  glm(as.numeric(as.character(epi_alarm2))~scaled_cumulative_incidence_c,family=binomial,data=.)

aggregated_logistic_output<-data.frame(OR = aggregated_logistic_regression%>%coef(.)%>%exp(.),
                                       OR_adjusted_SE = sqrt((aggregated_logistic_regression%>%coef(.)%>%exp(.))^2 * diag(vcov(aggregated_logistic_regression))),
                                       p_value = aggregated_logistic_regression%>%summary(.)%>%coef(.)%>%.[,4])

aggregated_logistic_output$holm_adjusted_p<-p.adjust(aggregated_logistic_output$p_value,method = "holm")
aggregated_logistic_output$term<-row.names(aggregated_logistic_output)


# save plots --------------------------------------------------------------


ggsave(plot = epi_size_cumulative_size_same_variant_plot,"./figures/main/figure_4.png",
       width=20, height=8,limitsize=FALSE)

ggsave(plot = epi_size_cumulative_size_same_variant_plot,"./figures/main/figure_4.eps",
       device = cairo_ps,
       width=20, height=8,limitsize=FALSE)


ggsave(plot = prob_successful_epi_cumulative_size_same_variant_plot,"./figures/supp/figure_S8.png",
       width=20, height=8,limitsize=FALSE)

write.csv(subtype_logistic_output%>%dplyr::mutate_if(is.numeric,signif,digits=3),"./tables/table_S5.csv",row.names = FALSE)



stop()

# aggregated plots --------------------------------------------------------

ggsave(plot=aggregated_cumulative_size_plot,"./figures/aggregated_possibles/agg_cumulative_size.png",
       width=8, height=8,limitsize=FALSE)

ggsave(plot = aggregated_prob_successful_epi_plot,"./figures/aggregated_possibles/agg_prob_success.png",
       width=8, height=8,limitsize=FALSE)

write.csv(aggregated_logistic_output%>%dplyr::mutate_if(is.numeric,signif,digits=3),"./figures/aggregated_possibles/prob_success.csv",row.names = FALSE)
