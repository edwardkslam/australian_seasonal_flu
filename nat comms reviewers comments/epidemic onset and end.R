library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)

library(surveillance)
library(sleekts)
library(lubridate)


# loading data ------------------------------------------------------------

raw_table<-read.csv("./dat/raw/raw_data.csv")
raw_table<-raw_table%>%dplyr::mutate(weeks_since_start_of_year = specimen_date%>%ymd(.)%>%week(.),
                                     weeks_since_start_of_year = ifelse(weeks_since_start_of_year==53,52,weeks_since_start_of_year),
                                     fortnights_since_start_of_year = (as.numeric(weeks_since_start_of_year)/2)%>%ceiling(.))

raw_table<-raw_table%>%subset(.,!(assumed_antigenic_variant=="B/Florida/4/2006-like" & year%in%c(2004,2005)))

#test<-raw_table%>%subset(.,city=="ADELAIDE" & assumed_antigenic_variant =="A/Moscow/10/99-like" & year==2000)

city_erp<-read.csv("./dat/raw/ERP.csv")

# detection algorithms ----------------------------------------------------

start_time_finder<-function(c,y,ag_variant,base = "fortnight",
                            window=3,a=0.12,smooth = TRUE,
                            debug=FALSE,return_start_threshold=FALSE){
  temp_ts<-raw_table%>%subset(.,city==c & year==y& assumed_antigenic_variant==ag_variant)
  if(base=="fortnight"){
    temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
    
    if(smooth==TRUE){
      temp_sleek<-sleek(temp_ts)
      #print(temp_sleek)
      temp_sleek[which(temp_ts==0)]<-0
    }
    
    forward_temp_sts<-sts(temp_sleek,frequency = 26)
    temp_alarm<-algo.bayes(forward_temp_sts%>%sts2disProg(),control = list(range=c(window+4:26),b=0,w=window,alpha=a))
    if(debug==TRUE){
      print(temp_ts)
      print(temp_sleek)
      return(temp_alarm)
    }
    alarm_ts<-c(rep(0,window+3),temp_alarm$alarm)
    alarm_ts[1:6]<-0
    
    alarm_times<-which(alarm_ts==1)
    
    if(length(alarm_times)==0){
      return(NA)
    }
       
    if(length(alarm_times)>1 & (alarm_times[1]-alarm_times[2])!=1){
      if(sum(temp_ts[alarm_times[1]:alarm_times[2]]==0)!=0){
        if(return_start_threshold=="TRUE"){
          return(c(rep(0,window+3),temp_alarm$upperbound)[alarm_times[2]])
        }
        return(alarm_times[2])
      }
    }
    output<-ifelse(sum(alarm_ts,na.rm = TRUE)==0,NA,alarm_ts%>%which.max)
    
    if(!is.na(output) & return_start_threshold=="TRUE"){
      return(c(rep(0,window+3),temp_alarm$upperbound)[output])
    }
    return(output)
  }
  if(base=="weekly"){
    temp_ts<-temp_ts%>%.$weeks_since_start_of_year%>%tabulate(.,52)
    
    if(smooth==TRUE){
      temp_sleek<-sleek(temp_ts)
      temp_sleek[which(temp_ts==0)]<-0
    }
    
    forward_temp_sts<-sts(temp_sleek,frequency = 52)
    temp_alarm<-algo.bayes(forward_temp_sts%>%sts2disProg(),control = list(range=c(window+1:52),b=0,w=window,alpha=a))
    
    alarm_ts<-c(rep(0,window),temp_alarm$alarm)
    alarm_ts[1:8]<-0
  }
  
  output<-ifelse(sum(alarm_ts,na.rm = TRUE)==0,NA,alarm_ts%>%which.max)
  return(output)
}

end_time_finder<-function(c,y,ag_variant,base = "fortnight",
                            window=3,a=0.12,smooth = TRUE){
  temp_ts<-raw_table%>%subset(.,city==c & year==y& assumed_antigenic_variant ==ag_variant)
  if(base=="fortnight"){
    temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
    temp_ts<-rev(temp_ts)
    
    if(smooth==TRUE){
      temp_sleek<-sleek(temp_ts)
      temp_sleek[which(temp_sleek<0)]<-0
      #temp_sleek[which(temp_ts==0)]<-0
    }
    
    rev_temp_sts<-sts(temp_sleek,frequency = 26)
    temp_alarm<-algo.bayes(rev_temp_sts%>%sts2disProg(),control = list(range=c(window+1:26),b=0,w=window,alpha=a))
    
    alarm_ts<-c(rep(0,window),temp_alarm$alarm)
    
    output<-ifelse(sum(alarm_ts,na.rm = TRUE)==0,NA,27-(alarm_ts%>%which.max))
    return(output)
  }
  if(base=="weekly"){
    temp_ts<-temp_ts%>%.$weeks_since_start_of_year%>%tabulate(.,52)
    temp_ts<-rev(temp_ts)
    
    if(smooth==TRUE){
      temp_sleek<-sleek(temp_ts)
      temp_sleek[which(temp_sleek<0)]<-0
      #temp_sleek[which(temp_ts==0)]<-0
    }
    
    rev_temp_sts<-sts(temp_sleek,frequency = 52)
    temp_alarm<-algo.bayes(rev_temp_sts%>%sts2disProg(),control = list(range=c(window+1:52),b=0,w=window,alpha=a))
    
    alarm_ts<-c(rep(0,window),temp_alarm$alarm)
    
    output<-ifelse(sum(alarm_ts,na.rm = TRUE)==0,NA,53-(alarm_ts%>%which.max))
    return(output)
  }

}

epi_finder<-function(c,y,s,ag_variant,base = "fortnight",
                     window1=3,a1=0.12,smooth1 = TRUE,
                     window2=3,a2=0.12,smooth2 = TRUE){
  temp_start<-start_time_finder(c,y,ag_variant,base,window1,a1,smooth1)
  temp_end<-end_time_finder(c,y,ag_variant,base,window2,a2,smooth2)
  
  temp_ts<-raw_table%>%subset(.,city==c & year==y& assumed_antigenic_variant ==ag_variant)
  if(base=="fortnight"){
    temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
  }
  if(base=="weekly"){
    temp_ts<-temp_ts%>%.$weeks_since_start_of_year%>%tabulate(.,52)
  }
  
  temp_output<-data.frame(city=c,
                          year=y%>%as.numeric,
                          subtype=s%>%as.character(),
                          reference_strain=ag_variant,
                          start=temp_start,
                          end=temp_end,
                          year_counts=sum(temp_ts))

  temp_output<-temp_output%>%dplyr::mutate(raw_year_fraction=year_counts/nrow(raw_table%>%subset(.,city==c & year==y)))
  
  peak<-which.max(temp_ts)
  
  if(!is.na(temp_output$start) & is.na(temp_output$end) & temp_output$raw_year_fraction>0.1){
    temp_output$fail_end<-"Y"
    
    #first id the first zero after peak
    zeroes<-which(temp_ts==0)

    pos_end1<-zeroes[which(zeroes>temp_start & zeroes>which.max(temp_ts))]
    if(length(pos_end1)!=0){
      pos_end1<-pos_end1%>%min()
      pos_end1<-pos_end1-1
    }else{
      pos_end1<-NA
    }
    
    #use start as threshold for end, if smoothening has prevented detection of end of epidemic
    temp_start_threshold<-start_time_finder(c,y,ag_variant,base,window1,a1,smooth1,return_start_threshold = TRUE)
    sleek_ts<-temp_ts%>%sleek()
    sleek_ts[which(temp_ts==0)]<-0
    
    rev_sleek_ts<-temp_ts%>%rev()%>%sleek()
    above_times<-which(rev_sleek_ts>=temp_start_threshold)#sleek_ts[temp_start])

    if(base=="fortnight"){
      above_times<-27-above_times
    }
    if(base=="weekly"){
      above_times<-53-above_times
    }
    
    #above_times<-which(temp_ts>=temp_ts[temp_start])#temp_start_threshold)
    if(!is.na(pos_end1)){
      pos_end2<-above_times[above_times>temp_start & above_times<=pos_end1]%>%max(.)
    }else{
      pos_end2<-above_times[above_times>temp_start & above_times>=peak]%>%max(.)
    }
    
    #first point at which raw counts drop back to start
    below_times<-which(temp_ts<=temp_ts[temp_start])
    pos_end3<-below_times[below_times>temp_start & below_times >=peak]%>%min()
    
    
    below_times_smooth<-which(rev_sleek_ts<=sleek_ts[temp_start])
    pos_end4<-below_times[below_times>temp_start & below_times >=peak]%>%min()
    
    # temp_output$pos_end1<-pos_end1
    # temp_output$pos_end2<-pos_end2
    # temp_output$pos_end3<-pos_end3
    # temp_output$pos_end4<-pos_end4
    pos_end_vec<-c(pos_end1,pos_end2,pos_end3,pos_end4)
    pos_end_vec<-pos_end_vec[is.finite(pos_end_vec)&!is.na(pos_end_vec)]
    min_end<-min(pos_end_vec)
    

    temp_output$end<-ifelse(min_end<=peak,NA,min_end)
  }else{
    temp_output$fail_end<-"N"
  }
  
  temp_erp<-city_erp%>%subset(.,city==c & year==y)%>%.$ERP
  
  temp_output<-temp_output%>%
    dplyr::mutate(peak=peak,
                  epi_alarm=ifelse(year==2009 | is.na(start)|is.na(end)| end<start | raw_year_fraction < 0.05,"N","Y"),
                  epi_counts = ifelse(epi_alarm=="Y",sum(temp_ts[start:end]),NA),
                  incidence_per_mil = epi_counts/temp_erp*10^6)
                                           
  return(temp_output)
}

ts_plotter<-function(c,y,ag_variant,base = "fortnight",window=3,a=0.12){
  temp_ts<-raw_table%>%subset(.,city==c & year==y& assumed_antigenic_variant ==ag_variant)
  if(base=="fortnight"){
    temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
    temp_ts<-data.frame(fortnight=c(1:26),count=temp_ts)
    forward_temp_sts<-sts(sleek(temp_ts$count),frequency = 52)
    
    temp_alarm<-algo.bayes(forward_temp_sts%>%sts2disProg(),control = list(range=c(window+1:26),b=0,w=window,alpha=a))
    alarm_ts<-data.frame(fortnight= c(1:26),count=c(rep(0,window),temp_alarm$upperbound)[1:26])

    temp_plot<-ggplot(data=temp_ts,aes(x=fortnight,y=count))+
      geom_point()+
      geom_line(data=temp_ts,aes(x=fortnight,y=sleek(count)))+
      
      geom_line(data=alarm_ts,aes(x=fortnight,y=count),colour="red")+
      scale_x_continuous(breaks = seq(1,26,1), limits = c(0.8,26.2))+
      ggtitle(paste(c,y,ag_variant,sep=" "))+
      theme(strip.background = element_blank())
    
    show(temp_plot)
  }
  return(temp_alarm)
}

manual_corrector<-function(summary_df,corrections_df,base="fortnight"){
  for(i in 1:nrow(corrections_df)){

    #modify start and end times
    temp_index<-which(summary_df$city==as.character(corrections_df$city[i]) &
                        summary_df$year==corrections_df$year[i] &
                        summary_df$reference_strain==as.character(corrections_df$reference_strain[i]))
    print(temp_index)
    if(!is.na(corrections_df$start[i])){
      summary_df$start[temp_index]<-corrections_df$start[i]
    }
    if(!is.na(corrections_df$end[i])){
      summary_df$end[temp_index]<-corrections_df$end[i]
    }

    summary_df$manual_mod[temp_index]<-1
    summary_df$epi_alarm[temp_index]<-"Y"
    
    #recalculate epi incidence
    temp_ts<-raw_table%>%subset(.,city==as.character(corrections_df$city[i]) &
                                  year==corrections_df$year[i]& 
                                  assumed_antigenic_variant ==as.character(corrections_df$reference_strain[i]))
    if(base=="fortnight"){
      temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
    }
    if(base=="weekly"){
      temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,52)
    }
    
    temp_erp<-city_erp%>%subset(.,city==as.character(corrections_df$city[i]) & year==corrections_df$year[i])%>%.$ERP
    
    summary_df$epi_counts[temp_index]<-sum(temp_ts[summary_df$start[temp_index]:summary_df$end[temp_index]])
    summary_df$incidence_per_mil[temp_index]<-summary_df$epi_counts[temp_index]/temp_erp*10^6
  }
  return(summary_df)
}

count_prior_activity<-function(c,y,ag_variant,st,d,base="fortnight"){
  if(d=="N"){
    return(0)
  }else{
    temp_others<-unique_list%>%subset(city==c & year==y & reference_strain!=ag_variant & epi_alarm=="Y" & start<st)
    if(nrow(temp_others)!=0){
      other_counts<-c()
      for(i in 1:nrow(temp_others)){
        temp_strain<-temp_others[i,]
        
        temp_ts<-raw_table%>%subset(.,city==c & year==y& assumed_antigenic_variant ==temp_strain$reference_strain)
        if(base=="fortnight"){
          temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,26)
        }
        if(base=="weekly"){
          temp_ts<-temp_ts%>%.$fortnights_since_start_of_year%>%tabulate(.,52)
        }
        other_counts<-c(other_counts,sum(temp_ts[temp_strain$start:min(st,temp_strain$end)]))
        
      }
    }else{
      return(0)
    }
    
    return(sum(other_counts))
  }
  
}



# unique year-ag variants --------------------------------------------------------
# take into account assumed phylogenetic corrections
unique_list<-raw_table%>%dplyr::group_by(city,year,subtype,assumed_antigenic_variant)%>%
  dplyr::summarise(year_counts=n())%>%
  dplyr::mutate(raw_year_fraction = year_counts/sum(year_counts))

colnames(unique_list)<-c("city","year","subtype","reference_strain","yc","raw_year_fraction")
unique_list$subtype<-factor(unique_list$subtype,levels=c("B/Vic","B/Yam","H1sea","H1pdm09","H3"))

unique_list<-do.call("rbind",apply(unique_list,1,function(x){return(epi_finder(c = x['city'],y=x['year'],s=x['subtype'],ag_variant = x['reference_strain']))}))
unique_list$subtype<-factor(unique_list$subtype,levels=c("B/Vic","B/Yam","H1sea","H1pdm09","H3"))

# manual corrections ------------------------------------------------------
unique_list$orig_start<-unique_list$start
unique_list$orig_end<-unique_list$end
unique_list$manual_mod<-0


manual_corrections<-data.frame(city=c("BRISBANE"),
                               year=c(2011),
                               reference_strain=c("A/California/7/2009-like"),
                               start=c(11),
                               end=c(NA))

unique_list<-manual_corrector(unique_list,manual_corrections)

unique_list<-unique_list%>%dplyr::group_by(city,year,epi_alarm)%>%
  dplyr::mutate(epi_fractional_counts = epi_counts/sum(epi_counts) )


# marking first detection of a new ag variant -----------------------------
unique_list<-unique_list%>%dplyr::group_by(city,reference_strain,epi_alarm)%>%
  dplyr::mutate(first_detection_of_new_ag = ifelse(year==min(year) ,1,0) )

new_ag_problem<-data.frame(year=c(2000,2000,2000,
                                  2009,2009),
                           reference_strain = c("A/Moscow/10/99-like","B/Sichuan/379/99-like","A/New Caledonia/20/99-like",
                                                "A/California/7/2009-like","A/Perth/16/2009-like"))
for(i in 1:nrow(new_ag_problem)){
  unique_list$first_detection_of_new_ag[which(as.character(unique_list$reference_strain)==as.character(new_ag_problem$reference_strain[i]) &  unique_list$year==new_ag_problem$year[i])]<-NA
  unique_list$first_detection_of_new_ag[which(as.character(unique_list$reference_strain)==as.character(new_ag_problem$reference_strain[i]) &  unique_list$year>new_ag_problem$year[i])]<-0
}




# marking the first epi of each newly emerged ag variant ------------------

unique_list<-unique_list%>%dplyr::group_by(city,reference_strain,epi_alarm)%>%
  dplyr::mutate(new_ag_marker = ifelse(year==min(year) & epi_alarm=="Y" ,1,0) )

#lack of data availability for 1999 prevents us from knowing if the 99 viruses are causing their first proper epi in 2000 or not
# and surveillance suspension in 2009 prevents us from accurately estimating epi size
new_ag_problem<-data.frame(year=c(2000,2000,2000,
                                  2009,2009),
                           reference_strain = c("A/Moscow/10/99-like","B/Sichuan/379/99-like","A/New Caledonia/20/99-like",
                                                "A/California/7/2009-like","A/Perth/16/2009-like"))

for(i in 1:nrow(new_ag_problem)){
  unique_list$new_ag_marker[which(as.character(unique_list$reference_strain)==as.character(new_ag_problem$reference_strain[i]) &  unique_list$year==new_ag_problem$year[i])]<-NA
  unique_list$new_ag_marker[which(as.character(unique_list$reference_strain)==as.character(new_ag_problem$reference_strain[i]) &  unique_list$year>new_ag_problem$year[i])]<-0
}


# marking the earliest and largest of each season -------------------------

unique_list<-unique_list%>%dplyr::group_by(city,year,epi_alarm)%>%
  dplyr::mutate(delayed = ifelse(epi_alarm=="Y",ifelse(start==min(start),"N","Y"),NA))

unique_list<-unique_list%>%dplyr::group_by(city,year,delayed,epi_alarm)%>%
  dplyr::mutate(first_n_biggest = ifelse(epi_alarm=="Y",ifelse(epi_counts==max(epi_counts)& delayed=="N","Y","N"),NA))

unique_list<-unique_list%>%dplyr::group_by(city,year,epi_alarm)%>%
  dplyr::mutate(relative_to_first_n_biggest = ifelse(epi_alarm=="Y",
                                                     ifelse(year!=2009,
                                                            epi_counts/epi_counts[first_n_biggest=="Y"],
                                                            NA),
                                                     NA),
                delay = ifelse(year!=2009 & epi_alarm=="Y",start-start[which(delayed=="N" & epi_alarm=="Y")],NA))

unique_list<-unique_list%>%
  rowwise(.)%>%
  dplyr::mutate(prior_everything = ifelse(epi_alarm=="Y",
                                          count_prior_activity(city,year,reference_strain,start,delayed),
                                          NA))

unique_list<-unique_list%>%
  dplyr::group_by(city)%>%
  dplyr::mutate(mean_epi_size=mean(epi_counts,na.rm=TRUE),
                prior_everything_scaled = prior_everything/mean_epi_size)



# formatting --------------------------------------------------------------

unique_list<-unique_list%>%arrange(city,year,epi_alarm,desc(epi_fractional_counts))
unique_list$strain_year<-paste(unique_list$year,unique_list$reference_strain,sep="_")

unique_list<-unique_list%>%select(.,city,strain_year,year,subtype,reference_strain,
                                  epi_alarm,start,end,epi_counts,incidence_per_mil,epi_fractional_counts,
                                  new_ag_marker,first_detection_of_new_ag,
                                  first_n_biggest,relative_to_first_n_biggest,
                                  delay,prior_everything_scaled,mean_epi_size)

write.csv(unique_list,file = "./dat/raw/epi_table.csv",row.names=FALSE)
