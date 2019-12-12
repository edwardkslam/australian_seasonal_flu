library(lubridate)
library(lme4)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)


# load data ---------------------------------------------------------------
raw_table<-read.csv("C:/Users/el382/Dropbox/PhD/code for manuscript/australian_seasonal_flu/raw_data.csv")


# initial processing ------------------------------------------------------
#add fortnights since start of the year
raw_table<-raw_table%>%
  dplyr::mutate(fortnights_since_start_of_year = yday(specimen_date)%/%14+1)

#list all the years and cities for which an ag variant has cases detected in
#and number of detected cases in that year
unique_seasons<-raw_table%>%
  dplyr::group_by(city,year,subtype,assumed_antigenic_variant)%>%
  dplyr::summarise(year_count = n())

#rename assumed antigenic variant as reference_strain from this point forwards
unique_seasons <- unique_seasons%>%
  dplyr::rename(reference_strain = assumed_antigenic_variant)

#create strain_year column for ease of searching for ag variant in a specific year across all cities
unique_seasons<-unique_seasons%>%
  dplyr::mutate(strain_year = paste(year,reference_strain,sep="_"))


unique_seasons<-unique_seasons%>%
  dplyr::group_by(city,year)%>%
  dplyr::mutate(year_fraction = year_count/sum(year_count))


# define functions --------------------------------------------------------

#takes a reference strain, year, city and returns the counts by fortnight
make_year_time_series<-function(x,print_plot="N"){
  x<-as.data.frame(x)
  temp<-raw_table%>%subset(.,city==x$city & year==x$year & assumed_antigenic_variant== x$reference_strain)
  temp<-temp%>%
    dplyr::group_by(year,fortnights_since_start_of_year)%>%
    dplyr::summarise(fortnight_count = n())
  
  fortnight_table<-data.frame(fortnights_since_start_of_year=c(1:26))
  
  temp<-left_join(fortnight_table,temp%>%subset(.,select=c(fortnights_since_start_of_year,fortnight_count)))
  
  na_row<-which(is.na(temp$fortnight_count))
  temp$fortnight_count[na_row]<-0
  
  if(print_plot=="Y"){
    temp_plot<-temp%>%ggplot(.,aes(x=fortnights_since_start_of_year,y=fortnight_count))+
      geom_point()+
      geom_line()+
      ggtitle(label = paste(x$city," " , x$year, "\n", x$subtype,"\n",x$reference_strain,sep=""))
    print(temp_plot)
  }
  return(temp)
}


# start detection ---------------------------------------------------------



unique_seasons<-raw_table%>%
  subset(.,select=c(city,year,subtype,assumed_antigenic_variant))%>%distinct(.,city,year,subtype,assumed_antigenic_variant)

#rename assumed antigenic variant as reference_strain from this point forwards
unique_seasons <- unique_seasons%>%
  dplyr::rename(reference_strain = assumed_antigenic_variant)

#create strain_year column for ease of searching for ag variant in a specific year across all cities
unique_seasons<-unique_seasons%>%
  dplyr::mutate(strain_year = paste(year,reference_strain,sep="_"))


#total number of cases per year
raw_table%>%dplyr::group_by(city,s)