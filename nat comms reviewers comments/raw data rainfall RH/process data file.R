all_cities_climate<-read.csv("./nat comms reviewers comments/all_cities_climate.csv")
mean_fortnightly_climate_30years<-all_cities_climate%>%
  subset(.,year%in%c(1985:2015))%>%
  subset(.,fortnights_since_start_of_year!=27)%>%
  dplyr::group_by(city,year,fortnights_since_start_of_year)%>%
  dplyr::summarise(mean_AH = mean(Absolute_humidity,na.rm=TRUE),
                   mean_SH = mean(Specific_humidity,na.rm=TRUE),
                   mean_RH = mean(Average_relative_humidity,na.rm=TRUE),
                   mean_rainfall = mean(Total_precipitation,na.rm=TRUE),
                   mean_temp = mean(Mean_Temp,na.rm=TRUE))

mean_fortnightly_climate_30years<-mean_fortnightly_climate_30years%>%
  subset(.,year%in%c(1985:2015))%>%
  dplyr::group_by(city,fortnights_since_start_of_year)%>%
  dplyr::mutate(mean_AH_for_that_fortnight_of_year = mean(mean_AH,na.rm = TRUE),
                mean_SH_for_that_fortnight_of_year = mean(mean_SH,na.rm = TRUE),
                mean_RH_for_that_fortnight_of_year = mean(mean_RH,na.rm = TRUE),
                mean_rainfall_for_that_fortnight_of_year = mean(mean_rainfall,na.rm = TRUE),
                mean_temp_for_that_fortnight_of_year = mean(mean_temp,na.rm = TRUE))

mean_fortnightly_climate_30years<-mean_fortnightly_climate_30years%>%
  dplyr::mutate(d.AH = mean_AH - mean_AH_for_that_fortnight_of_year,
                d.SH = mean_SH - mean_SH_for_that_fortnight_of_year,
                d.RH = mean_RH - mean_RH_for_that_fortnight_of_year,
                d.rainfall = mean_rainfall - mean_rainfall_for_that_fortnight_of_year,
                d.temp = mean_temp - mean_temp_for_that_fortnight_of_year)

write.csv(mean_fortnightly_climate_30years,"C:/Users/el382/Dropbox/PhD/shaman_bootstrap/mean_fortnightly_climate_30years.csv",
          row.names = FALSE)



mean_weekly_climate_30years<-all_cities_climate%>%
  subset(.,year%in%c(1985:2015))%>%
  subset(.,weeks_since_start_of_year!=53)%>%
  dplyr::group_by(city,year,weeks_since_start_of_year)%>%
  dplyr::summarise(mean_AH = mean(Absolute_humidity,na.rm=TRUE),
                   mean_SH = mean(Specific_humidity,na.rm=TRUE),
                   mean_RH = mean(Average_relative_humidity,na.rm=TRUE),
                   mean_rainfall = mean(Total_precipitation,na.rm=TRUE),
                   mean_temp = mean(Mean_Temp,na.rm=TRUE))

mean_weekly_climate_30years<-mean_weekly_climate_30years%>%
  subset(.,year%in%c(1985:2015))%>%
  dplyr::group_by(city,weeks_since_start_of_year)%>%
  dplyr::mutate(mean_AH_for_that_week_of_year = mean(mean_AH,na.rm = TRUE),
                mean_SH_for_that_week_of_year = mean(mean_SH,na.rm = TRUE),
                mean_RH_for_that_week_of_year = mean(mean_RH,na.rm = TRUE),
                mean_rainfall_for_that_week_of_year = mean(mean_rainfall,na.rm = TRUE),
                mean_temp_for_that_week_of_year = mean(mean_temp,na.rm = TRUE))

mean_weekly_climate_30years<-mean_weekly_climate_30years%>%
  dplyr::mutate(d.AH = mean_AH - mean_AH_for_that_week_of_year,
                d.SH = mean_SH - mean_SH_for_that_week_of_year,
                d.RH = mean_RH - mean_RH_for_that_week_of_year,
                d.rainfall = mean_rainfall - mean_rainfall_for_that_week_of_year,
                d.temp = mean_temp - mean_temp_for_that_week_of_year)

write.csv(mean_weekly_climate_30years,"C:/Users/el382/Dropbox/PhD/shaman_bootstrap/mean_weekly_climate_30years.csv",
          row.names = FALSE)
