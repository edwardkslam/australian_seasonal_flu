unique_list02%>%
  subset(.,year_counts<100)%>%
  ggplot(.,aes(x=year_counts,color=epi_alarm))+
  geom_histogram(alpha=0.1,binwidth = 5)


unique_list005%>%
  subset(.,year_counts<100)%>%
  ggplot(.,aes(x=year_counts,color=epi_alarm))+
  geom_histogram(alpha=0.1,binwidth = 5)

sensitivity_df<-rbind(unique_list005%>%dplyr::mutate(alpha=0.05),
                      unique_list012%>%dplyr::mutate(alpha=0.12),
                      unique_list02%>%dplyr::mutate(alpha=0.2))


sensitivity_df%>%
  subset(.,year_counts<100)%>%
  ggplot(.,aes(x=year_counts,color=epi_alarm,fill=epi_alarm))+
  geom_histogram(alpha=0.1,binwidth = 5)+
  ggtitle("Annual counts")+
  facet_grid(~alpha)


sensitivity_df%>%
  subset(.,epi_alarm=="Y")%>%
  ggplot(.,aes(x=start))+
  geom_histogram(binwidth = 1)+
  ggtitle("Start time")+
  facet_grid(~alpha)
