''' import main dataset and country-specific ones extraction  '''
events <- read.csv('coronavirus_dec.csv', header=T)
events <- subset(events, EVENT_TYPE!="Strategic developments")
events$id<-seq.int(nrow(events))

mexico <- events[events$COUNTRY=='Mexico',]


mexico_events<-mexico[c("LATITUDE", "LONGITUDE", "EVENT_DATE", "id")]

## k-means
# Elbow method
fviz_nbclust(mexico_events[c("LATITUDE", "LONGITUDE")], kmeans, method = "wss", nboot=100) +
  ggtitle(label="")+
  theme(axis.text = element_text(size = 24, color = "black"))+
  theme(axis.title.x = element_text(size=25))+
  theme(axis.title.y = element_text(size=20))+
  geom_vline(xintercept = 4, linetype = 2)
#labs(subtitle = "Elbow method")



mexico_kmeans <- kmeans(mexico_events[c("LATITUDE", "LONGITUDE")], 4, nstart=25) # RUN PROPER K-MEANS
print(mexico_kmeans)


### visualize cluster in a 2d space

fviz_cluster(mexico_kmeans, data = mexico_events[c("LATITUDE", "LONGITUDE")], geom="point")+
  theme_bw()+
  ggtitle(label="")+
  theme(axis.text = element_text(size = 15, color = "black"))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
  
# connect retrieved cluster with the mexico dataset  
  
mexico_events <- cbind(mexico_events, cluster=mexico_kmeans$cluster)


# retrieve and plot map with events colored by cluster 

mexico_map <- get_map(c(left = -120, bottom =10, right = -86, top = 33.5), maptype="roadmap")
mexico_map <- ggmap(mexico_map)

cbp2 <- c("#000000", "#FF0000", "#0000FF", "#008080",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mexico_map + geom_point(data=mexico_events,aes(x=LONGITUDE, y=LATITUDE, color=as.factor(cluster)), size=1.8, alpha=0.8)+
  scale_colour_manual(values=cbp2, name = "Cluster")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title=element_text(size=22),
        legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18))+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_colour_discrete("Cluster")
