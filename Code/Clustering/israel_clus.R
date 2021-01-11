''' import main dataset and country-specific ones extraction  '''
events <- read.csv('coronavirus_dec.csv', header=T)
events <- subset(events, EVENT_TYPE!="Strategic developments")
events$id<-seq.int(nrow(events))

israel <- events[events$COUNTRY=='Israel',]

israel_events<-israel[c("LATITUDE", "LONGITUDE", "EVENT_DATE", "id")]


## k-means
# Elbow method
fviz_nbclust(israel_events[c("LATITUDE", "LONGITUDE")], kmeans, method = "wss") +
  ggtitle(label="")+
  theme(axis.text = element_text(size = 24, color = "black"))+
  theme(axis.title.x = element_text(size=25))+
  theme(axis.title.y = element_text(size=20))+
  geom_vline(xintercept = 4, linetype = 2)
  
  
  ''' Run k-means. IMPORTANT NOTE: due to the randomness of the centroid selection process, the cluster labels will be different every time the 
script is run, therefore one has to check that the cluster assignments match some specific characteristics that are
used in other scripts. '''

isr_kmeans <- kmeans(israel_events[c("LATITUDE", "LONGITUDE")], 4, nstart=25)
print(isr_kmeans)


### visualize cluster in a 2d space

fviz_cluster(isr_kmeans, data = israel_events[c("LONGITUDE", "LATITUDE")],geom="point")+
  theme_bw()+
  ggtitle(label="")+
  theme(axis.text = element_text(size = 15, color = "black"))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
  
  
# connect retrieved cluster with the israel dataset 
israel_events <- cbind(israel_events, cluster=isr_kmeans$cluster)


# retrieve and plot map with events colored by cluster 

israel_map2 <- get_map(c(left = 33, bottom =29.2, right = 36.5, top = 34), maptype="roadmap")
israel_map <-ggmap(israel_map2)

israel_map + geom_point(data=israel_events,aes(x=LONGITUDE, y=LATITUDE, color=as.factor(cluster)), size=1.5, alpha=0.8)+
  scale_colour_manual(values=cbp2, name = "Cluster")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title=element_text(size=22),
        legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18))+
  xlab("Longitude")+
  ylab("Latitude")
