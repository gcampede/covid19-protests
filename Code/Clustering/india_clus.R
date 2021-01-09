''' import main dataset and country-specific ones extraction  '''
events <- read.csv('coronavirus_dec.csv', header=T)
events <- subset(events, EVENT_TYPE!="Strategic developments")
events$id<-seq.int(nrow(events))

india <- events[events$COUNTRY=='India',]

india_events<-india[c("LATITUDE", "LONGITUDE", "EVENT_DATE", "id")]


## k-means
# Elbow method
fviz_nbclust(india_events[c("LATITUDE", "LONGITUDE")], kmeans, method = "wss") +
  ggtitle(label="")+
  theme(axis.text = element_text(size = 24, color = "black"))+
  theme(axis.title.x = element_text(size=25))+
  theme(axis.title.y = element_text(size=20))+
  geom_vline(xintercept = 4, linetype = 2)
  
  
  ''' Run k-means. IMPORTANT NOTE: due to the randomness of the centroid selection process, the cluster labels will be different every time the 
script is run, therefore one has to check that the cluster assignments match some specific characteristics that are
used in other scripts. '''

india_kmeans <- kmeans(india_events[c("LATITUDE", "LONGITUDE")], 4, nstart=25)
print(india_kmeans)

### visualize cluster in a 2d space

fviz_cluster(india_events[c("LATITUDE", "LONGITUDE")], data = india_events, geom="point")+
  theme_bw()+
  ggtitle(label="")+
  theme(axis.text = element_text(size = 15, color = "black"))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))
  
  
# connect retrieved cluster with the india dataset 

india_events <- cbind(india_events, cluster=india_kmeans$cluster)


# retrieve and plot map with events colored by cluster 

india_map2 <- get_map(c(left = 65, bottom =7, right = 97, top = 35.5), maptype="roadmap")
india_map <- ggmap(india_map2)
india_map


india_map + geom_point(data=india_events,aes(x=LONGITUDE, y=LATITUDE, color=as.factor(cluster)), size=1.5, alpha=0.8)+
  scale_colour_manual(values=cbp2, name = "Cluster")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title=element_text(size=22),
        legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 18))+
  xlab("Longitude")+
  ylab("Latitude")
