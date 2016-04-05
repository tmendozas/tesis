library(ggplot2)

wd <- "/Users/tania/tesis-itam/tesis/"

### Karate Club ###
karate_data <- read.table(paste(wd,"karate/resultskarate-multi.csv",sep=""),header = TRUE,sep = ",")
karate_data

ggplot(data=karate_data, aes(x=MarkovTime, y=meanVI)) +
  geom_line(colour="#3399FF") +
  geom_point(colour="#004C99") +
  ggtitle("Karate Club")


ggplot(data=karate_data, aes(x=MarkovTime, y=meanVI)) +
  geom_line(colour="#3399FF") +
  geom_point(colour="#004C99") +
  ggtitle("Karate Club") + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

ggplot(data=karate_data, aes(x=MarkovTime, y=numClusters)) +
  geom_point(colour="#004C99") +
  ggtitle("Karate Club")

### 721x280

### Les Miserables ###
lesmis_data <- read.table(paste(wd,"lesmis/results/first-result.txt",sep=""),header = TRUE,sep = "\t")
lesmis_data

ggplot(data=lesmis_data, aes(x=MarkovTime, y=meanVI)) +
  geom_line(colour="#3399FF") +
  geom_point(colour="#004C99") +
  ggtitle("Les Miserables")

categs_data <- lesmis_data

categs_data$ComsCount <- ifelse(categs_data$numClusters>12 ,'13 to 59', categs_data$numClusters)
categs_data$ClusterCount <- factor(categs_data$ComsCount, levels = c("13 to 59","12","10","9","8","7","6","5","4","3","2"))
categs_data$ID<-seq.int(nrow(categs_data))
categs_data

require(scales)
n <- length(levels(categs_data$ClusterCount))
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))]

ggplot(data=categs_data, aes(x=MarkovTime, y=meanVI,colour=ClusterCount)) +
  geom_line(colour="#D3D3D3") +
  geom_point(size=3) +
  scale_color_manual(values = cols)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position="left",legend.title=element_blank())



 # number of colors
 # color palette in random order
ggplot(data,aes(x,y,colour=category))+stat_smooth(alpha=0, size=2) + 
  scale_color_manual(values = cols) 
