#####################Calculation of subeseted plot to avoid edge effects
####Calculations


###Calulate patches within boundary for each plot
####Make column in coordinates that indicates wether its winthin new boundarys
###Subset community dataframe to include only patches in the subplot area
###Rerun ZIP models using this subsetted data

coordinates<-  read.csv("ManData/coordinates.csv")
pc_minmax <- read.csv("Data/plot_corners_minmax.csv")

meandist=12.5 #2014 monthly data, removing first check (transient from release)

pc_minmax$subXmin <- pc_minmax$Xmin + meandist
pc_minmax$subXmax <- pc_minmax$Xmax - meandist
pc_minmax$subYmin <- pc_minmax$Ymin + meandist
pc_minmax$subYmax <- pc_minmax$Ymax - meandist

inmax <- subset(pc_minmax, Plot==1)

coordinates$coordInc <- NA
head(coordinates)

subcoordInd = list()

for(i in 1:17){  
  plot.i<-subset(coordinates, plot==i)
  c_minmax <- subset(pc_minmax, Plot==i)
  subcoordInd[[i]]<- ifelse(plot.i$East > c_minmax$subXmin & plot.i$East < c_minmax$subXmax & 
           plot.i$North > c_minmax$subYmin & plot.i$North < c_minmax$subYmax,1,0)
  
}
################
coordinates$coordInc <- unlist(subcoordInd)

subsetmerge<- coordinates[,c(2,3,16)]
names(subsetmerge)[1:2] <- c("Plot","Patch")

#write.csv(subsetmerge,"ManData/subsetmerge.csv")
  
#communityNew <-  merge(community2015.nocont,subsetmerge , by=c("Plot","Patch"), all=F)
#communitySub <- subset(communityNew, coordInc == 1)



