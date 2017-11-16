#####################Calculation of subeseted plot to avoid edge effects
####Calculations


###Calulate patches within boundary for each plot
####Make column in coordinates that indicates wether its winthin new boundarys
###Subset community dataframe to include only patches in the subplot area
###Rerun ZIP models using this subsetted data

pc_minmax <- read.csv("Data/plot_corners_minmax.csv")

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

subsetmerge<- coordinates[,c(1,2,15)]
names(subsetmerge)[1:2] <- c("Plot","Patch")
  
communityNew <-  merge(community2015.nocont,subsetmerge , by=c("Plot","Patch"), all=F)
communitySub <- subset(communityNew, coordInc == 1)


######################

library(pscl)

zip.CheliAdN.patchsize2<-zeroinfl(Cheli.adults ~ patcharea | 1, data = communitySub)
zip.CheliAdN.bufferarea2<-zeroinfl(Cheli.adults ~ bufferarea | 1, data = communitySub)
zip.CheliAdN.tbufferarea2<-zeroinfl(Cheli.adults ~ tbufferarea | 1, data = communitySub)
zip.CheliAdN.patchsizeXTypeXLoss<-zeroinfl(Cheli.adults ~ scale(patcharea)*TrtType*PropLost| 1, data = communitySub)
zip.CheliAdN.patchsizeTypeLoss<-zeroinfl(Cheli.adults ~ patcharea+TrtType+PropLost| 1, data = communitySub)

zip.CheliAdN.bufferfararea<-zeroinfl(Cheli.adults ~ bufferfararea | 1, data = communitySub)
zip.CheliAdN.tbufferfararea<-zeroinfl(Cheli.adults ~ tbufferfararea | 1, data = communitySub)


AIC(zip.CheliAdN.patchsize,zip.CheliAdN.bufferarea,zip.CheliAdN.tbufferarea, zip.CheliAdN.tbufferfararea, zip.CheliAdN.bufferfararea)

summary(zip.CheliAdN.patchsize)
summary(zip.CheliAdN.patchsize2)


summary(zip.CheliAdN.bufferarea)
summary(zip.CheliAdN.bufferarea2)

summary(zip.CheliAdN.tbufferarea)
summary(zip.CheliAdN.tbufferarea2)