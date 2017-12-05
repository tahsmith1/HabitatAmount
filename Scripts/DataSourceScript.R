#################Script for calculations###########################

library(reshape2)
library(plyr)
library(ggplot2)
library(survival) 
library(MuMIn)
library(lme4)
library(glmmADMB)#zero-inflated mixed model


################

load("Data/habitat_amount_2014_2016.RData")

community<- unique(community)

#make lists for each plot for distance/area
dist.matrix.2015=list()#from patch centers
distedge.matrix.2015=list()#from patch edges (approximate)
distfaredge.matrix.2015=list() #center to patch far edge distance matrix
area.2015=list()
xy.2015=list()
plotID.2015=list()
patchID.2015=list()

#2015 distance matrix calculations
for(i in 1:17){  
  plot.i<-subset(coordinates, plot==i & posttrt=="Y")
  #remove NA locations posttrt
  plot.i<-plot.i[complete.cases(plot.i[,"East"]),]
  #plot.i<-subset(plot.i, dead.!="Y")#remove dead patches in 2015 (dead by fall; area 2014 measurements at end of season)
  patchID.2015[[i]]<-plot.i$patch
  plotID.2015[[i]]<-plot.i$plot
  
  #calculate/store dist/area
  dist.matrix.2015[[i]]<-as.matrix(dist(cbind(plot.i$East,plot.i$North))) #Euclidean distance matrix
  area.2015[[i]]<-plot.i$AreaMean2
  
  #approximate edge to edge distance
  radius.i<-(plot.i$AreaMean/pi)^0.5
  distedge.matrix.2015[[i]]<-dist.matrix.2015[[i]]-outer(radius.i,radius.i,"+")
  xy.2015[[i]]<-cbind(plot.i$East,plot.i$North)
  
  #Distance matrices to far edge
  distfaredge.matrix.2015[[i]] <- sweep(dist.matrix.2015[[i]], 2, radius.i, "+")
  diag(distfaredge.matrix.2015[[i]])<-0 #doesn't include focal patch
}

#2014 distance matrix calculations
dist.matrix.2014=list()
distfaredge.matrix.2014=list() #center to patch far edge distance matrix
area.2014=list()
xy.2014=list()
plotID.2014=list()
patchID.2014=list()

for(i in 2:17){  
  if (i==12) next
  plot.i<-subset(coordinates, plot==i)
  #remove NA locations
  plot.i<-plot.i[complete.cases(plot.i[,c("area","East")]),]
  patchID.2014[[i]]<-plot.i$patch
  plotID.2014[[i]]<-plot.i$plot
  
  #calculate/store dist/area
  dist.matrix.2014[[i]]<-as.matrix(dist(cbind(plot.i$East,plot.i$North))) #Euclidean distance matrix
  area.2014[[i]]<-plot.i$AreaMean2
  xy.2014[[i]]<-cbind(plot.i$East,plot.i$North)
  
  #Distance matrices to far edge
  radius.i<-(plot.i$AreaMean2/pi)^0.5
  distfaredge.matrix.2014[[i]] <- sweep(dist.matrix.2014[[i]], 2, radius.i, "+")
  diag(distfaredge.matrix.2014[[i]])<-0 #doesn't include focal patch
}

##############################

###Calculate the buffer area based on the far edges of patches (most conservative way to calculate)
##2015
bufferfararea2015=list()
for(i in 2:17){  
  if (i==12) next
  
  #area/dist
  
  dist.i<-distfaredge.matrix.2015[[i]]
  area.i<-area.2015[[i]]
  
  #distance metrics
  dist.inv<-1/dist.i
  diag(dist.inv)<-0
  
  #buffer area
  dist.alpha<-ifelse(dist.i>meandist,0,1)#movement dist
  diag(dist.alpha)<-0 #doesn't include focal patch
  bufferfararea2015[[i]]<-rowSums(sweep(dist.alpha, 2, area.i, "*"))  
}

bufferfararea2014=list()
##2014
for(i in 2:17){  
  if (i==12) next
  
  #area/dist
  
  dist.i<-distfaredge.matrix.2014[[i]]
  area.i<-area.2014[[i]]
  
  #distance metrics
  dist.inv<-1/dist.i
  diag(dist.inv)<-0
  
  #buffer area
  dist.alpha<-ifelse(dist.i>meandist,0,1)#movement dist
  diag(dist.alpha)<-0 #doesn't include focal patch
  bufferfararea2014[[i]]<-rowSums(sweep(dist.alpha, 2, area.i, "*"))  
}


#-------------------------------------------#
#link metrics with survey data
#-------------------------------------------#

landscape2015.df<-data.frame(Plot=unlist(plotID.2015),Patch=unlist(patchID.2015),
                             patcharea=unlist(area.2015),
                             nnd=unlist(nnd2015), 
                             bufferarea=unlist(bufferarea2015),
                             bufferfararea=unlist(bufferfararea2015),
                             ifm=unlist(ifm2015))
landscape2014.df<-data.frame(Plot=unlist(plotID.2014),Patch=unlist(patchID.2014),
                             patcharea=unlist(area.2014),
                             nnd=unlist(nnd2014), 
                             bufferarea=unlist(bufferarea2014), 
                             bufferfararea=unlist(bufferfararea2014),
                             ifm=unlist(ifm2014))

landscape2014.df$year <- 2014
landscape2015.df$year <- 2015

#write.csv(landscape2014.df,"ManData/landscape2014.csv")
#write.csv(landscape2015.df,"ManData/landscape2015.csv")

landscape.df<-rbind(landscape2015.df, landscape2014.df)

community<- community[,c(-33,-32, -31, -30, -29)]
community <- unique(community)

communityt<-merge(community, landscape.df, by=c("Plot","Patch","year"), all=F)


communityt$tbufferfararea<-communityt$patcharea + communityt$bufferfararea
communityt$tbufferarea<-communityt$patcharea + communityt$bufferarea

community <- communityt

##### Add month column
community$Survey.Date <- as.Date(community$Survey.Date, format = "%m/%d/%Y")
community$month <- months(community$Survey.Date)
community$month <-factor(community$month, levels= c("January","February","March","April", "May", "June", "July", "August","September", "October", "November", "December"))
community$month<- as.numeric(community$month)

#####################Calculation of subseseted plot to avoid edge effects
####Calculations

meandist=12.5 #2014 monthly data, removing first check (transient from release)
alpha<-1/meandist

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
################   Merges subsetted coordinate data with community date toi include only patchs > 12.5m from edge of plots
coordinates$coordInc <- unlist(subcoordInd)

subsetmerge<- coordinates[,c(2,3,16)]
names(subsetmerge)[1:2] <- c("Plot","Patch")

communityNew <-  merge(community,subsetmerge , by=c("Plot","Patch"), all = F)
communitySub <- subset(communityNew, coordinates$coordInc == 1)



####writes already calculated community and coordinate data frames#######


write.csv(community,"ManData/community.csv")
write.csv(communitySub, "ManData/communitySub.csv")

community2015<-subset(community, SurveyNumber>5 & SurveyNumber<23)
community2015.nocont<-subset(community2015, TrtType!="cont")

write.csv(community2015,"ManData/community2015.csv")
write.csv(community2015.nocont, "ManData/community2015.nocont.csv")

#################################################################################
