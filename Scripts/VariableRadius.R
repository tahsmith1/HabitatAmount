###############Calculates buffer at different radii############
community <- read.csv("ManData/community.csv")
coordinates <-read.csv("ManData/coordinates.csv")

######

landscape2015.df<-data.frame(Plot=unlist(plotID.2015),Patch=unlist(patchID.2015),
                             patcharea=unlist(area.2015))
landscape2014.df<-data.frame(Plot=unlist(plotID.2014),Patch=unlist(patchID.2014),
                             patcharea=unlist(area.2014))


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
meandist=12.5

bufferfar_calc15 <- function(meandist) {
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
  return(bfararea2015=unlist(bufferfararea2015))
}  
bufferfar_calc14 <- function(meandist) {  
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
  return(bfararea2014=unlist(bufferfararea2014))
}

################################Runs function to calculate for all 49 in seq

bufferradii <- seq(1,25,by=0.5)
bufferfar2015<- lapply(bufferradii,FUN=bufferfar_calc15)
bufferfar2014<- lapply(bufferradii,FUN=bufferfar_calc14)

####Saves as dataframe with columns for each of the 49 calcs
bf2015 <- as.data.frame(bufferfar2015)
names(bf2015) <- paste("bfararea", bufferradii, sep= "_")

bf2014 <- as.data.frame(bufferfar2014)
names(bf2014) <- paste("bfararea", bufferradii, sep= "_")

lands_2014<- cbind(landscape2014.df, bf2014, year= 2014)
lands_2015<- cbind(landscape2015.df, bf2015, year = 2015)

lland<- rbind(lands_2014,lands_2015)

communityt<-merge(community, lland, by=c("Plot","Patch","year"), all=F)

community <- communityt

##### Add month column
community$Survey.Date <- as.Date(community$Survey.Date, format = "%m/%d/%Y")
community$month <- months(community$Survey.Date)
community$month <-factor(community$month, levels= c("January","February","March","April", "May", "June", "July", "August","September", "October", "November", "December"))
community$month<- as.numeric(community$month)


       