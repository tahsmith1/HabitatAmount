###############Calculates buffer at different radii############
community <- read.csv("ManData/community.csv")
coordinates <-read.csv("ManData/coordinates.csv")

######

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
#####


landscape2015.df<-data.frame(Plot=unlist(plotID.2015),Patch=unlist(patchID.2015),
                             patcharea=unlist(area.2015))
landscape2014.df<-data.frame(Plot=unlist(plotID.2014),Patch=unlist(patchID.2014),
                             patcharea=unlist(area.2014))


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


############ subset by Aug-Nov

comm2014<- subset(community, subset = community$year == 2014)
comm2015<- subset(community, subset = community$year == 2015)

comm2014an <- subset(comm2014, comm2014$month >7)
comm2015an <- subset(comm2015, comm2015$month >7)


####Takes the 1st, 3rd, and 5th observation in 2015 to balance it with the 3 monthly obs in 2014
comhold = list()
chold = list()
library(data.table)

for(i in 2:17){  
  if (i==12) next
  comm.i <- subset(comm2015an, comm2015an$Plot == i)
  chold = list()
  for(ii in  1:length(unique(comm.i$Patch))){
    pcomm.ii <- subset(comm.i, comm.i$Patch ==  unique(comm.i$Patch)[ii])
    pcomm.ii <- pcomm.ii[order(pcomm.ii$SurveyNumber),]
    chold[[ii]]<- pcomm.ii[c(1,3,5),]
  }
  comhold[[i]] <- rbindlist(chold)
  
}

comm2015an3<- rbindlist(comhold)

################# Removing patches in 2014 that don't remain in 2015 after treatment
#communityYear<- rbind(comm2014an,comm2015an3)

leftpatches <- comm2015an3[,c("Patch","Plot")]
leftpatches<- unique(leftpatches)
comm2014ms <- merge(comm2014an, leftpatches, by = c("Plot", "Patch"), all= F)

leftpatches1<- comm2014ms[,c("Patch","Plot")]
leftpatches1 <-unique(leftpatches1)
comm2015ms<- merge(comm2015an3, leftpatches1, by = c("Plot", "Patch"), all= F)

###Make combined commmunity data frame from 2014 and 2015 data

communityYear<- rbind(comm2014ms,comm2015ms)

########Plot, patch and year made into factors for analysis

communityYear$Plot<- as.factor(communityYear$Plot)
communityYear$Patch<- as.factor(communityYear$Patch)
communityYear$year<- as.factor(communityYear$year)


###Creates data frame without plots in which no habitat was removed (3 plots not used)
comm.nocont<- subset(communityYear, TrtType != "cont")


###############################################  Analysis

col_buffnames <- names(comm.nocont)[40:88]

zipmPB = list()

for(ii in 40:88){
  zipmPB[[ii-39]] <-glmmadmb(Cheli.adults~scale(patcharea.x)+scale(comm.nocont[,ii])+(1|Plot), data=comm.nocont, family="Poisson", zeroInflation=T)
}

names(zipmPB)<- col_buffnames


