###########################################################
###########################################################
#habitat amount hypothesis tests:
#buffer area predicts pop size/occupancy
#the configuration/isolation within this area doesn't matter
###########################################################

#######
#calculate: 
#1) patch size
#2) buffer area (with and without patch size)
#3) patch isolation (she focuses on nnd)
#4) landscape isolation

#unequal sampling effort on patches:
#for abundance, the B for patch size should be >1 if a bigger effect than the sampling issue alone

#predictions:
#HA: for the same habitat amount, no effect of patch size (positive effect of patch size: island effect)

library(reshape2)
library(plyr)
library(ggplot2)
library(survival) 
library(MuMIn)
library(lme4)
library(glmmADMB)#zero-inflated mixed model

#working directory
load("Data/habitat_amount_2014_2016.RData")

################################################################
################################################################
#survey data
################################################################
################################################################

#plotsummary=read.csv("C:\\Research\\Cactus bugs\\EAGER\\Analysis\\plot_trt_summary.csv",head=T)
#surveys=read.csv("C:\\Research\\Cactus bugs\\EAGER\\Analysis\\survey_dates_2014_2016.csv",head=T) 
#community=read.csv("C:/Research/Cactus bugs/EAGER/Analysis/community/Community_2014_2016.csv",head=T) 
#patch.coordinates.2015<-read.csv("C:/Research/Cactus bugs/EAGER/Analysis/cactus/cactus_quality.csv",head=T) 
#patch.coordinates.2014<-read.csv("C:/Research/Cactus bugs/EAGER/Analysis/patch_2014.csv",header=T)

#patch.coordinates.2015$posttrt<-"Y"
#coordinates<-merge(patch.coordinates.2014, patch.coordinates.2015, by.x=c("plot","patch"),by.y=c("Plot","Patch"), all=T)

#surveys<-subset(surveys, select=c("Plot", "Survey.Date","SurveyNumber"))
#community<-merge(community,plotsummary, by="Plot",all=TRUE)
#community<-merge(community,surveys, by=c("Plot","Survey.Date"), all=FALSE)

#format date
#community$date<-as.POSIXct(strptime(community$Survey.Date, "%m/%d/%Y"))
#community$year<-format(community$date,'%Y')
#community$day<-format(community$date,'%d')

#two plots (1,12) in 2014 not surveyed in 2015: remove those
#community$Plot<-factor(community$Plot)
#community<-community[! community$Plot %in% c('1', '12'), ]


################################################################
################################################################
#landscape metrics
################################################################
################################################################

#make new area metric
coordinates$AreaMean2<-coordinates$area
coordinates$AreaMean2<-ifelse(is.na(coordinates$AreaMean2),coordinates$AreaMean,coordinates$AreaMean2)

#----------------------------------#
#distance matrices
#----------------------------------#

#make lists for each plot for distance/area
dist.matrix.2015=list()#from patch centers
distedge.matrix.2015=list()#from patch edges (approximate)
distfaredge.matrix.2015=list() #center to patch far edge distance matrix
area.2015=list()
xy.2015=list()
plotID.2015=list()
patchID.2015=list()

#2015
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

#2014
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


#----------------------------------#
#patch metrics
#----------------------------------#

meandist=12.5 #2014 monthly data, removing first check (transient from release)
alpha<-1/meandist

#patch metrics
#currently: code only considers buffer for buffer area
#metrics could be changed to only consider patches within buffer area

ifm2015=list()
bufferarea2015=list()
nnd2015=list()
pc2015=list()


ifm2014=list()
bufferarea2014=list()
nnd2014=list()
pc2014=list()


#2015
library(igraph)

for(i in 2:17){  
  if (i==12) next
  
  #area/dist
  dist.i<-dist.matrix.2015[[i]]
  area.i<-area.2015[[i]]
  
  #distance metrics
  dist.inv<-1/dist.i
  diag(dist.inv)<-0
  dist.nnd.i<-dist.i
  diag(dist.nnd.i)<-NA
  nnd2015[[i]]<-apply(dist.nnd.i,1,min,na.rm=T)
  
  #ifm, pc
  g.i = exp(-alpha*dist.i)
  diag(g.i)<-0
  ifm2015[[i]]<-rowSums(sweep(g.i, 2, area.i, "*"))#excludes incidence
  pc2015[[i]]<-prob.connectivity(area=area.i,prob.matrix=g.i, landarea=2500)
  
  #buffer area
  dist.alpha<-ifelse(dist.i>meandist,0,1)#movement dist
  diag(dist.alpha)<-0 #doesn't include focal patch
  bufferarea2015[[i]]<-rowSums(sweep(dist.alpha, 2, area.i, "*"))  
}

#2014
for(i in 2:17){  
  if (i==12) next
  
  #area/dist
  dist.i<-dist.matrix.2014[[i]]
  area.i<-area.2014[[i]]
  
  #distance metrics
  dist.inv<-1/dist.i
  diag(dist.inv)<-0
  dist.nnd.i<-dist.i
  diag(dist.nnd.i)<-NA
  nnd2014[[i]]<-apply(dist.nnd.i,1,min,na.rm=T)
  
  #ifm, pc
  g.i = exp(-alpha*dist.i)
  diag(g.i)<-0
  ifm2014[[i]]<-rowSums(sweep(g.i, 2, area.i, "*"))#excludes incidence
  pc2014[[i]]<-prob.connectivity(area=area.i,prob.matrix=g.i, landarea=2500)
  
  #buffer area
  dist.alpha<-ifelse(dist.i>meandist,0,1)#movement dist
  diag(dist.alpha)<-0 #doesn't include focal patch
  bufferarea2014[[i]]<-rowSums(sweep(dist.alpha, 2, area.i, "*"))  
}


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

#get total population size in bufferarea
bufferN=vector(length=nrow(community))

#-------------------------------------------#
#link metrics with survey data
#-------------------------------------------#
#get mean by plot
lapply(nnd2015,mean,na.rm=T)

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
landscape.df<-rbind(landscape2015.df, landscape2014.df)

head(community,1)
head(landscape.df,1)

community<-merge(community, landscape.df, by=c("Plot","Patch"), all=F)
names(community)[names(community) == "patcharea.x"] <- "patcharea"
names(community)[names(community) == "nnd.x"] <- "nnd"
names(community)[names(community) == "bufferarea.x"] <- "bufferarea"
names(community)[names(community) == "ifm.x"] <- "ifm"
head(community,5)

#combine patch size with surrounding bufferarea for total area
#community$tbufferarea<-community$patch + community$bufferarea
community$tbufferfararea<-community$patcharea + community$bufferfararea

################################################################
################################################################
#statistical models
################################################################
################################################################

str(community)
community2015<-subset(community, SurveyNumber>5 & SurveyNumber<23)
community2015.nocont<-subset(community2015, TrtType!="cont")


#ZIP 2015 only
zipm.CheliAdN.int<-glmmadmb(Cheli.adults~(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.patchsize<-glmmadmb(Cheli.adults~scale(patcharea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.bufferarea<-glmmadmb(Cheli.adults~scale(bufferarea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.tbufferarea<-glmmadmb(Cheli.adults~scale(tbufferarea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.TypeXLoss<-glmmadmb(Cheli.adults~TrtType*PropLost+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)


zipm.CheliAdN.bufferfararea<-glmmadmb(Cheli.adults~scale(bufferfararea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.tbufferfararea<-glmmadmb(Cheli.adults~scale(tbufferfararea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)


zipm.CheliAdN.patchsizeXbufferarea<-glmmadmb(Cheli.adults~scale(patcharea)*scale(bufferarea)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.patchsizeXType<-glmmadmb(Cheli.adults~scale(patcharea)*TrtType+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.patchsizeXLoss<-glmmadmb(Cheli.adults~scale(patcharea)*PropLost+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.tbufferareaXType<-glmmadmb(Cheli.adults~scale(tbufferarea)*TrtType+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)
zipm.CheliAdN.tbufferareaXLoss<-glmmadmb(Cheli.adults~scale(tbufferarea)*PropLost+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)

zipm.CheliAdN.nndXtbufferarea<-glmmadmb(Cheli.adults~scale(tbufferarea)*scale(nnd)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)

model.sel(zipm.CheliAdN.int,zipm.CheliAdN.patchsize,zipm.CheliAdN.bufferarea,
          zipm.CheliAdN.tbufferarea,
          zipm.CheliAdN.nndXtbufferarea,
          zipm.CheliAdN.TypeXLoss,
          zipm.CheliAdN.patchsizeXbufferarea,
          zipm.CheliAdN.patchsizeXType,zipm.CheliAdN.patchsizeXLoss,
          zipm.CheliAdN.tbufferareaXType,zipm.CheliAdN.tbufferareaXLoss)

summary(zipm.CheliAdN.tbufferareaXType)
summary(zipm.CheliAdN.patchsizeXType)
summary(zipm.CheliAdN.nndXtbufferarea)

#ZIP pschy
library(pscl)

zip.CheliAdN.patchsize<-zeroinfl(Cheli.adults ~ patcharea | 1, data = community2015.nocont)
zip.CheliAdN.bufferarea<-zeroinfl(Cheli.adults ~ bufferarea | 1, data = community2015.nocont)
zip.CheliAdN.tbufferarea<-zeroinfl(Cheli.adults ~ tbufferarea | 1, data = community2015.nocont)
zip.CheliAdN.patchsizeXTypeXLoss<-zeroinfl(Cheli.adults ~ scale(patcharea)*TrtType*PropLost| 1, data = community2015.nocont)
zip.CheliAdN.patchsizeXTypeXLossXbuffer<-zeroinfl(Cheli.adults ~ scale(patcharea)*bufferarea*TrtType*PropLost| 1, data = community2015.nocont)

zipm.CheliAdN.nndXtbufferarea<-glmmadmb(Cheli.adults~scale(tbufferarea)*scale(nnd)+(1|Plot), data=community2015.nocont, family="Poisson", zeroInflation=T)


zip.CheliAdN.bufferXndd<-zeroinfl(Cheli.adults ~ tbufferarea*nnd | 1, data = community2015.nocont)
zip.CheliAdN.bufferXpatcharea<-zeroinfl(Cheli.adults ~ bufferarea*patcharea | 1, data = community2015.nocont)

zip.CheliAdN.bufferfararea<-zeroinfl(Cheli.adults ~ bufferfararea | 1, data = community2015.nocont)
zip.CheliAdN.tbufferfararea<-zeroinfl(Cheli.adults ~ tbufferfararea | 1, data = community2015.nocont)

AIC(zip.CheliAdN.patchsize,zip.CheliAdN.bufferarea,zip.CheliAdN.tbufferarea, zip.CheliAdN.tbufferfararea, zip.CheliAdN.bufferfararea)
summary(zip.CheliAdN.patchsize)
summary(zip.CheliAdN.bufferarea)
summary(zip.CheliAdN.tbufferarea)
summary(zip.CheliAdN.tbufferfararea)


summary(zip.CheliAdN.patchsizeXTypeXLossXbuffer)
summary(zip.CheliAdN.bufferXpatcharea)
summary(zip.CheliAdN.bufferarea)
summary(zip.CheliAdN.patchsize)

#save.image("C:/Research/Cactus bugs/EAGER/Analysis/Habitat amount test/habitat_amount_2014_2016.RData")

#########################################################################
#########################################################################
#functions
#########################################################################
#########################################################################

prob.connectivity <- function(prob.matrix,area,landarea){
  
  #use log-transformed weights and get shortest path
  pc.graph<-graph.adjacency(prob.matrix, mode="undirected", weighted=TRUE)
  pstar.mat<-shortest.paths(pc.graph, weights=-log(E(pc.graph)$weight))
  pstar.mat<-exp(-pstar.mat)#back-transform to probabilities
  
  PCmat<-outer(area, area)*pstar.mat
  PC<-sum(PCmat)/landarea^2
  
  N=nrow(prob.matrix)
  dPC<-rep(NA,N)
  
  for (i in 1:N) {
    prob.matrix.i=prob.matrix[-i,-i] 
    area.i<-area[-i]
    
    pc.graph.i<-graph.adjacency(prob.matrix.i, mode="undirected", weighted=TRUE)
    pstar.mat.i<-shortest.paths(pc.graph.i, weights=-log(E(pc.graph.i)$weight))
    pstar.mat.i<-exp(-pstar.mat.i)
    
    PCmat.i<-outer(area.i, area.i)*pstar.mat.i
    PC.i<-sum(PCmat.i)/landarea^2
    
    dPC[i] <- (PC-PC.i)/PC*100
  }
  
  dPCintra<-area^2/sum(PCmat)*100 
  dPCflux<- 2*(rowSums(PCmat)-area^2)/sum(PCmat)*100
  dPCconn<-dPC-dPCintra-dPCflux
  
  patch.metrics=data.frame(dPC,dPCintra,dPCflux,dPCconn)
  pc.list = list(PC=PC,patch=patch.metrics)
  return(pc.list)   
}

