plot.2<-subset(coordinates, plot==2 & posttrt=="Y")
plot.2<-plot.2[complete.cases(plot.2[,"East"]),]
dist.matrix.2<-as.matrix(dist(cbind(plot.2$East,plot.2$North))) #Euclidean distance matrix
area.2<-plot.2$AreaMean2

radius.2<-(plot.2$AreaMean/pi)^0.5      ###Radius of each patch
dist.edge.cent.2 <- sweep(dist.matrix.2, 2, radius.2, "-")  ###Sweeps across columns subracting radius from distance matix

dist.alpha.2<-ifelse(dist.edge.cent.2>meandist,0,1)#movement dist
diag(dist.alpha.2)<-0 #doesn't include focal patch



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



######Prop area



##2014
for(i in 1:17){  
  #if (i==12) next
  
  #area/dist
  
  dist.2 <-distfaredge.matrix.2014[[12]]
  area.2<-area.2014[[12]]
  
  #distance metrics
  dist.inv2<-1/dist.2
  diag(dist.inv2)<-0
  
  #buffer area
  dist.alpha2<-ifelse(dist.2>meandist,0,1)#movement dist
  diag(dist.alpha2)<-0 #doesn't include focal patch
  bufferfararea2014.2<-rowSums(sweep(dist.alpha2, 2, area.2, "*"))  
}



