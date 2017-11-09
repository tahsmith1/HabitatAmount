####Loop deconstruction using plot 2 in 2015


dist.nnd<- dist.matrix.2015[[2]]    ##Nearest neighbor distance
diag(dist.nnd)<-NA                  ##removes 0s from diagnol so min works in next line
nnd2015t<-apply(dist.nnd,1,min,na.rm=T) ##takes the minimum distance from each row in the matrix i.e. the shorest distance between patched for each patch
 



gi = exp(-alpha*dist.matrix.2015[[2]])  ##probabity based on distance and mean travel
diag(gi)<-0                             

a2015<- area.2015[[2]]

ifm2 <-rowSums(sweep(gi, 2, a2015, "*"))              ###ifm? area times the distance and mean travel

pc2<-prob.connectivity(area=a2015,prob.matrix=gi, landarea=2500) ##using connectivity function


dist.al<-ifelse(dist.nnd>meandist,0,1)#movement dist    ### codes 1 for all patches within the mean dispersal dist
diag(dist.al)<-0 #doesn't include focal patch
bufferarea2<-rowSums(sweep(dist.al, 2, a2015, "*"))  #Sums up the area of all these patches


##################
