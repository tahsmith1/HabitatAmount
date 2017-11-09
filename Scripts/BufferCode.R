plot.2<-subset(coordinates, plot==2 & posttrt=="Y")
plot.2<-plot.2[complete.cases(plot.2[,"East"]),]
dist.matrix.2<-as.matrix(dist(cbind(plot.2$East,plot.2$North))) #Euclidean distance matrix
area.2<-plot.2$AreaMean2

radius.2<-(plot.2$AreaMean/pi)^0.5      ###Radius of each patch
dist.edge.cent.2 <- sweep(dist.matrix.2, 2, radius.2, "-")  ###Sweeps across columns subracting radius from distance matix

dist.alpha.2<-ifelse(dist.edge.cent.2>meandist,0,1)#movement dist
diag(dist.alpha.2)<-0 #doesn't include focal patch


######Prop area

