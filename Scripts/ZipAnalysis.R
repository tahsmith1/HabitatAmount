####Reads in already calculated community and coordinate data frames#######


community<- read.csv("ManData/community.csv")
#community2015 <- read.csv('ManData/community2015.csv')
#community2015.nocont <- read.csv("ManData/community2015.nocont.csv")
#communitySub <- read.csv("ManData/communitySub.csv")
landscape2014 <- read.csv("ManData/landscape2014.csv")
landscape2015 <- read.csv("ManData/landscape2015.csv")

#################

treatprop<- unique(community[,c("Plot","Patch","TrtType","PropLost")])

############ subset by Aug-Nov

comm2014<- subset(community, subset = community$year == 2014)
comm2015<- subset(community, subset = community$year == 2015)

comm2014an <- subset(comm2014, comm2014$month >7)
comm2015an <- subset(comm2015, comm2015$month >7)

#commlate <- rbind(comm2014an,comm2015an)

##########
#str(treatprop)
lnd2014<- landscape2014[,-1]

landsc2014<- merge(lnd2014,treatprop, by= c("Plot","Patch"), all = F)

#mergedcomm2014 <- aggregate(Cheli.adults~Plot+Patch, data= comm2014an, FUN= sum)
#comm2014m<- merge(mergedcomm2014,landsc2014, by = c("Plot", "Patch"), all = F)

lnd2015<- landscape2015[,-1]
landsc2015<- merge(lnd2015,treatprop, by= c("Plot","Patch"), all = F)

#mergedcomm2015 <- aggregate(Cheli.adults~Plot+Patch, data= comm2015an, FUN= sum)
#comm2015m<- merge(mergedcomm2015,landsc2015, by = c("Plot", "Patch"), all = F)

####an3 in subsetobervation.R
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
#################

#communityYear<- rbind(comm2014an,comm2015an3)

leftpatches <- comm2015an3[,c("Patch","Plot")]
leftpatches<- unique(leftpatches)
comm2014ms <- merge(comm2014an, leftpatches, by = c("Plot", "Patch"), all= F)

communityYear<- rbind(comm2014ms,comm2015an3)

communityYear$Plot<- as.factor(communityYear$Plot)
communityYear$Patch<- as.factor(communityYear$Patch)
communityYear$year<- as.factor(communityYear$year)

comm.nocont<- subset(communityYear, TrtType != "cont")
#########
subsetmerge <- read.csv("ManData/subsetmerge.csv")
communityNew <-  merge(comm.nocont,subsetmerge , by=c("Plot","Patch"), all=F)
comm.sub <- subset(communityNew, coordInc == 1)

#######

library(pscl)
library(glmmADMB)

#############

cor(comm.nocont$bufferarea,comm.nocont$bufferfararea)

cor(comm.nocont$bufferarea,comm.nocont$patcharea)
cor(comm.nocont$bufferfararea,comm.nocont$patcharea)

cor(comm.nocont$bufferarea,comm.nocont$patcharea+comm.nocont$bufferarea)
####################Testing patch area and nnd effects against bufferarea across treatments
zPatch<-zeroinfl(Cheli.adults ~ patcharea*year*TrtType | 1, data = comm.nocont)
zBuffer<-zeroinfl(Cheli.adults ~ bufferarea*year*TrtType | 1, data = comm.nocont)
zPatchBuffer<-zeroinfl(Cheli.adults ~ patcharea*bufferarea*year*TrtType | 1, data = comm.nocont)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferarea*nnd*year*TrtType | 1, data = comm.nocont)
zNND<-zeroinfl(Cheli.adults ~ nnd*year*TrtType | 1, data = comm.nocont)
zNNDPatch<-zeroinfl(Cheli.adults ~ nnd*patcharea*year*TrtType | 1, data = comm.nocont)
#zNNDPatchBuffer <-zeroinfl(Cheli.adults ~ nnd*patcharea*bufferarea*year*TrtType | 1, data = comm.nocont)  NA groups

summary(zPatch)
summary(zBuffer)
summary(zNND)
summary(zPatchBuffer)
summary(zNNDBuffer)

AIC(zPatch,zBuffer,zNND,zPatchBuffer,zNNDBuffer, zNNDPatch,zNNDPatchBuffer)


###########################Farbuffer area differences

zPatch<-zeroinfl(Cheli.adults ~ patcharea*year*TrtType | 1, data = comm.nocont)
zBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*year*TrtType | 1, data = comm.nocont)
zPatchBuffer<-zeroinfl(Cheli.adults ~ patcharea*bufferfararea*year*TrtType | 1, data = comm.nocont)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*nnd*year*TrtType | 1, data = comm.nocont)
zNND<-zeroinfl(Cheli.adults ~ nnd*year*TrtType | 1, data = comm.nocont)
zNNDPatch<-zeroinfl(Cheli.adults ~ nnd*patcharea*year*TrtType | 1, data = comm.nocont)
#zNNDPatchBuffer <-zeroinfl(Cheli.adults ~ nnd*patcharea*bufferarea*year*TrtType | 1, data = comm.nocont)  NA groups

summary(zPatch)
summary(zBuffer)
summary(zNND)
summary(zPatchBuffer)
summary(zNNDBuffer)

AIC(zPatch,zBuffer,zNND,zPatchBuffer,zNNDBuffer, zNNDPatch)

#############################################

zPatch<-zeroinfl(Cheli.adults ~ patcharea*year*TrtType | 1, data = comm.sub)
zBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*year*TrtType | 1, data = comm.sub)
zPatchBuffer<-zeroinfl(Cheli.adults ~ patcharea*bufferfararea*year*TrtType | 1, data = comm.sub)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*nnd*year*TrtType | 1, data = comm.sub)
zNND<-zeroinfl(Cheli.adults ~ nnd*year*TrtType | 1, data = comm.sub)
zNNDPatch<-zeroinfl(Cheli.adults ~ nnd*patcharea*year*TrtType | 1, data = comm.sub)

AIC(zPatch,zBuffer,zNND,zPatchBuffer,zNNDBuffer, zNNDPatch)

summary(zPatch)
summary(zBuffer)
summary(zNND)
summary(zPatchBuffer)
summary(zNNDBuffer)

zmPatchBuffer<-glmmadmb(Cheli.adults~scale(bufferarea)*scale(patcharea)*TrtType*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
summary(zmPatchBuffer)


zPatch<-zeroinfl(Cheli.adults ~ patcharea*year*TrtType | 1, data = comm.sub)
zBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*year*TrtType | 1, data = comm.sub)
zPatchBuffer<-zeroinfl(Cheli.adults ~ patcharea*bufferfararea*year*TrtType | 1, data = comm.sub)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*nnd*year*TrtType | 1, data = comm.sub)
zNND<-zeroinfl(Cheli.adults ~ nnd*year*TrtType | 1, data = comm.sub)
zNNDPatch<-zeroinfl(Cheli.adults ~ nnd*patcharea*year*TrtType | 1, data = comm.sub)

zPatch<-zeroinfl(Cheli.adults ~ patcharea | 1, data = comm.sub)
zBuffer<-zeroinfl(Cheli.adults ~ bufferarea+patcharea | 1, data = comm.sub)

####################New Actual testing patch area and nnd effects against bufferarea across treatments
zPatch<-zeroinfl(Cheli.adults ~ patcharea | 1, data = comm.nocont)
zBuffer<-zeroinfl(Cheli.adults ~ patcharea + bufferarea | 1, data = comm.nocont)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferarea*nnd | 1, data = comm.nocont)
zNND<-zeroinfl(Cheli.adults ~ nnd | 1, data = comm.nocont)
zPropTrt<-zeroinfl(Cheli.adults ~PropLost*TrtType | 1, data = comm.nocont)
comm.nocont$PropLost

AICtable <- AIC(zPatch,zBuffer,zNNDBuffer,zNND,zPropTrt)
AICtable[order(AICtable$AIC),]


zmPatchBuffer<-glmmadmb(Cheli.adults~scale(bufferarea)+offset(scale(patcharea))+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmPatch<-glmmadmb(Cheli.adults~scale(patcharea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmNND<-glmmadmb(Cheli.adults~scale(nnd)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmNNDBuffer<-glmmadmb(Cheli.adults~scale(nnd)*scale(bufferarea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmPropTrt<-glmmadmb(Cheli.adults~PropLost*TrtType+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)

AICtable <- AIC(zmPatch,zmPatchBuffer,zmNND,zmNNDBuffer ,zmPropTrt)
AIC(zmPatch)
summary(zmNNDBuffer)
summary(zmNND)
summary(zmPatchBuffer)


