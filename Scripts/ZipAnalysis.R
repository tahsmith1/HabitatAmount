
########Begin analyses after sourcing pre analysis data manipulations 

library(pscl)
library(glmmADMB)

source("Scripts/PreAnalysisSource.R")

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

zmBuffer<-glmmadmb(Cheli.adults~scale(bufferfararea)*scale(patcharea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmPatchBuffer<-glmmadmb(Cheli.adults~scale(patcharea)+scale(bufferfararea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmPatch<-glmmadmb(Cheli.adults~scale(patcharea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmNND<-glmmadmb(Cheli.adults~scale(nnd)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmNNDBuffer<-glmmadmb(Cheli.adults~scale(nnd)*scale(bufferfararea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
zmPropTrt<-glmmadmb(Cheli.adults~PropLost*TrtType+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)

AICtable <- AIC(zmPatch,zmPatchBuffer,zmNND,zmNNDBuffer ,zmPropTrt)
summary(zmPropTrt)

AICtable$deltaAIC<- AICtable$AIC- min(AICtable$AIC)

plot(patcharea~bufferarea, data= comm.sub)

summary(zmNNDBuffer)

zmNNDBuffer<-glmmadmb(Cheli.adults~scale(nnd)*scale(bufferarea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T, mcmc = TRUE, mcmc.opts = mcmcControl(mcmc=1000))

boxplot(log(comm.sub$bufferarea)~comm.sub$year)  ####Box plots of bufferarea
boxplot(log(comm.sub$patcharea)~comm.sub$year)

library(ggplot2)

ggplot(comm.sub, aes(x = year, y= log(bufferarea)))+
  geom_boxplot( fill= "forestgreen") +
  theme_minimal()+
  theme(axis.title = element_text(size=20)) +
  xlab("Year") +
  ylab("log(Local Habitat Amount)")

ggplot(comm.sub, aes(x = year, y= log(patcharea)))+
  geom_boxplot( fill= "forestgreen") +
  theme_minimal()+
  theme(axis.title = element_text(size=20)) +
  xlab("Year") +
  ylab("log(Patch Area)")

ggplot(comm.sub, aes(x = year, y= log(nnd)))+
  geom_boxplot(fill= "forestgreen") +
  theme_minimal()+
  theme(axis.title = element_text(size=20)) +
  xlab("Year") +
  ylab("Nearest neighbor distance")

mPatchBuffer<-glmmadmb(Cheli.adults~scale(patcharea)+scale(bufferarea)*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPatch<-glmmadmb(Cheli.adults~scale(patcharea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mNND<-glmmadmb(Cheli.adults~scale(patcharea)+scale(nnd)*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mNNDBuffer<-glmmadmb(Cheli.adults~scale(patcharea)+scale(nnd)*scale(bufferarea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPropTrt<-glmmadmb(Cheli.adults~scale(patcharea)+PropLost*TrtType+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPropTrt<-glmmadmb(Cheli.adults~scale(patcharea)+PropLost*TrtType*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)


AICtable <- AIC(mPatch,mPatchBuffer,mNND,mNNDBuffer ,mPropTrt)

summary(mNNDBuffer)

mean(scale(comm.sub$patcharea))


pred_nndbuffer = data.frame(bufferarea = seq(min(scale(comm.sub$bufferarea)), max(scale(comm.sub$bufferarea)), length.out = 200)) 
pred_nndbuffer$nnd = median(comm.sub$nnd)
pred_nndbuffer$patcharea = median(scale(comm.sub$patcharea))
predict(mNNDBuffer, pred_nndbuffer, type = "response")


