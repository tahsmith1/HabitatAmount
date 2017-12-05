
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

zPatch<-zeroinfl(Cheli.adults ~ patcharea*year*TrtType | 1, data = comm.sub)
zBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*year*TrtType | 1, data = comm.sub)
zPatchBuffer<-zeroinfl(Cheli.adults ~ patcharea*bufferfararea*year*TrtType | 1, data = comm.sub)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferfararea*nnd*year*TrtType | 1, data = comm.sub)
zNND<-zeroinfl(Cheli.adults ~ nnd*year*TrtType | 1, data = comm.sub)
zNNDPatch<-zeroinfl(Cheli.adults ~ nnd*patcharea*year*TrtType | 1, data = comm.sub)

AIC(zPatch,zBuffer,zNND,zPatchBuffer,zNNDBuffer, zNNDPatch)


zmPatchBuffer<-glmmadmb(Cheli.adults~scale(bufferarea)*scale(patcharea)*TrtType*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
summary(zmPatchBuffer)


####################Testing patch area and nnd effects against bufferarea across treatments using unsubsetted data
zPatch<-zeroinfl(Cheli.adults ~ patcharea | 1, data = comm.nocont)
zBuffer<-zeroinfl(Cheli.adults ~ patcharea + bufferarea | 1, data = comm.nocont)
zNNDBuffer<-zeroinfl(Cheli.adults ~ bufferarea*nnd | 1, data = comm.nocont)
zNND<-zeroinfl(Cheli.adults ~ nnd | 1, data = comm.nocont)
zPropTrt<-zeroinfl(Cheli.adults ~PropLost*TrtType | 1, data = comm.nocont)


AICtable <- AIC(zPatch,zBuffer,zNNDBuffer,zNND,zPropTrt)
AICtable[order(AICtable$AIC),]

######ZIP mixed models using final dataset 

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

#####Created boxplots in ggplot2

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

####Adding year interaction into the analysis (picks up differences esp for trtmenttype)

mPatchBuffer<-glmmadmb(Cheli.adults~scale(patcharea)+scale(bufferarea)*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPatch<-glmmadmb(Cheli.adults~scale(patcharea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mNND<-glmmadmb(Cheli.adults~scale(patcharea)+scale(nnd)*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mNNDBuffer<-glmmadmb(Cheli.adults~scale(patcharea)+scale(nnd)*scale(bufferarea)+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPropTrt<-glmmadmb(Cheli.adults~scale(patcharea)+PropLost*TrtType+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)
mPropTrt<-glmmadmb(Cheli.adults~scale(patcharea)+PropLost*TrtType*year+(1|Plot), data=comm.sub, family="Poisson", zeroInflation=T)


AICtable <- AIC(mPatch,mPatchBuffer,mNND,mNNDBuffer ,mPropTrt)

summary(mNNDBuffer)
summary(mPropTrt)



