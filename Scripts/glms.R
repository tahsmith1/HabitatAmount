####Reads in already calculated community and coordinate data frames#######


community<- read.csv("ManData/community.csv")
community2015 <- read.csv('ManData/community2015.csv')
community2015.nocont <- read.csv("ManData/community2015nocont.csv")
communitySub <- read.csv("ManData/communitySub.csv")
coordinates<-  read.csv("ManData/coordinates.csv")


##############Basic Binomial regression

community2015.nocont$CheliBinary <-ifelse(community2015.nocont$Cheli.adults >0, 1,0)

glm_patchareaXbufferarea <- glm(CheliBinary~bufferarea*patcharea,family = binomial(link = logit),data = community2015.nocont)

glm_nndXbufferarea <- glm(CheliBinary~bufferarea*nnd,family = binomial(link = logit),data = community2015.nocont)



summary(glm_patchareaXbufferarea)
summary(glm_nndXbufferarea)


cor(community2015.nocont$patcharea, community2015.nocont$bufferarea)  ### correlation


plot(glm_patchareaXbufferarea)

