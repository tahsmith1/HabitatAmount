####Reads in already calculated community and coordinate data frames#######


community<- read.csv("ManData/community.csv")
community2015 <- read.csv('ManData/community2015.csv')
community2015nocont <- read.csv("ManData/community2015nocont.csv")
communitySub <- read.csv("ManData/communitySub.csv")
coordinates<-  read.csv("ManData/coordinates.csv")


##############Basic Binomial regression

community2015nocont$CheliBinary <-ifelse(community2015nocont$Cheli.adults >0, 1,0)

glm_patchareaXbufferarea <- glm(CheliBinary~bufferarea*patcharea,family = binomial(link = logit),data = community2015nocont)

summary(glm_patchareaXbufferarea)
