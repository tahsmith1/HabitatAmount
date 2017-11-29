

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

for(ii in  1:length(unique(comm.1$Patch))){
  pcomm.ii <- subset(comm.1, comm.1$Patch ==  unique(comm.1$Patch)[ii])
  pcomm.ii <- pcomm.ii[order(pcomm.ii$SurveyNumber),]
  chold[[ii]]<- pcomm.ii[c(1,3,5),]
  print(unique(comm.1$Patch)[ii])
}

chold <- rbindlist(chold)

pcomm1 <- subset(comm.1, comm.1$Patch == unique(comm.1$Patch)[1])
pcomm1 <- pcomm1[order(pcomm1$SurveyNumber),]
 chold[[1]]<- pcomm1[c(1,3,5),]
str(pcomm1$Survey.Date)
