####Reads in already calculated community and coordinate data frames#######
####Subsets to observations from Aug-Nov, using monthly observations in 2015 and includes only patches remaining from 2014 to 2015###
####Plots that experienced no habitat loss are removed and each plot is subsetted by mean dispersal distance to remove edge effects###

community<- read.csv("ManData/community.csv")
#community2015 <- read.csv('ManData/community2015.csv')
#community2015.nocont <- read.csv("ManData/community2015.nocont.csv")
#communitySub <- read.csv("ManData/communitySub.csv")
#landscape2014 <- read.csv("ManData/landscape2014.csv")
#landscape2015 <- read.csv("ManData/landscape2015.csv")

#################

#treatprop<- unique(community[,c("Plot","Patch","TrtType","PropLost")])

############ subset by Aug-Nov

comm2014<- subset(community, subset = community$year == 2014)
comm2015<- subset(community, subset = community$year == 2015)

comm2014an <- subset(comm2014, comm2014$month >7)
comm2015an <- subset(comm2015, comm2015$month >7)

#commlate <- rbind(comm2014an,comm2015an)

##########  Seperate calculations of only treated patches
#str(treatprop)
#lnd2014<- landscape2014[,-1]

#landsc2014<- merge(lnd2014,treatprop, by= c("Plot","Patch"), all = F)

#mergedcomm2014 <- aggregate(Cheli.adults~Plot+Patch, data= comm2014an, FUN= sum)
#comm2014m<- merge(mergedcomm2014,landsc2014, by = c("Plot", "Patch"), all = F)

#lnd2015<- landscape2015[,-1]
#landsc2015<- merge(lnd2015,treatprop, by= c("Plot","Patch"), all = F)

#mergedcomm2015 <- aggregate(Cheli.adults~Plot+Patch, data= comm2015an, FUN= sum)
#comm2015m<- merge(mergedcomm2015,landsc2015, by = c("Plot", "Patch"), all = F)

####Takes the 1st, 3rd, and 5th observation in 2015 to balance it with the 3 monthly obs in 2014
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

################# Removing patches in 2014 that don't remain in 2015 after treatment
#communityYear<- rbind(comm2014an,comm2015an3)

leftpatches <- comm2015an3[,c("Patch","Plot")]
leftpatches<- unique(leftpatches)
comm2014ms <- merge(comm2014an, leftpatches, by = c("Plot", "Patch"), all= F)

leftpatches1<- comm2014ms[,c("Patch","Plot")]
leftpatches1 <-unique(leftpatches1)
comm2015ms<- merge(comm2015an3, leftpatches1, by = c("Plot", "Patch"), all= F)

###Make combined commmunity data frame from 2014 and 2015 data

communityYear<- rbind(comm2014ms,comm2015ms)

########Plot, patch and year made into factors for analysis

communityYear$Plot<- as.factor(communityYear$Plot)
communityYear$Patch<- as.factor(communityYear$Patch)
communityYear$year<- as.factor(communityYear$year)


###Creates data frame without plots in which no habitat was removed (3 plots not used)
comm.nocont<- subset(communityYear, TrtType != "cont")

######### Reading in subsetmerge from subsetcalcs.R to only get patch >12,5m from edge
subsetmerge <- read.csv("ManData/subsetmerge.csv")
communityNew <-  merge(comm.nocont,subsetmerge , by=c("Plot","Patch"), all=F)
comm.sub <- subset(communityNew, coordInc == 1)

print("Ready for analysis:")
print("Data frame 'communityYear' for non-subsetted data (all plots, no edge effect subset)")
print("Data frame 'comm.sub' for subsetted data (of primary interest)")
