# EXAMINE HOW Re VARIES THROUGH TIME (AND SPACE), ACCORDING TO DIFFERENT FACTORS:
# Vaccination coverage
# Population density
# Previous incidence

rm(list=ls())
library(RColorBrewer)
library(plyr)
library(pscl)
require(lme4)
library(MASS) # library for negative binomial
library(rgdal)
library(codetools)
library(lubridate)
library(ggplot2)
library(reshape)
start = 1; end = 133 #Jan 2010 - Jan 2021

source("R/calcRe.R") # Function to examine how Re varies through time (and space) - maybe change to do moving avg?

# Import case data to match the tree data
rabid <- fread("Output/trees/case_dt_cleaned.csv") # rabid <- read.csv("Output/rabid.csv") # Created in Run_transmission_trees.R
dogs <- read.csv("Output/dogPopMatPembaWard.csv", header = FALSE); dim(dogs)
coverage <- read.csv("Output/vcMonthPembaWards.csv", header = FALSE); dim(coverage)
names(coverage)[2:ncol(coverage)]<-1:(ncol(coverage)-1)
names(coverage)[1] <- "shehia_id" # Rename the village column

rows = match(tolower(rabid$Loc_ID), tolower(coverage$shehia_id), nomatch=0)
list(rows)
which(rows==0)


# Re organise the vaccination coverage data and rename the columns
coverage_new <- melt(coverage, id="shehia_id")
dim(coverage_new)
names(coverage_new)[names(coverage_new)=="value"] <- "vacc_cov"
names(coverage_new)[names(coverage_new)=="variable"] <- "month"
coverage_new$month <- as.numeric(coverage_new$month)
head(coverage_new)


trees <- readRDS("output/sources_lnorm_combo.rda")
LL=readRDS("output/ll_lnorm_combo.rda")

# Examine where tree building algorithm failed:
nas = which(is.na(apply(trees, 1, sum)))
NAS = rep(NA, length(nas))
for(i in 1:length(nas)){
  NAS[i] <- length(which(!is.na(trees[nas[i],])))
}
rabid$ID[nas[which(NAS == 0)]] # The IDs that have no ancestor and need sorting out!

#Find most likely progenitors for each case (and their bootstrap support):
cases <- nrow(rabid)
runs <- ncol(trees)

#' _: Warning appears when i=185 - two biters (520 & 526) are possible_
biter <- bootstrap <- progll <- numeric(cases)
for(i in 1:cases){
  if(length(which(is.na(trees[i,])))!=runs){   # If not all results are NA, assign most likely progenitor
    biters = tabulate(t(trees[i,]))
    biter[i] <- which(biters == max(biters))
    bootstrap[i] <- max(biters)
#     print(i)
  }
}

length(which(bootstrap!=0))/cases   # % of assigned progenitors

# Assign Re to each rabid animal
nc = hist(biter, breaks=-1:max(rabid$ID, na.rm=TRUE), plot=T)$counts[-1]
rabid$Re <- nc[match(rabid$ID, 1:max(rabid$ID))]
mean(rabid$Re)

# Examine how Re varies with space and time
times = seq(start, end, 3)
Re = calcRe(start, end, 3, rabid$Re, rabid$ms)
Re$times = times

pdf("figures/Re_timeseries.pdf")
par(mfrow=c(1,1), mar=c(4,3,3,3))
plot(times, Re$mean, type="l", lwd=2, ylim = c(0,8),xlab="", ylab="Mean Re",axes=F)
axis(2)
axis(1, at=seq(1, 133, by=12),labels=2010:2021)
axis(1, at =times, labels=rep("", length(times)), lwd=0.5, tck=-0.01)
points(rabid$ms, jitter(rabid$Re), pch=20, ylim=c(0,20), cex=0.6, col="dark gray") #
lines(times, Re$UCI, col="red")
lines(times, Re$LCI, col="red")
lines(times, Re$CI_90, col="red", lwd=0.5)
lines(times, Re$CI_10, col="red", lwd=0.5)
lines(times, rep(1, length(times)), col="gray")
dev.off()

rabies_ts = hist(rabid$ms, breaks=(start-1):end, plot=FALSE)$counts
plot(start:end, rabies_ts, type="l", xlab="",axes=F)
axis(2)
axis(1, at=seq(1, 133, by=12),labels=2010:2021)
axis(1, at =times, labels=rep("", length(times)), lwd=0.5, tck=-0.01)
hist(Re$mean, breaks=seq(0,3,0.5))
mean(Re$mean, na.rm=TRUE)
range(Re$mean, na.rm=TRUE)


#' Examine how Re varies between vaccination campaigns
vc = c(17,32,50,60,80,93,106,117) # Vaccination campaigns Apr-May 2011 (17); July-Aug 2012 (32); 
# Oct 2013 plus jan-feb 2014 (50), dec 2014 (60),2016 july 80,2017 july 93, 2018 july 106, 
# 2019 aug 117

vacc_interval <-rep(NA, nrow(rabid))
vacc_interval[which(rabid$ms > 0 & rabid$ms <=17)]  <- "Pre vaccination"
vacc_interval[which(rabid$ms > 17 & rabid$ms <= 32)]<-"1st Round" # Aug 2012
vacc_interval[which(rabid$ms > 32 & rabid$ms <=50)] <- "2nd Round"  # Feb 2014
vacc_interval[which(rabid$ms > 50 & rabid$ms <=60)] <- "3rd Round"  # Dec 2014
vacc_interval[which(rabid$ms > 60 & rabid$ms <=80)] <- "4th Round"  # July 2016
vacc_interval[which(rabid$ms > 80 & rabid$ms <=93)] <- "5th Round"  # July 2017
vacc_interval[which(rabid$ms > 93 & rabid$ms <=106)] <-"6th Round" # July 2018
vacc_interval[which(rabid$ms > 106 & rabid$ms <=117)]<-"7th Round"  # Aug 2019
vacc_interval[which(rabid$ms > 117 )] <- "8th Round"
vacc_interval = as.factor(vacc_interval)
print(levels(vacc_interval))

rabid$vac<-vacc_interval
rabid$vac<- as.factor(rabid$vac)

# Create a function for caluclating the mean, stadndard deviation and error together with their confidence intervals
meancom <-function(x, k)
{ 
  Mn<-tapply(x, list(k), mean)
  sD<-tapply(x, list(k), sd)
  lD<-tapply(x, list(k), length)
  se<-sD/sqrt(lD)
  print("95% CI")
  Lower = Mn - qnorm(0.975)*se
  Upper = Mn + qnorm(0.975)*se
  All<-cbind(Mn, Lower, Upper)
  print(All)

}

 meancom(rabid$Re, rabid$vac)

allmtr<-as.data.frame(meancom(rabid$Re, rabid$vac))
allmtr$rounds<-rownames(allmtr)
allmtr$rounds<-as.factor(allmtr$rounds)

ggplot(allmtr, aes(x=rounds, y=Mn)) +
geom_bar(aes(x=rounds, y=Mn), stat="identity", fill="skyblue", alpha=0.5) +
geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, colour="orange", alpha=0.9, size=1.3)+
 xlab("Rounds") + ylab("Mean")


#'*Examine how Re (for individual rabies cases) varies with vaccination coverage, 1st*
#'*add vaccination coverage to the individual rabies cases (RE)*
#'
rabid$coverage <- NA
for (n in 1:nrow(rabid)){
  
  #' _RS: Remember to comment what your code is doing - it can help you identify issues_
  #' _You had tolower() in here, which meant that no location names matched_
  index = which(coverage_new$shehia_id == rabid$Loc_ID[n] & coverage_new$month == rabid$ms[n])
  # which(coverage_new$shehia_id == rabid$shehia[n] & coverage_new$month == rabid$ms[n])
  rabid$coverage[n] <- ifelse(length(index)==0, NA, coverage_new$vacc_cov[index])
} 

#  Write the rabid cases (RE)  and coverage  to a csv file
 write.csv(rabid, "Output/Re.coverage1.csv")

#  Running a glm to examine how vaccination coverage affects Re
#  Fisrt categorize vaccination coverage into either low, medium or high


cfactor <-rep(NA, nrow(rabid))
cfactor[which(rabid$coverage <= 20)] <- "Low"
cfactor[which(rabid$coverage > 20 & rabid$coverage < 70)] <- "Medium"
cfactor[which(rabid$coverage >= 70 )] <- "High"

cfactor = as.factor(cfactor)
print(levels(cfactor))
rabid$cfactor <- cfactor


#  Run and select the best model (GLMM)
#  USE THE DROP ONE FUNCTION TO OMIT THE LEAST SIGNIFICANT VALUE FROM THE MODEL UNTIL I REMAIN WITH THE BEST MODEL
 rabid$re = as.factor(1:nrow(rabid)) # residual error term

#  m1 <- glmer(Re~ cfactor*ms + dogpop + ms +
#              (1|ms), family = poisson, data = rabid)
#   summary(m1) #
#  drop1 (m1, test="Chi") # drop the interaction btn months and vaccination coverage ..its the least significant
# 
#  m2 <- glmer(Re~ cfactor+ ms + dogpop  +(1|re)+
#              (1|ms), family = poisson, data = rabid)
#  summary(m2) # Not significant at all
#   drop1 (m2, test="Chi") # drop months
# 
# # m3 <- glmer(Re~ cfactor+ dogpop  +(1|re) +
#             (1|ms), family = poisson, data = rabid)
# summary(m3) # Not signficant
# drop1 (m3, test="Chi") # drop dog population that is the least significant
#
#
# m4<-glmer(Re~ cfactor+(1|re)+
#             (1|ms), family = poisson, data = rabid)
# summary(m5)# Not all significant
#
#
# # EXAMINE THE RELATIONSHIP BETWEEN RE AND THE INTERVAL BTN CAMPAIGNS
# m6 <- glm(Re~ vacc_interval, family = poisson, data = rabid)
# summary(m6) # DOES NOT SHOW ANY RELATIONSHIP, over dispoersed
#
#
# m7<-glm.nb(Re~vacc_interval, data = rabid)
# summary(m7) # changed to negative binomial still insignificant
#
# m8= glm(Re~ cfactor, family=poisson, data=rabid)
# summary(m8)
#
# m9<-glm.nb(Re~cfactor, data = rabid)
# summary(m9) # changed to negative binomial still insignificant
#
# # Examine how Re is affected with time
# m10 <- glm(Re~ms,family= poisson,data=rabid)
# summary(m10) # a significant and over dispispersed, change to n.binomial error
#
# m11 <- glm.nb (Re~ms,data=rabid)
# summary(m11) # a significant and over dispispersed, change to n.binomial error
#
#
# m12 <- glm.nb(Re ~ coverage, data = rabid)
# summary(m12)
#
# # uSING THE ZERO INFLATION MODEL TO ACCOUNT FOR THE MANY ZERO VALUES
# m13 <- zeroinfl(Re ~ coverage, data = rabid)
# summary(m13) # everything  not significant
#
#
#
