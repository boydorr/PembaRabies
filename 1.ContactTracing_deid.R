#---------------------------------------------------------
# Author: Kennedy Lushasi
# Explores exposure and rabies data - compare contact tracing vs hospital records, PEP use & adherance
#---------------------------------------------------------
rm(list=ls())

## load all libraries
library(lubridate)
library(tidyverse)
library(Hmisc) # binomial probability of death if no PEP
library(fitdistrplus)
library(MASS)
library(boot)
options(stringsAsFactors=F)

#---------------------------------------------------------
#  Set start dates
startDate <- as.Date("2010-01-01")

## Read in data 
animals <- readRDS(file = "output/animals_deid.rda")  # animal cases & animals bitten
cases <- readRDS(file = "output/rabid_deid.rda") # animal cases
humans <- readRDS(file = "output/humans_deid.rda") # human exposures & bite patients
humanpop <- 472958 # from 2012 census 

#---- Process data -------------------------------------------------------------
# Identify if received PEP
humans$no.PEP <- (humans$If.no.PEP..why.not.Sought.but.not.available=="true" |
                    humans$If.no.PEP..why.not.Not.aware=="true" |
                    humans$If.no.PEP..why.not.Too.expensive=="true" |
                    humans$If.no.PEP..why.not.Not.advised.by.health.worker=="true" |
                    humans$If.no.PEP..why.not.Not.advised.by.dog.owner=="true" |
                    humans$If.no.PEP..why.not.Didn.t.think.was.needed=="true" |
                    humans$If.no.PEP..why.not.Currently.seeking....outcome.unknown=="true")
# or got late PEP or were searching
humans$late.PEP <- (humans$PEP.Status.Raising.funds=="true"| humans$PEP.Status.Seeking=="true")

# Table pre-/post-outbreak
table(humans$Year.bitten) # None in 2015
humans$outbreak <- ifelse(humans$Year.bitten < 2015, "Pre-outbreak",
                          ifelse(humans$Year.bitten > 2015, "Post-outbreak", NA))
table(humans$outbreak, useNA="always") # bite records (healthy & rabid)
exposures <- subset(humans, Rabid == "Yes"); nrow(exposures) # Subset of exposures (known rabid)

#---------- Calculate annual incidence -----------------------------------
# exposures, deaths and patients (went for PEP) for rabid & healthy bites by year
bites_y <- humans %>%
  mutate(Year = factor(Year.bitten, levels = 2010:2021)) %>%
  group_by(Year, .drop=FALSE)  %>%
  dplyr::summarize(
    exposures = length(which(Rabid=="Yes")), # likely exposures
    poss_exposures = length(which(Rabid == "Unknown")), # possible exposures 
    healthy = length(which(Rabid=="No" & PEP.1=="true")), # Healthy bites that received PEP (excludes healthy bite victims found through contact tracing)
    deaths = length(which(Patient.outcome=="Died")), # Rabies exposures who died
    exp_inc = exposures*100000/humanpop, # Exposure incidence
    death_inc = deaths*100000/humanpop, # Death incidence
    PEP1 = length(which(PEP.1=="true")), # Patients (rabid & healthy) that received PEP dose 1
    PEP2 = length(which(PEP.2=="true")), # dose 2
    PEP3 = length(which(PEP.3=="true")), # dose 3
    noPEP_rabid = length(which(Rabid=="Yes" & no.PEP == TRUE)), # Exposures that did not receive PEP
    noPEP_rabid_poss = length(which(Rabid=="Unknown" & no.PEP == TRUE)), # Possible exposures that did not receive PEP
    latePEP_rabid = length(which(Rabid=="Yes" & late.PEP == TRUE)), # Exposures that received PEP late
    latePEP_rabid_poss = length(which(Rabid=="Unknown" & late.PEP == TRUE)), # Possible exposures that received PEP late
    rabid_noPEPsought = length(which(Rabid=="Yes" & PEP.Status.Not.sought=="true")), # Exposures that did not seek PEP
    possrabid_noPEPsought = length(which(Rabid=="Unknown" & PEP.Status.Not.sought=="true")), # Exposures that did not seek PEP
    rabid_PEPcomplete = length(which(Rabid=="Yes" & PEP.Status.Completed=="true")), # Exposures that completed PEP
    rabid_PEPIncomplete =length(which(Rabid=="Yes" & PEP.Status.Incomplete=="true"))) # Exposures that did not complete PEP
# View table
View(bites_y)
write.csv(bites_y, "output/bite_summary.csv", row.names=FALSE)

# Probability of completing next dose
sum(bites_y$PEP2)/sum(bites_y$PEP1) # 91% chance of completing 2nd dose if start 
sum(bites_y$PEP3)/sum(bites_y$PEP1)

# Alternative table: pre-outbreak and post-outbreak
bites_outbreak_gr <- humans %>%
  mutate(outbreak = factor(outbreak, levels = c("Pre-outbreak", "Post-outbreak"))) %>%
  group_by(outbreak)  %>%
  dplyr::summarize(
    exposures = length(which(Rabid=="Yes")), # likely exposures
    poss_exposures = length(which(Rabid == "Unknown")), # possible exposures - removed Rabid=="Yes" | as other "poss" columns exclude this
    healthy = length(which(Rabid=="No" & PEP.1=="true")), # Healthy bites that received PEP
    deaths = length(which(Patient.outcome=="Died")), # Exposures that died of probable rabies
    exp_inc = round(exposures*100000/humanpop, digits=3), # Exposure incidence
    death_inc = round(deaths*100000/humanpop, digits=3), # Death incidence
    PEP1 = length(which(PEP.1=="true")), # patients (rabid & healthy) that received PEP dose 1
    PEP2 = length(which(PEP.2=="true")), # dose 2
    PEP3 = length(which(PEP.3=="true")), # dose 3
    noPEP_rabid = length(which(Rabid=="Yes" & no.PEP == TRUE)), # rabid exposures that did not receive PEP
    noPEP_rabid_poss = length(which(Rabid=="Unknown" & no.PEP == TRUE)), # possible exposures that did not receive PEP
    noPEP_healthy = length(which(Rabid=="No" & no.PEP == TRUE)), # Healthy bites that did not receive PEP
    latePEP_rabid = length(which(Rabid=="Yes" & late.PEP == TRUE)), # rabid exposures that received PEP late
    latePEP_rabid_poss = length(which(Rabid=="Unknown" & late.PEP == TRUE)), # possible exposures that received PEP late
    latePEP_healthy = length(which(Rabid=="No" & late.PEP == TRUE)), # Healthy bites that received PEP late
    noPEPsought_rabid = length(which(Rabid=="Yes" & PEP.Status.Not.sought=="true")), # rabid exposures that did not seek PEP
    noPEPsought_possrabid = length(which(Rabid=="Unknown" & PEP.Status.Not.sought=="true")), # possible exposures that did not seek PEP
    noPEPsought_healthy = length(which(Rabid=="No" & PEP.Status.Not.sought=="true")), # Healthy bites that did not seek PEP
    PEPcomplete_rabid = length(which(Rabid=="Yes" & PEP.Status.Completed=="true")), # rabid exposures that completed their PEP treatment
    PEPcomplete_possrabid = length(which(Rabid=="Unknown" & PEP.Status.Completed=="true")), # possible exposures that completed their PEP treatment
    PEPcomplete_healthy = length(which(Rabid=="No" & PEP.Status.Completed=="true")), # Healthy bites that completed their PEP treatment
    PEPIncomplete_rabid =length(which(Rabid=="Yes" & PEP.Status.Incomplete=="true")), # rabid exposures that did not complete their PEP treatment
    PEPIncomplete_possrabid =length(which(Rabid=="Unknown" & PEP.Status.Incomplete=="true")), # possible exposures that did not complete their PEP treatment
    PEPIncomplete_healthy =length(which(Rabid=="No" & PEP.Status.Incomplete=="true"))) # Healthy bites that did not complete their PEP treatment
# View table
View(bites_outbreak_gr)

# Transform to long dataframe
bites_outbreak_gr_t <- t(bites_outbreak_gr)
colnames(bites_outbreak_gr_t) <- bites_outbreak_gr_t[1,]
bites_outbreak_gr_t <- bites_outbreak_gr_t[-1,]
write.csv(bites_outbreak_gr_t, "output/outbreak_table.csv") # Save as output

# Ratio of healthy: exposed bites
bites_outbreak_gr$exposures/bites_outbreak_gr$healthy

# Examine PEP access and deaths
n <- sum(bites_y$exposures); n # Total rabid exposures
n2 <- sum(bites_y$poss_exposures, bites_y$exposures); n2 # Total possible exposures
deaths <- subset(humans, Patient.outcome == "Died");
nDeaths <- nrow(deaths); nDeaths; sum(bites_y$deaths) # 6 deaths
noPEP <- sum(bites_y$noPEP_rabid); noPEP2 <- sum(bites_y$noPEP_rabid_poss) # 60 rabid/1 possible did not complete timely PEP
pPEP <- (n-noPEP)/n; round(pPEP, 2); # probability of PEP for definite exposures
pPEP2 <- (n2-(noPEP+noPEP2))/n2; round(pPEP2, 2) # 0.75-0.78 probability of PEP

# Deaths
deaths$Date.bitten; deaths$When.died;
d1 <- which(deaths$Year.bitten == 2010)
d2 <- which(deaths$Year.bitten == 2017)
as.Date(deaths$When.died)-as.Date(deaths$Date.bitten) # Uncertainty about case 5 - mid Dec 2010 but no exact dates (died 2 weeks after bite)

# Examine bite locations for death cases
message("Head: ", paste0(deaths$Bite.site.Head...neck, collapse=", ")); 
message("Arms: ", paste0(deaths$Bite.site.Arms, collapse=", ")); 
message("Hands: ", paste0(deaths$Bite.site.Hands, collapse=", ")); 
message("Trunk: ", paste0(deaths$Bite.site.Trunk, collapse=", ")); 
message("Legs: ", paste0(deaths$Bite.site.Legs, collapse=", ")); 
message("Feet: ", paste0(deaths$Bite.site.Feet, collapse=", "))
message("PEP 1: ", paste0(deaths$PEP.1.Date, collapse=", ")); 
message("PEP 2: ", paste0(deaths$PEP.2.Date, collapse=", ")); 
message("PEP 3: ", paste0(deaths$PEP.3.Date, collapse=", "))
deaths$PEP.Status.Not.sought[d1]; deaths$PEP.Status.Completed[d1]; deaths$PEP.Status.Incomplete[d1]
deaths$PEP.Status.Not.sought[d2]; deaths$PEP.Status.Completed[d2]; deaths$PEP.Status.Incomplete[d2]

# Examine proportion of patients due to healthy vs rabid dogs
bites_y$Year <- as.numeric(levels(bites_y$Year))[bites_y$Year]
p1 <- which(bites_y$Year < 2015)
p2 <- which(bites_y$Year > 2015)

# Note - rabid bite victims that did not seek care are NOT patients!
total_patients_p1 <- sum(bites_y$healthy[p1]) + sum(bites_y$exposures[p1]) + sum(bites_y$poss_exposures[p1]) - sum(bites_y$rabid_noPEPsought[p1]); total_patients_p1
total_patients_p2 <- sum(bites_y$healthy[p2]) + sum(bites_y$exposures[p2]) + sum(bites_y$poss_exposures[p2]) - sum(bites_y$rabid_noPEPsought[p2]); total_patients_p2
pExposures1 <- (sum(bites_y$exposures[p1]) - sum(bites_y$rabid_noPEPsought[p1]))/total_patients_p1; pExposures1
pExposures_uk <- (sum(bites_y$poss_exposures[p1]) - sum(bites_y$possrabid_noPEPsought[p1]))/total_patients_p1; pExposures_uk + pExposures1
pExposures2 <- (sum(bites_y$exposures[p2]) - sum(bites_y$rabid_noPEPsought[p2]))/total_patients_p2; pExposures2

# Incidence of bite by healthy dogs (includes possible exposures that add uncertainty)
mean(bites_y$healthy[p1]);
mean(bites_y$healthy[p1] + bites_y$poss_exposures[p1]) # bites_y$unknowns[p1])
mean(bites_y$healthy[p2]); 
healthy <- data.frame(
  BI1 = mean(bites_y$healthy[p1])*100000/humanpop,
  BI1min = min(bites_y$healthy[p1])*100000/humanpop,
  BI1max = max(bites_y$healthy[p1]+bites_y$poss_exposures[p1])*100000/humanpop,
  BI1uk = mean(bites_y$healthy[p1] + bites_y$poss_exposures[p1])*100000/humanpop, #bites_y$unknowns[p1])*100000/humanpop
  BI2 = mean(bites_y$healthy[p2])*100000/humanpop,
  BI2min = min(bites_y$healthy[p2])*100000/humanpop,
  BI2max = max(bites_y$healthy[p2])*100000/humanpop
)
write.csv(healthy, "output/healthy_bites.csv", row.names = FALSE)

# Numbers of patients that likely did not get PEP
params <- read.csv("data/bio_data.csv") # parameters for protective efficacy of vaccine from Changalucha et al; WHO modelling consortium
pRabies <- params$p_rabies_transmission
nDeaths/pRabies # would expect ~36 to not obtain PEP
pDeath_range <- binconf(sum(bites_y$exposures[p1])*pRabies, sum(bites_y$exposures[p1]))
pDeath_range_uk <- binconf(sum(bites_y$poss_exposures[p1])*pRabies, sum(bites_y$poss_exposures[p1])) # This now correctly only uses unknowns, previously included rabid exposures
pDeath_range; pDeath_range_uk

sum(bites_y$noPEP_rabid[p1]) 
sum(bites_y$deaths[p1])/c(pDeath_range,pDeath_range_uk) # expect this many people to not seek treatment (prior to 2015)!
sum(bites_y$rabid_noPEPsought[p2]) # during the 2016 outbreak all sought PEP
sum(bites_y$noPEP_rabid[p2]) # but 39 didn't obtain (adequate) PEP - SO WHAT DID THEY GET? either late or incomplete PEP

#' _Exposures who did not get PEP with reasons_
n1 <- sum(bites_y$exposures[p1]); n1
n2 <- sum(bites_y$exposures[p2]); n2
pre <- which(exposures$Year.bitten<2015)
post <- which(exposures$Year.bitten>2015)
shortages1 <- table(exposures$If.no.PEP..why.not.Sought.but.not.available[pre])
shortages2 <- table(exposures$If.no.PEP..why.not.Sought.but.not.available[post])
round(shortages1 *100/n1, 2)
round(shortages2 *100/n2, 2)
# No pep shortage before the outbreak, but later 0.56% ran into shortages

unaware1 <- table(exposures$If.no.PEP..why.not.Not.aware[pre])
unaware2 <- table(exposures$If.no.PEP..why.not.Not.aware[post])
round(unaware1*100/n1,2)
round(unaware2*100/n2,2)
# 7.94% of exposures were not aware of PEP before outbreak, but after the outbreak it decreased to 1.67%

expense1 <- table(exposures$If.no.PEP..why.not.Too.expensive[pre])
expense2 <- table(exposures$If.no.PEP..why.not.Too.expensive[post])
round(expense1*100/n1,2)
round(expense2*100/n2,2)
# Before the outbreak cost was not a problem as PEP was provided at no cost, but during the outbreak, 0.56% failed to get PEP due to cost

illadvised1 <- table(exposures$If.no.PEP..why.not.Not.advised.by.health.worker[pre])
illadvised2 <- table(exposures$If.no.PEP..why.not.Not.advised.by.health.worker[post])
round(illadvised1*100/n1,2)
round(illadvised2*100/n2,2)
# 19.05% of exposures were not advised by HW before the outbreak vs 2.78 during the outbreak

dogowner1 <- table(exposures$If.no.PEP..why.not.Not.advised.by.dog.owner[pre])
dogowner2 <- table(exposures$If.no.PEP..why.not.Not.advised.by.dog.owner[post])
round(dogowner1*100/n1,2)
round(dogowner2*100/n2,2)
# During endemic period, no exposures failed to get PEP due not being advised by dog owners, but during the outbreak 0.56 did not get pep as were not advised by dog owners

#### SAME AS NOT BEING AWARE!
table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[pre])
round(table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[pre])*100/n1,2)
table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[post])
round(table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[post])*100/n2,2)
# 12.7% did not think PEP was need during endemic period Vs 2.78 during the outbreak

seeking1 <- table(exposures$If.no.PEP..why.not.Currently.seeking....outcome.unknown[pre])
seeking2 <- table(exposures$If.no.PEP..why.not.Currently.seeking....outcome.unknown[post])
round(seeking1*100/n1,2)
round(seeking2*100/n2,2)
# 14.44% were still searching for pep at the time of investigation during the outbreak, and none during endemic period

table("PEP complete before 2015:"=exposures$PEP.Status.Complete[pre])
table("PEP incomplete before 2015:"=exposures$PEP.Status.Incomplete[pre])
table("PEP complete after 2015:"=exposures$PEP.Status.Complete[post])
table("PEP not sought after 2015:"=exposures$PEP.Status.Not.sought[post])
table("PEP late after 2015:"=exposures$late.PEP[post])

# Probability of starting PEP
PEP_params <- data.frame(
  denom = c(n1, n2),
  pStart = c(n1- sum(bites_y$noPEP_rabid[p1]), n2 - sum(bites_y$noPEP_rabid[p2]))/c(n1,n2), # low probability
  pStart2 = c(sum(bites_y$rabid_PEPcomplete[p1], bites_y$rabid_PEPIncomplete[p1]),
              sum(bites_y$rabid_PEPcomplete[p2], bites_y$rabid_PEPIncomplete[p2]))/c(n1,n2), # high probability
  pComplete = c(n1- sum(bites_y$rabid_PEPIncomplete[p1]), n2 - sum(bites_y$rabid_PEPIncomplete[p2]))/c(n1,n2), # Probability of completing PEP
  denom_PEP = c(sum(bites_y$rabid_PEPIncomplete[p1],bites_y$rabid_PEPcomplete[p1]), sum(bites_y$rabid_PEPIncomplete[p2],bites_y$rabid_PEPcomplete[p2])),
  pComplete2 = c(sum(bites_y$rabid_PEPcomplete[p1])/sum(bites_y$rabid_PEPIncomplete[p1],bites_y$rabid_PEPcomplete[p1]),
                sum(bites_y$rabid_PEPcomplete[p2])/sum(bites_y$rabid_PEPIncomplete[p2],bites_y$rabid_PEPcomplete[p2]))
  )
PEP_params
write.csv(PEP_params, "output/PEP_params_Pemba.csv")


#----- Table of individual human deaths, ages, and reason for not obtaining PEP -----
deaths_table <- data.frame("Year"=substr(deaths$When.died, 1, 4), # year(as.Date(deaths$When.died)) # better way to pull out year? (as.Date)
                           "Age"=deaths$Age..in.years.,
                           "Head"=deaths$Bite.site.Head...neck,
                           "Arms"=deaths$Bite.site.Arms,
                           "Hands"=deaths$Bite.site.Hands,
                           "Trunk"=deaths$Bite.site.Trunk,
                           "Legs"=deaths$Bite.site.Legs,
                           "Feet"=deaths$Bite.site.Feet,
                           "Severe_broken_bones"=deaths$Bite.details.Severe..broken.bones.,
                           "Severe_hospitalisation"=deaths$Bite.details.Severe..hospitalization.,
                           "Large_wound"=deaths$Bite.details.Large.wound.s.,
                           "Minor_wound"=deaths$Bite.details.Minor.wound.s.,
                           "Scratch"=deaths$Bite.details.Scratch,
                           "Currently_seeking"=deaths$If.no.PEP..why.not.Currently.seeking....outcome.unknown,
                           "Thought_unecessary"=deaths$If.no.PEP..why.not.Didn.t.think.was.needed,
                           "Not_advised_by_owner"=deaths$If.no.PEP..why.not.Not.advised.by.dog.owner,
                           "Not_advised_by_HW"=deaths$If.no.PEP..why.not.Not.advised.by.health.worker,
                           "Not_aware"=deaths$If.no.PEP..why.not.Not.aware,
                           "Not_available"=deaths$If.no.PEP..why.not.Sought.but.not.available,
                           "Too_expensive"=deaths$If.no.PEP..why.not.Too.expensive)

# Transform "true" to the column heading
deaths_table$Head <- ifelse(deaths_table$Head=="true", "Head", NA)
deaths_table$Arms <- ifelse(deaths_table$Arms=="true", "Arm", NA)
deaths_table$Hands <- ifelse(deaths_table$Hands=="true", "Hand", NA)
deaths_table$Trunk <- ifelse(deaths_table$Trunk=="true", "Trunk", NA)
deaths_table$Legs <- ifelse(deaths_table$Legs=="true", "Leg", NA)
deaths_table$Feet <- ifelse(deaths_table$Feet=="true", "Foot", NA)
deaths_table$Severe_broken_bones <- ifelse(deaths_table$Severe_broken_bones=="true", "Severe (broken bones)", NA)
deaths_table$Severe_hospitalisation <- ifelse(deaths_table$Severe_hospitalisation=="true", "Severe (hospitalisation)", NA)
deaths_table$Large_wound <- ifelse(deaths_table$Large_wound=="true", "Large wound", NA)
deaths_table$Minor_wound <- ifelse(deaths_table$Minor_wound=="true", "Minor wound", NA)
deaths_table$Scratch <- ifelse(deaths_table$Scratch=="true", "Scratch", NA)
deaths_table$Currently_seeking <- ifelse(deaths_table$Currently_seeking=="true", "Currently seeking", NA)
deaths_table$Thought_unecessary <- ifelse(deaths_table$Thought_unecessary=="true", "Thought unecessary", NA)
deaths_table$Not_advised_by_owner <- ifelse(deaths_table$Not_advised_by_owner=="true", "Not advised by dog owner", NA)
deaths_table$Not_advised_by_HW <- ifelse(deaths_table$Not_advised_by_HW=="true", "Not advised by Healthworker", NA)
deaths_table$Not_aware <- ifelse(deaths_table$Not_aware=="true", "Not aware", NA)
deaths_table$Not_available <- ifelse(deaths_table$Not_available=="true", "Not available", NA)
deaths_table$Too_expensive <- ifelse(deaths_table$Too_expensive=="true", "Too expensive", NA)

# Paste results together in new column, ignoring NAs
deaths_table <- deaths_table %>%
  unite("Bite_location", Head:Feet, sep=", ", remove = F, na.rm = T)
deaths_table <- deaths_table %>%
  unite("Bite_injury", Severe_broken_bones:Scratch, sep=", ", remove = F, na.rm = T)
deaths_table <- deaths_table %>%
  unite("Reason_for_no_PEP", Currently_seeking:Too_expensive, sep=", ", remove = F, na.rm = T)
View(deaths_table) # View table to check

# Remove surplus columns, and arrange by year then age
deaths_table <- deaths_table %>%
  dplyr::select(Year, Age, Bite_location, Bite_injury, Reason_for_no_PEP) %>%
  arrange(Year, Age)

# Save output
write.csv(deaths_table, "output/Table_deaths.csv", row.names = F)


#-----Disease incidence in dogs per year--------
# Cases per year and between different species
rabid <- subset(animals, Suspect=="Yes") # only rabid animals
table(Year = rabid$Year, Species=rabid$Species) # does not list all years (missing 2015, 2019-2021 because no cases)
table(rabid$Species)*100/nrow(rabid) # proportion of each rabid species

dogpopulation <- 3156 # dog population has been adjusted from 4000 to the current estimates
cases_y <- rabid %>%
  mutate(Year = factor(Year, levels = 2010:2020)) %>%
  group_by(Year, .drop=FALSE)  %>%
  dplyr::summarize(n=n(),
                   dogs = length(which(Species=="Domestic dog")),
                   cows = length(which(Species=="Livestock: Cow")),
                   goats = length(which(Species=="Livestock: Goat"))) # sample positivity incidence - actually tractable as %
cases_y$VaxRound = c(0,1,2,3,4,5,0,1,2,3,4)
cases_y$Year <- as.numeric(levels(cases_y$Year))[cases_y$Year]
cases_y
write.csv(cases_y, "output/cases_y.csv", row.names = FALSE) 

# Look at dog rabies cases per year and case detection extrapolation
detectQs <- read.csv("Output/DetectQs.csv")
q50 <- which(detectQs$X== "50%")
pDetect <- c(detectQs$TD1[q50], detectQs$TD2[q50]) # detection probability of cases in Pemba before 2015 and after # 
sum(cases_y$n); sum(cases_y$dogs);

p1_exposures <- c(sum(bites_y$exposures[p1]), sum(bites_y$exposures[p1])+sum(bites_y$poss_exposures[p1]))
p2_exposures <- sum(bites_y$exposures[p2])

rabid_dogs <- c(sum(cases_y$dogs[p1]), sum(cases_y$dogs[p2], na.rm=TRUE))
est_rabid_dogs <- rabid_dogs/pDetect; est_rabid_dogs # Adjust rabid dog estimates given case detection

# Print confidence intervals on case detection/ outbreak size
rabid_dogs[1]/detectQs$TD1
rabid_dogs[2]/detectQs$TD2

# Ratio of exposures per rabid dog (based on adjusted rabid dog cases according to case detection)
p1_exposures/est_rabid_dogs[1]
p2_exposures/est_rabid_dogs[2]
p1_exposures; p2_exposures

#----- Cases/ bites/ per rabid dog -----------
# How many human exposures per rabid dogs on average
n.exposures <- sum(bites_y$exposures); n.exposures
n.cases <- sum(cases_y$n) ; n.cases  
n.exposures/n.cases # 1.3 - NOT ADJUSTED FOR CASE DETECTION  - think this can be deleted!
length(unique(exposures$Biter.ID)) # How many dogs are responsible for all the bites?
n.exposures/length(unique(exposures$Biter.ID)) # Biting dogs bite ~how many people each? 
which(is.na(exposures$Biter.ID)) # No exposed persons with NA Biter.ID 

# Function that takes the rabid animals and the bite victims and calculates numbers bitten per rabid animal
bites <- function(biters, bite_victims){
  rmax = max(biters, na.rm=TRUE) + 1 # Find the range of the IDs of the bite victims
  bitten_index = hist(bite_victims, breaks = -1:rmax, plot=F)$counts[-1] # aggregate bite victims by biter ID
  NB = rep(NA, length(biters)) # create a vector to store the numbers bitten (NB)
  for (i in 1:length(biters)){
    NB[i] = bitten_index[biters[i]] # match the numbers bitten to each biter
  }
  NB
}

# Alternative function to check biting!
bites2 <- function(biters, bite_victims){
  bites_table = table(bite_victims) 
  NB = rep(0, length(biters))
  index = match(names(bites_table), biters)
  NB[index] =bites_table
  NB
}

n.exposed <- bites(cases$ID, exposures$Biter.ID) # people exposed per rabid ANIMAL (because cases are ALL spp!)
test = bites2(cases$ID, exposures$Biter.ID)
sum(n.exposed); sum(test)

##################################################################
n.exp <- bites(subset(cases, Species == "Domestic dog")$ID, exposures$Biter.ID) # exposures per rabid dog
n.exposed <- bites(cases$ID, exposures$Biter.ID) # number of people exposed per rabid ANIMAL (because cases are ALL spp!)
n.bitten <- bites2(animals$ID, !is.na(animals$Biter.ID)) # Number of dogs BITTEN per RABID DOG

max(animals$Biter.ID)
max(animals$ID)
test1 = hist(animals$Biter.ID, breaks = -1:607, plot=F)$counts[-1] 
test2 = rep(NA, length(rabid$ID))

n.rabid <- bites(rabid$ID, rabid$Biter.ID) # Number of dogs INFECTED (i.e. secondary cases) per RABID DOG

hist(n.exposed, -1:15) # THE DISTRIBUTION OF EXPOSED PEOPLE! 
hist(n.bitten, -1:15) # THE DISTRIBUTION OF BITTEN ANIMALS! 
hist(n.rabid, -1:15) # THE DISTRIBUTION OF SECONDARY CASES (RABID ANIMALS) - I think this will be quite biased but we should look at the transmission trees!

# Export animal cases (with IDs) and persons bitten:
rabid_exposed = data.frame(id_case = rabid$ID, species = rabid$Species, exposed = n.exposed)
write.csv(rabid_exposed, "Output/ID_exposed.csv", row.names = FALSE) # to be used for Fig 4

# Calculate rabid dogs under case detection probability
est_rabid_dogs 
sum(rabid_dogs)
sum(est_rabid_dogs)
cases1 <- subset(cases, Species == "Domestic dog" & Year < 2016)
cases2 <- subset(cases, Species == "Domestic dog" & Year > 2015)
exposures1 <- subset(exposures, Year.bitten < 2016)
exposures2 <- subset(exposures, Year.bitten > 2015)
n.exposed1 <- bites(cases1$ID, exposures1$Biter.ID) 
n.exposed2 <- bites(cases2$ID, exposures2$Biter.ID) 
nonbiters <- est_rabid_dogs - rabid_dogs; sum(nonbiters)
biters <- hist(n.exposed, -1:15, plot=FALSE)$counts
biters1 <- hist(n.exposed1, -1:15, plot=FALSE)$counts
biters2 <- hist(n.exposed2, -1:15, plot=FALSE)$counts
biters[1] <- biters[1] + round(sum(nonbiters))
biters1[1] <- biters1[1] + round(nonbiters[1])
biters2[1] <- biters2[1] + round(nonbiters[2])

# - Estimate parameter values for mu & k (consider pDetect for cases!)
negbinom.params <- fitdistr(rep(0:15, biters),"negative binomial", method = "SANN")$estimate
NBpars <- fitdist(rep(0:15, biters),"nbinom")$estimate
NB1 <- fitdist(rep(0:15, biters1),"nbinom")$estimate
NB2 <- fitdist(rep(0:15, biters2),"nbinom")$estimate

# CIs on negative binomial parameters
statistic1 <- function(x, inds) {fitdist(x[inds],"nbinom")$estimate[1]}
statistic2 <- function(x, inds) {fitdist(x[inds],"nbinom")$estimate[2]}
bs1 <- boot(rep(0:15, biters), statistic1, R = 4000)
bs2 <- boot(rep(0:15, biters), statistic2, R = 4000)
print(boot.ci(bs1, conf=0.95, type="bca")$bca[c(4,5)])
print(boot.ci(bs2, conf=0.95, type="bca")$bca[c(4,5)])

# save parameter estimates
NBparams = data.frame(
  mu = NBpars["mu"], size = NBpars["size"],
  mu1 = NB1["mu"], size1 = NB1["size"],
  mu2 = NB2["mu"], size2 = NB2["size"]
  )
write.csv(NBparams, "output/NBparams.csv", row.names = FALSE)

