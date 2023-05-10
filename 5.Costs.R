#' Calculate number of dogs vaccinated annually and PEP administered (for rabies and healthy dog bites).
#' Estimate cost per dog vaccinated and approximate annual dog vaccination campaign
#' and cost per PEP and annual spend on PEP
rm(list=ls())

## load libraries
library(tidyverse)
library(lubridate)
library(Hmisc)
options(stringsAsFactors = FALSE)

#---------------------------------------------------------
# Import data
vacc <- read.csv("data/Vax_Pemba_2010_to_2019.csv") # vaccination data
vacc$ID_shehia <- paste(vacc$District, vacc$Shehia)
human <- readRDS(file = "output/humans_deid.rda") # human exposures & bite patients
params <- read.csv("data/bio_data.csv") # parameters for protective efficacy of vaccine from Changalucha et al; WHO modelling consortium

# Sort vaccination dates
vacc$date <- as.Date(vacc$vacc_date, "%d/%m/%Y")
vacc$mth <- (as.POSIXlt(vacc$date)$year-110)*12 + (as.POSIXlt(vacc$date)$mon+1)
vacc$year <- year(as.Date(vacc$vacc_date,format="%d/%m/%Y"))

# Annual dog vaccination data
vax.year <- vacc %>%
  mutate(year = factor(year, levels = 2010:2021)) %>%
  group_by(year, .drop=FALSE)  %>%  # group by year
  dplyr::summarize(vax = sum(dogs_vaccinated, na.rm=TRUE))
vax.year

# human bites and rabies exposures,  and patients (that went for PEP) for bites by healthy animals
bites_y <- human %>%
  mutate(Year = factor(Year.bitten, levels = 2010:2021)) %>%
  group_by(Year, .drop=FALSE)  %>%
  dplyr::summarize(
    all_bites = length(which(Rabid=="Yes" | Rabid == "Unknown"|Rabid=="No")),
    bites_suspects = length(which(Rabid=="Yes" | Rabid == "Unknown")),
    bites_healthy = length(which(Rabid=="No")),
    deaths = length(which(Patient.outcome == "Died")),
    PEP1 = length(which(PEP.1=="true")), # PEP doses 1, 2 and 3
    PEP2 = length(which(PEP.2=="true")),
    PEP3 = length(which(PEP.3=="true")),
    PEP4 = length(which(PEP.4=="true")),
    PEP_rabid_complete = length(which(Rabid=="Yes"|Rabid == "Unknown" & PEP.3 == "true")),
    PEP_rabid_complete2 = length(which(Rabid=="Yes" & PEP.3 == "true")),
    PEP_rabid_started = length(which(Rabid=="Yes"|Rabid == "Unknown" & PEP.1 == "true")),
    PEP_rabid_started2 = length(which(Rabid=="Yes" & PEP.1 == "true")),
    PEP_healthy_complete = length(which(Rabid=="No" & PEP.3 == "true")),
    exposures_survived = length(which(Rabid=="Yes"|Rabid=="Unknown" & Patient.outcome != "Died")))
bites_y

#' _HOW MANY RABIES DEATHS PREVENTED from Contact tracing data?_
pRabies <- params$p_rabies_transmission
# People who started vs completed PEP
sum(bites_y$PEP_rabid_complete)
sum(bites_y$PEP1); sum(bites_y$PEP3)
# Annually
bites_y$PEP3/bites_y$PEP1
bites_y$PEP1; bites_y$PEP3 

sum(bites_y$PEP_rabid_complete) * pRabies # deaths prevented by PEP
deaths_averted_PEP = round(binconf(sum(bites_y$PEP_rabid_complete) * pRabies,
                               sum(bites_y$PEP_rabid_complete))*sum(bites_y$PEP_rabid_complete))
deaths_averted_PEP2 = round(binconf(sum(bites_y$PEP_rabid_complete2) * pRabies,
                                   sum(bites_y$PEP_rabid_complete2))*sum(bites_y$PEP_rabid_complete2))
poss_deaths_averted_PEP = round(binconf(sum(bites_y$PEP_rabid_started) * pRabies,
                                   sum(bites_y$PEP_rabid_started))*sum(bites_y$PEP_rabid_started))
poss_deaths_averted_PEP2 = round(binconf(sum(bites_y$PEP_rabid_started2) * pRabies,
                                        sum(bites_y$PEP_rabid_started2))*sum(bites_y$PEP_rabid_started2))

#' _COST ANALYSES_
usd = 2296 # USD conversion to Tsh

#'*Control measures in humans*
#'depends on the PEP regimen (doses), vial cost, wound cleaning & health workers time etc
overhead = 25000/ usd # Consultation - covers cost of building, hospital infrastructure, consultation fees etc from the National Health Insurance schemes
vial_cost = 10.98 # usd - or ~30000 TZS a price of vials paid by bite patient
syringes_needle = (123 / usd) 
syringes_needle
avg_yearly_salary =  1293100*12/ usd # avg RCT nurse salary (public health nurse / nursing assistant): Tsh 1,293,100 per mnth = 15,517,200 Tsh or $6758.362 annually
working_mins_per_yr = 40*5*8*60 # 40 weeks per year, 5 days per week, 8 hours per day * 60 mins per hour
HW_timecost = (avg_yearly_salary/ working_mins_per_yr)* 5 # Assuming 5 minutes per injection
HW_timecost # $6758.362 * 40/ (1600*60) = $2.82 per patient assuming they receive the full course

# Type of PEP
IM_course = (vial_cost + syringes_needle + HW_timecost)*4 + overhead; IM_course  # up to 4 injections for IM PEP course
ID_course = vial_cost + (syringes_needle + HW_timecost)*8 + overhead; ID_course # up to 8 injections for ID PEP course
IM_inj = vial_cost + syringes_needle + HW_timecost; IM_inj # use 1 vial per IM injection
ID_inj = vial_cost/2 + (syringes_needle + HW_timecost)*2; ID_inj # use 1/2 vial per ID visit

# Indirect costs (not included in manuscript)
transport_victims = 3000 # assuming the victim is always accompanied by an adult
visits = 3 # time lost by victims during vaccine series (hrs)
avgdaily_wage_victim = 5000
pep_indirect = (transport_victims + avgdaily_wage_victim) * visits/usd; pep_indirect # $12 losses direct to patient 

# TOTAL COSTS
pep_IM_y = (bites_y$PEP1 + bites_y$PEP2 + bites_y$PEP3 + bites_y$PEP4)  * IM_inj + (bites_y$PEP1*overhead)
pep_ID_y = (bites_y$PEP1 + bites_y$PEP2 + bites_y$PEP3 + bites_y$PEP4)  * ID_inj + (bites_y$PEP1*overhead)
pep_IM_y; pep_ID_y
sum(pep_IM_y); sum(pep_ID_y) # Total cost of PEP if IM or ID
mean(pep_ID_y); mean(pep_IM_y) # Average annual cost of PEP if IM or ID

# Proportion of bites given IM vs ID vaccination
PEP_route <- table(human$PEP.administration)
IM_pc <- PEP_route["IM"]/sum(PEP_route[c("IM","IV")])  # 18% given as IM
govt_cost <- (sum(pep_IM_y) * IM_pc) + (sum(pep_ID_y)*(1-IM_pc)) # Total cost of PEP given IM & ID regimens
annual_gov_cost <- mean(pep_IM_y * IM_pc) + mean(pep_ID_y *(1-IM_pc)) # Average annual cost of PEP given IM or ID regimens
govt_cost; annual_gov_cost

#' _Annual government spending on PEP in the absence of rabies_
yrs <- as.numeric(as.character(bites_y$Year))
rabies_free <- which(yrs>2018)
rabies_free_presentations <- sum(bites_y$PEP_healthy_complete[rabies_free])/length(rabies_free)
rabies_free_cost_IM <- sum(bites_y$PEP_healthy_complete[rabies_free]) * IM_course/length(rabies_free) 
rabies_free_cost_ID <- sum(bites_y$PEP_healthy_complete[rabies_free]) * ID_course/length(rabies_free)
rabies_free_PEP <- ((rabies_free_cost_ID * (1-IM_pc)) + (rabies_free_cost_IM * IM_pc))
rabies_free_PEP # annual cost of precautionary PEP after rabies eliminated

#' _what is the cost-effectiveness of PEP per death averted ?_
cost_per_death_avert <- govt_cost/deaths_averted_PEP #
cost_per_death_avert # With 95% CIs 

#' _Costs for rabies control in dogs_
#'The costs of mass vaccination include costs of the vaccine, consumables such as needles, syringes, etc.
#'vaccinators, the information campaign/sensitization, and time from dog owners
dose_vaccine = 1500/usd #' _****CHECK****_ quite high - used to be $0.20 per dose  # Tsh per dose (ministry price guidance) 
syringes_needle = 100/usd
collar = 250/usd
dogs_vacc = vax.year$vax # Number of dogs vaccinated each year
per_dog_cost <- dose_vaccine + syringes_needle + collar

shehias = 120
shehias_per_day = 4 
adverts = 17100/usd * shehias/ shehias_per_day # printing and distribution of the leaflets and posters, broadcast of the radio advertisements/sensitization - per CP
stationary = 6000/usd * shehias # This includes mark pens, flipcharts, pens, pins/glue for a vaccinator
fuel = 17100/usd * shehias # per day for the vaccinator

# A vaccination team for each station is made up of three people (two LFOs, and 1 CAHW)
vaccinator_perdiem = 30000/usd # vaccinator cost (~£14/ day!)
vaccinator_daily_salary = 55572.97/usd # $24 vs $33 for a nurse
CAHW_allowance = 5000/usd
daily_team_cost = 2*(vaccinator_perdiem + vaccinator_daily_salary) + CAHW_allowance
shehias_vacc = ifelse(dogs_vacc < 1000, shehias/4, shehias) 
campaign = ((daily_team_cost * shehias_vacc) + adverts + stationary + fuel) * (dogs_vacc > 0) 
# HW & LFO costs were v different between mainland Tz & Pemba - limitation is used Tz mainland

# campaign costs
n <-  sum(vax.year$vax>0); n# campaigns
vax = per_dog_cost * dogs_vacc
vax_max = per_dog_cost * max(dogs_vacc) # The best coverage campaign
vax_mean = per_dog_cost * sum(dogs_vacc)/n # mean per campaign

annual_MDV <- vax + campaign
annual_MDV_mean <- sum(annual_MDV) / sum(vax>0); annual_MDV_mean # mean
island_MDV_avg <- sum(annual_MDV[annual_MDV>5000]) / sum(annual_MDV>5000); island_MDV_avg # mean
island_MDV_range <- range(annual_MDV[annual_MDV>5000]); island_MDV_range 

total_MDV <- sum(annual_MDV)
range(annual_MDV[which(vax>0)]) # 4233 - 13145 USD per campaign depending on numbers vaccinated PLUS VACCINATOR COSTS
cost_per_dog_vax = total_MDV/ sum(dogs_vacc); cost_per_dog_vax 
range(annual_MDV[annual_MDV>0]/dogs_vacc[annual_MDV>0]) # range of dog vaccination costs

# annual cost of precautionary PEP after rabies eliminated:
rabies_free_PEP # When dog rabies eliminated.... still get dog bite presentations from healthy dog bites...
# Govts pay this residual PEP cost, PLUS annual MDV  - TOGETHER = One Health preventative costs....
# One Health approach - SENSE CHECK
deaths = 3
time_horizon <- 10
OH_deaths_averted <- deaths_averted_PEP + (3 * time_horizon) # deaths averted per year by PEP plus additional deaths spared through dog vaccination
annual_precautionary_PEP <- ID_course/usd * 37
OH_costs <- (annual_MDV_mean + annual_precautionary_PEP) * time_horizon
OH_deaths_averted
OH_costs/OH_deaths_averted