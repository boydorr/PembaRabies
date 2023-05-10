# RUN DECISION TREE MODEL
# Using the contact tracing, estimate the combined costs of PEP & MDV over the decade
# The probabilities of seeking & completing PEP, rabid dog biting etc

rm(list = ls())
library(tidyverse)
library(scales)
library(reshape)

#' * Discounting *
#' * PEP costs - IM vs ID options *
#' * Scenarios - SQ (no MDV + charge for PEP), free PEP (no MDV), One Health (free PEP + MDV) *
#' * Vary - PEP regimen (IM Essen, ID UTRC, ID 1-week)*
#' * Run count - 1000 *

# Pemba health econ model
source("R/decision_tree.R")
source("R/HelperFun.R") # Simulate timeseries from pDetect & Scenario

# PARAMETERS
# from literature
PEP_protect_pars <- read.csv("data/bio_data.csv") # Protective efficacy of vaccine (from WHO modeling consortium 2018)
pRabies <- PEP_protect_pars$p_rabies_transmission # Probability of developing rabies following exposure if no PEP
pPrev_imperfect <- PEP_protect_pars$p_prevent_given_imperfect # Probability of developing rabies following exposure if late/ incomplete PEP

# from Pemba contact tracing:
detectQs <- read.csv("Output/DetectQs.csv") # Detection probabilities
q50 <- which(detectQs$X== "50%")
pDetect <- detectQs$TDoverall[q50] # subset(detectQs, X=="50")$TDoverall

biting_pars <- read.csv("Output/NBparams.csv") # Rabid dog biting 
PEP_seek_pars <- read.csv("Output/PEP_params_Pemba.csv") # PEP seeking & completion
PEP_course <- read.csv("Output/PEPcourse.csv") # PEP costs
healthy <- read.csv("Output/healthy_bites.csv") # Incidence of bites by healthy dogs
healthy_seek_ratio <- healthy$BI2/healthy$BI1
humanpop <- 472958 # from 2012 census
HDR <- 100 # Estimate

# Time series from Pemba
bites <- read.csv("Output/bite_summary.csv")
cases <- read.csv("Output/cases_y.csv")

# Prepare dataframe to sample from
pDetect_y <- rep(detectQs$TD1[q50], 11) # initialize detection probability 
pDetect_y[which(cases$Year>2015)] <- detectQs$TD2[q50] # detection probability over time
df = data.frame(cases = cases$dogs, pDet = pDetect_y, VR = cases$VaxRound)

# Prepare endemic vs MDV scenarios
endemic = df[which(df$VR<2),] # If less than 2 vaccination rounds
mdv <- df; mdv[] <- 0
mdv1 <- mdv2 <- mdv # create two MDV variants
mdv1[1:6,] <- df[1:6,] # from the first Pemba period
mdv2[1:5,] <- df[7:11,] # from the Pemba outbreak

# Run 1000 iterations for all scenarios
N = 1000

# STATUS QUO - endemic & charge for PEP, discount at 3%, low Pstart and Pcomplete, IM vaccination
out_SQ <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "endemic", PEPpolicy = "charge", 
                              mu = biting_pars$mu, k = biting_pars$size,
                              pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                              pDeath = pRabies, pPrevent = pPrev_imperfect, 
                              full_cost = PEP_course$IM_full, partial_cost = PEP_course$IM_incomplete, campaign_cost = 0)
# SQ variation - use ID regimen (without changing to WHO recommendations of 1-week IPC)
out_SQ_ID <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "endemic", PEPpolicy = "charge",  
                                 mu = biting_pars$mu, k = biting_pars$size,
                                 pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                                 pDeath = pRabies, pPrevent = pPrev_imperfect, 
                                 full_cost = PEP_course$ID_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 0)
# SQ variation - no discounting
out_SQ_0disc <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0, setting = "endemic", PEPpolicy = "charge", 
                                   mu = biting_pars$mu, k = biting_pars$size,
                                   pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                                   pDeath = pRabies, pPrevent = pPrev_imperfect,
                                   full_cost = PEP_course$IM_full, partial_cost = PEP_course$IM_incomplete, campaign_cost = 0)


# FREE PEP - endemic rabies, discount at 3%, high Pstart and Pcomplete, IM vaccination
out_freePEP <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "endemic", PEPpolicy = "free", 
                                mu = biting_pars$mu, k = biting_pars$size,
                                pStart = PEP_seek_pars$pStart[2], pComplete = PEP_seek_pars$pComplete[2], 
                                pDeath = pRabies, pPrevent = pPrev_imperfect, 
                                full_cost = PEP_course$IPC_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 0)
# free PEP variation - no discounting
out_free_0disc <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "endemic", PEPpolicy = "free", 
                      mu = biting_pars$mu, k = biting_pars$size,
                      pStart = PEP_seek_pars$pStart[2], pComplete = PEP_seek_pars$pComplete[2], 
                      pDeath = pRabies, pPrevent = pPrev_imperfect, 
                      full_cost = PEP_course$IPC_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 0)


# MDV - charge for PEP, discount at 3%, low Pstart and Pcomplete, IM vaccination
out_MDV <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "mdv", PEPpolicy = "charge", 
                    mu = biting_pars$mu, k = biting_pars$size,
                    pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                    pDeath = pRabies, pPrevent = pPrev_imperfect, 
                    full_cost = PEP_course$IM_full, partial_cost = PEP_course$IM_incomplete, campaign_cost = 12323.7)
# MDV variation - no discounting
out_MDV_0disc <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0, setting = "mdv", PEPpolicy = "charge", 
                               mu = biting_pars$mu, k = biting_pars$size,
                               pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                               pDeath = pRabies, pPrevent = pPrev_imperfect, 
                               full_cost = PEP_course$IM_full, partial_cost = PEP_course$IM_incomplete, campaign_cost = 12323.7)
# MDV variation - ID PEP
out_MDV_ID <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "mdv", PEPpolicy = "charge", 
                                     mu = biting_pars$mu, k = biting_pars$size,
                                     pStart = PEP_seek_pars$pStart[1], pComplete = PEP_seek_pars$pComplete[1], 
                                     pDeath = pRabies, pPrevent = pPrev_imperfect, 
                                     full_cost = PEP_course$ID_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 12323.7)


# ONE HEALTH - MDV & free PEP, discount at 3%, high Pstart and Pcomplete, ID vaccination
out_OH <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0.03, setting = "mdv", PEPpolicy = "free", 
                    mu = biting_pars$mu, k = biting_pars$size,
                    pStart = PEP_seek_pars$pStart[2], pComplete = PEP_seek_pars$pComplete[2], 
                    pDeath = pRabies, pPrevent = pPrev_imperfect, 
                    full_cost = PEP_course$IPC_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 12323.7)
# OH variation - no discounting
out_OH_0disc <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0, setting = "mdv", PEPpolicy = "free", 
                              mu = biting_pars$mu, k = biting_pars$size,
                              pStart = PEP_seek_pars$pStart[2], pComplete = PEP_seek_pars$pComplete[2], 
                              pDeath = pRabies, pPrevent = pPrev_imperfect, 
                              full_cost = PEP_course$IPC_full, partial_cost = PEP_course$ID_incomplete, campaign_cost = 12323.7)

# BASELINE - endemic rabies WITHOUT use of PEP 
out_baseline <- decision_tree_ndraw(ndraw = N, horizon=10, discount=0, setting = "endemic", PEPpolicy = "charge", 
                                    mu = biting_pars$mu, k = biting_pars$size,
                                    pStart = 0, pComplete = 0, 
                                    pDeath = pRabies, pPrevent = pPrev_imperfect, 
                                    full_cost = 0, partial_cost = 0, campaign_cost = 0)

# Check number of healthy bite patients (under SQ i.e. charged for PEP and free PEP)
quantile(out_SQ$healthy_bite_patients, c(0.025, 0.5, 0.975)) # 2010-2014 had 6-12 per year
quantile(out_freePEP$healthy_bite_patients, c(0.025, 0.5, 0.975)) # 2016+ had ~37 per year

# SUMMARIZE OUTCOMES: deaths; DA, PEPcost; MDVcost; Total Cost (and also look at pattern over time)
sum_base <- sum_horizon(out_baseline, "base", 0.03, "none") # baseline - no PEP or MDV
sum_SQ <- sum_horizon(out_SQ, "SQ", 0.03, "IM") # SQ & variants
sum_SQ_ID <- sum_horizon(out_SQ_ID, "SQ", 0.03, "ID")
sum_SQ_0disc <- sum_horizon(out_SQ_0disc, "SQ", 0, "IM")
sum_freePEP <- sum_horizon(out_freePEP, "free PEP", 0.03, "ID") # free PEP & variants
sum_free_0disc <- sum_horizon(out_free_0disc, "free PEP", 0, "ID")
sum_MDV <- sum_horizon(out_MDV, "MDV", 0.03, "IM") # MDV & variants
sum_MDV_0disc <- sum_horizon(out_MDV_0disc, "MDV", 0, "IM")
sum_MDV_ID <- sum_horizon(out_MDV_ID, "MDV", 0.03, "ID") 
sum_OH <- sum_horizon(out_OH, "OH", 0.03, "ID") # OH & variants
sum_OH_0disc <- sum_horizon(out_OH_0disc, "OH", 0, "ID")

summaries <- rbind.data.frame(
  horizon_CEA(sum_base, sum_base),
  horizon_CEA(sum_base, sum_SQ),
  horizon_CEA(sum_base, sum_SQ_ID),
  horizon_CEA(sum_base, sum_SQ_0disc),
  horizon_CEA(sum_base, sum_freePEP),
  horizon_CEA(sum_base, sum_free_0disc),
  horizon_CEA(sum_base, sum_MDV),
  horizon_CEA(sum_base, sum_MDV_0disc),
  horizon_CEA(sum_base, sum_MDV_ID),
  horizon_CEA(sum_base, sum_OH),
  horizon_CEA(sum_base, sum_OH_0disc))
write.csv(summaries, "Output/CEA_summaries.csv", row.names = TRUE)

setDT(summaries)
death_summ <- melt(summaries, id.vars = c("scenario","CIs","discount","PEP"), measure.vars = c("deaths"))
death_summ$CIs <- as.factor(death_summ$CI)
death_scen_sum <- cast(death_summ, scenario+variable+discount+PEP ~ CIs, value = "value")

cost_perDA_summ <- melt(summaries, id.vars = c("scenario","CIs","discount","PEP"), measure.vars = c("cost_per_deaths_avert"))
cost_perDA_summ$CIs <- as.factor(cost_perDA_summ$CI)
CEA_scen_sum <- cast(cost_perDA_summ, scenario+variable+discount+PEP ~ CIs, value = "value")

scen_sum <- rbind.data.frame(death_scen_sum, CEA_scen_sum)
names(scen_sum)[5:7] <- c("lci", "med","uci")


# Plot Deaths & Cost per DA for the different scenarios: SQ (IM vs ID), free PEP (IM vs ID), MDV, OH
# Summarize over time horizon
# SQ
SQ_mean <- year_avg(model_output = out_SQ, scenario = "SQ")
SQ_median <- year_summary(model_output = out_SQ, q = 0.5, scenario = "SQ")
SQ_lci <- year_summary(model_output = out_SQ, q = 0.0275, scenario = "SQ")
SQ_uci <- year_summary(model_output = out_SQ, q = 0.975, scenario = "SQ")

# free PEP
freePEP_mean <- year_avg(model_output = out_freePEP, scenario = "free PEP")
freePEP_median <- year_summary(model_output = out_freePEP, q = 0.5, scenario = "free PEP")
freePEP_lci <- year_summary(model_output = out_freePEP, q = 0.0275, scenario = "free PEP")
freePEP_uci <- year_summary(model_output = out_freePEP, q = 0.975, scenario = "free PEP")

# One Health
OH_mean <- year_avg(model_output = out_OH, scenario = "OH")
OH_median <- year_summary(model_output = out_OH, q = 0.5, scenario = "OH")
OH_lci <- year_summary(model_output = out_OH, q = 0.0275, scenario = "OH")
OH_uci <- year_summary(model_output = out_OH, q = 0.975, scenario = "OH")

ts_summary <- rbind.data.frame(SQ_mean, SQ_median, SQ_lci, SQ_uci, 
                               freePEP_mean, freePEP_median, freePEP_lci, freePEP_uci, 
                               OH_mean, OH_median, OH_lci, OH_uci)
setDT(ts_summary)
exp_long <- melt(ts_summary, id.vars = c("years","scenario","metric"), measure.vars = c("exposures"))
exp_scen <- cast(exp_long, scenario+years+variable ~ metric, value = "value")
deaths_long <- melt(ts_summary, id.vars = c("years","scenario","metric"), measure.vars = c("deaths"))
deaths_scen <- cast(deaths_long, scenario+years+variable ~ metric, value = "value")
DA_long <- melt(ts_summary, id.vars = c("years","scenario","metric"), measure.vars = c("deaths_averted"))
DA_scen <- cast(DA_long, scenario+years+variable ~ metric, value = "value")
# PEPcost_long <- melt(ts_summary, id.vars = c("years","scenario","metric"), measure.vars = c("PEP_cost"))
# PEPcost_scen <- cast(PEPcost_long, scenario+years+variable ~ metric, value = "value")
# tsCIs <- rbind.data.frame(deaths_scen, exp_scen, PEPcost_scen) # DA_scen, 
tsCIs <- rbind.data.frame(deaths_scen, exp_scen) # DA_scen, 
names(tsCIs)[5:7] <- names(deaths_scen)[5:7] <- names(exp_scen)[5:7] <- c("med","lci","uci")

facet_labels <- c(deaths="deaths", exposures="rabies exposures", PEP_cost = "PEP costs")
palette <- c("firebrick2", "#E69F00", "#0072B2")

scen_ts <- ggplot(data=tsCIs) + geom_line(aes(x=years, y=avg, col=scenario), lwd=1) +
  facet_wrap(~variable, ncol=3, strip.position ="top", scales="free_y", labeller = as_labeller(facet_labels)) +
  theme_bw(base_size = 10) +
  labs(y="Values", x="Year") +
  scale_y_continuous(limits = c(0,NA)) + scale_x_continuous(limits = c(0,10), breaks=1:10) +
  geom_ribbon(aes(x=years, ymin=lci, ymax=uci, fill=scenario), alpha=0.15) +
  scale_colour_manual(values=palette, breaks=c("SQ", "free PEP", "OH"), 
                      labels=c("Status Quo", "free PEP", "One Health")) +
  scale_fill_manual(values=palette, breaks=c("SQ", "free PEP", "OH"), 
                    labels=c("Status Quo", "free PEP", "One Health")) +
  theme(legend.title=element_blank(),
        legend.position="right",
        strip.background = element_blank(),
        strip.text.y = element_text(size=10),
        plot.margin = unit(c(.1,0.5,0.1,0.5), "cm"))
scen_ts
ggsave("figures/scenario_timeseries.jpeg", height = 5, width = 10)

# Plot scenario arrangements
scen_sum$scenario <- factor(scen_sum$scenario, levels = c("SQ", "free PEP", "MDV", "OH"))
scen_sum <- subset(scen_sum, discount == 0.03)
all_palette <- c("firebrick2", "#E69F00", "deepskyblue4", "#0072B2")
outcome_names <- as_labeller(c("deaths"="deaths","cost_per_DA_OH"="cost per death averted"))
facet_labeller <- function(variable,value){return(outcome_names[value])}

scen_summary <-ggplot(data=scen_sum, aes(x=scenario, y=med, group=variable)) +
  theme_classic() +
  geom_point(aes(colour=scenario, shape=PEP, fill= PEP), position=position_dodge2(.25)) +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour=scenario, lty=PEP, width=0.3), position=position_dodge2(.5)) +
  facet_wrap( ~ variable, scales="free_y", labeller = outcome_names) + 
  labs(colour="Scenario", y="Value", x="") +
  theme(plot.title = element_text(size=10, face="bold", hjust = 0.5)) +
  theme(legend.title=element_blank()) + theme_classic(base_size = 12) +
  scale_y_continuous(labels=comma, limits=c(0,NA)) +
  scale_colour_manual(values=all_palette,
                      labels=c("Status Quo", "Free PEP", "MDV", "One Health")) +
  scale_x_discrete(breaks=c("SQ","free PEP","MDV","OH"),
                   labels=c("Status quo","free PEP","Dog Vaccination","One Health")) +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=10), legend.position="right") + guides(color="none")
scen_summary
ggsave("figures/scenario_summary.pdf", height = 4, width = 8)

##########################################################
# simplified plot!
scen_plot <- subset(scen_sum, discount == 0.03 & 
                      (scenario == "SQ"|scenario == "free PEP"|scenario == "OH"))
scen_plot <- scen_plot[-which(scen_plot$PEP == "ID" & scen_plot$scenario == "SQ"),]

scen_plot$scenario <- factor(scen_plot$scenario, levels = c("SQ", "free PEP", "OH"))

outcome_names <- as_labeller(c("deaths"="deaths","cost_per_deaths_avert"="cost per death averted"))
facet_labeller <- function(variable,value){return(outcome_names[value])}

scen_summary_simple <-ggplot(data=scen_plot, aes(x=scenario, y=med, group=variable)) +
  theme_classic() +
  geom_point(aes(colour=scenario, shape=PEP, fill= PEP), position=position_dodge2(.25)) +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour=scenario, lty=PEP, width=0.3), position=position_dodge2(.5)) +
  facet_wrap( ~ variable, scales="free_y", labeller = outcome_names) + 
  labs(colour="Scenario", y="Value", x="") +
  theme(plot.title = element_text(size=10, face="bold", hjust = 0.5)) +
  theme(legend.title=element_blank()) + theme_classic(base_size = 10) +
  scale_y_continuous(labels=comma, limits=c(0,NA)) +
  scale_colour_manual(values=palette,
                      labels=c("Status Quo", "Free PEP", "One Health")) +
  scale_x_discrete(breaks=c("SQ","free PEP","OH"),
                   labels=c("Status quo","free PEP","One Health")) +
  theme(legend.title=element_blank(), 
        strip.background = element_blank(),
        strip.text.y = element_text(size=10),
        axis.text.x  = element_text(angle=45, hjust=1, size=10), 
        plot.margin = unit(c(.1,0.5,0.1,0.5), "cm"), 
        legend.position="right") + guides(color="none")
scen_summary_simple
        
# Combine health econ plots
plot_CEA <-
  (scen_ts / scen_summary_simple) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
plot_CEA 
ggsave("figures/CEA_summaries.jpeg", plot_CEA , height = 6, width = 6)


# Plot the cost-effectiveness plane for the difference in costs 
# y = difference in costs, x = difference in costs
ICER_median <- subset(summaries, iter == 500.500 & discount == 0.03)
range(ICER_median$cost_diff)
range(ICER_median$death_diff)
yaxis=120000
xaxis=100
ICER_median$scenario <- factor(ICER_median$scenario, levels=c("SQ", "free PEP", "MDV", "OH")) 
ICER_median$Route <- factor(c("IM","ID","ID","IM","ID", "ID"))

ICER_plane <- ggplot(data=ICER_median, aes(x=death_diff, y=cost_diff)) +
  geom_vline(xintercept = 5, col="gray70") +
  geom_hline(yintercept = 0, col="gray70") +
  geom_point(aes(col=scenario, shape=Route), size=3) +
  theme_classic() +
  scale_x_reverse(labels=comma, limits=c(xaxis, -xaxis), breaks=c(100,50,0,-50,-100)) +
  scale_y_continuous(labels= comma, limits=c(-yaxis, yaxis)) +
  labs(x="Difference in rabies deaths", y="Difference in cost (USD)") +
  scale_colour_manual(values=all_palette) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12),
        legend.title = element_text(size=14), legend.text = element_text(size=12))
ICER_plane
ggsave("figures/ICER_plane.jpeg", height = 5, width = 6)




