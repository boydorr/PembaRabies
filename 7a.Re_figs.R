# Make figure with Re time series & vs. cov----

# pkgs ----
library(data.table)
library(tidyverse)
library(magrittr)
library(cowplot)
library(patchwork)
library(lubridate)

# load in outputs ----
ttrees <- fread("Output/trees/trees_all.gz")
links_consensus <- fread("Output/trees/links_consensus_consistent.csv") 
scenarios <- fread("Output/trees/scenarios.csv")
case_dt <- fread("Output/trees/case_dt_cleaned.csv")
coverage <- read.csv("Output/dogVCWaningPembaWard.csv", header = FALSE)
case_exposures <- read.csv("Output/ID_exposed.csv") 
case_date_locs <- case_dt[, .(id_case = id, symptoms_started, district, ward, 
                              village, utm_easting, 
                              utm_northing)]
case_date_locs <- case_date_locs[case_exposures, on = c("id_case")] # join with exposure data

# clean up labels and pull in colors ----
scenarios[, c("data_used", 
              "prune_type") := .(ifelse(!use_gen, "Epi data only", "Epi & genetic data"), 
                                 fcase(!prune, "Unpruned", 
                                       prune & dist_cutoff == 0.99, 
                                       "Pruned by time & distance", 
                                       prune & dist_cutoff == 1, 
                                       "Pruned by time"))]
scenarios <- scenarios[, c("scenario", "data_used", "prune_type", "use_gen")]
scenarios$prune_type <- factor(scenarios$prune_type, 
                               levels = c("Unpruned", "Pruned by time", 
                                          "Pruned by time & distance"))
scenarios$data_used <- factor(scenarios$data_used, 
                              levels = c("Epi data only", "Epi & genetic data"))
data_used_col <- c("#024B79", "#FFAD48")

# Look at individual Re estimates (across scenarios) ----
mean_re <- ttrees[, .(reff = mean(reff)), by = c("id_case", "scenario")]

# join with links consensus to get the other data
links_re <- links_consensus[mean_re, on = c("id_case", "scenario")]

# join with case data
links_re <- case_date_locs[links_re, on = c("id_case")]

# time series: moving average of cases forward looking over 3 mos.
grps <- expand.grid(month = seq.Date(ymd("2010-01-01"), ymd("2019-01-01"), by = "month"), 
                    join = -2:3, scenario = 1:6)
grps$month_group <- floor_date(as_date(grps$month + dmonths(grps$join)), unit = "months")

links_re %>%
  mutate(month = floor_date(ymd(symptoms_started), unit = "month")) %>%
  right_join(grps) %>%
  group_by(month_group, scenario) %>%
  dplyr::summarize(mean_re = mean(reff, na.rm = TRUE), 
            upper = quantile(reff, 0.975, na.rm = TRUE), 
            lower = quantile(reff, 0.025, na.rm = TRUE), nobs = n()) -> re_ts

# looks pretty much the same across all the different algs.
ggplot() +
  geom_line(data = re_ts, aes(x = month_group, y = mean_re)) +
  facet_wrap(~scenario)

# join up with the scenarios
links_re <- scenarios[links_re, on = "scenario"]
re_ts <- left_join(re_ts, scenarios, by = "scenario")

# Filter to the chosen scenario & also plot cases colored by lineage ---
# to make a final decision, but for now going with 
# using genetic data only to diff chains
links_re <- links_re[prune_type == "Unpruned" & use_gen == TRUE]
re_ts <- filter(re_ts, prune_type == "Unpruned", use_gen == TRUE)

# lineage colors
cols_chains <- c("grey50",
                 "#ebac23", "#b80058", "#006e00", "#00bbad",
                 "#ff9287", "#5954d6", "#d163e6", "#b24502", "#00c6f8",
                 "#878500", "#008cf9")
names(cols_chains) <- c(0:7)
cols_chains <- cols_chains[2:8] # none that aren't assigned to a chain

plot_ts <-
  ggplot(data = links_re) +
  geom_ribbon(
    data = re_ts,
    aes(x = month_group, ymin = lower, ymax = upper),
    color = "NA", fill = "grey50", alpha = 0.25
  ) +
  geom_point(aes(x = symptoms_started, y = reff, fill = factor(lineage_chain)),
    size = 2.5, alpha = 0.5, shape = 21, color = "NA"
  ) +
  geom_point(aes(x = symptoms_started, y = reff, color = factor(lineage_chain)),
    size = 2.5, alpha = 1, shape = 21, fill = "NA"
  ) +
  geom_line(
    data = re_ts,
    aes(x = month_group, y = mean_re), size = 1.2,
    color = "grey50"
  ) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  scale_x_date(breaks = "16 months", date_labels = "%b %Y") +
  scale_fill_manual(
    values = cols_chains, labels = names(cols_chains),
    name = "Chain", drop = TRUE, aesthetics = c("color", "fill")
  ) +
  labs(x = "", y = bquote("Moving average of"~R[e]~"(6 mos)")) +
  theme_minimal_hgrid(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_ts

# Individual level Re vs. cov ----
# join with coverage & apply negbinom regression to get coefficient
cov_long <- 
  coverage %>%
  pivot_longer(-V1, names_to = "month", values_to = "coverage") %>%
  mutate(month = as.numeric(gsub("V","", month)))

cov_long <- cov_long %>% 
  mutate(month = as_date(floor_date(ymd("2010-12-01") + dmonths(month), 
                                    unit = "months"))) %>%
  group_by(month, V1) %>%
  slice_max(coverage, n = 1, with_ties = FALSE) %>%
  ungroup()

cov_join <- expand.grid(V1 = unique(cov_long$V1), 
                        month = seq.Date(ymd("2010-01-01"), ymd("2019-01-01"), 
                                         by = "month"))

cov_join %>%
  left_join(cov_long) %>%
  tidyr::replace_na(list(coverage = 0)) -> cov_long

links_re %>%
  mutate(month = floor_date(symptoms_started, unit = "months"),
         year = year(symptoms_started),
         V1 = paste(district, ward, sep = "_")) %>%
  left_join(cov_long) -> reff_cov

plot_cov <-
  ggplot(reff_cov) +
  geom_point(aes(x = coverage, y = reff, fill = factor(lineage_chain)),
             size = 2.5, alpha = 0.5, shape = 21, color = "NA"
  ) +
  geom_point(aes(x = coverage, y = reff, color = factor(lineage_chain)),
             size = 2.5, alpha = 1, shape = 21, fill = "NA"
  ) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  scale_fill_manual(
    values = cols_chains, labels = names(cols_chains),
    name = "Chain", drop = TRUE, aesthetics = c("color", "fill")
  ) +
  labs(x = "Vaccination coverage", y = bquote(R[e])) +
  theme_minimal_hgrid(font_size = 12)
plot_cov
ggsave("figures/supplement/re_cov.jpeg", plot_cov, height = 5, width = 5)


# how about the mcc tree
ttrees[scenario == 1 & mcc == TRUE] %>%
  left_join(case_date_locs) %>%
  mutate(month = floor_date(symptoms_started, unit = "months"),
         year = year(symptoms_started),
         V1 = paste(district, ward, sep = "_")) %>%
  left_join(cov_long) -> reff_mcc

plot(reff_mcc$coverage*100, jitter(reff_mcc$reff))

mod_nb <- MASS::glm.nb(reff ~ coverage, data = reff_mcc) # not at all significant
plot(reff_mcc$coverage, reff_mcc$reff)
summary(mod_nb)

mod <- glmm(reff ~ coverage, data = reff_mcc, class="Poisson") # not at all significant

re_plot_main <-
  (plot_ts / plot_cov) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(legend.position = "bottom") 

ggsave("figures/main/re_plot.jpeg", re_plot_main, height = 8, width = 6)

# Add the exposures by case coloured by lineage
lin_exps <- links_re %>% 
  filter(species == "Domestic dog")  %>% 
  group_by(exposed, lineage_chain) %>% 
  tally()

plot_exp <- 
  ggplot(lin_exps, aes(fill=factor(lineage_chain), y=n, x=exposed)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Persons exposed \n per rabid dog", y = "Frequency") + 
  scale_fill_manual(values = alpha(cols_chains[1:7], 0.6), labels = names(cols_chains),
    name = "Chain", drop = TRUE, aesthetics = c("color", "fill")) + 
  theme_minimal_hgrid(font_size = 12) +
  theme(legend.position="none", axis.text.x = element_text(hjust = 0))
plot_exp


# Combine
plot_lineages <-
  (plot_ts + plot_exp) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") # & theme(legend.position = "bottom") 
plot_lineages 

ggsave("figures/Figure_4_chains.jpeg", plot_lineages , height = 5, width = 8)

