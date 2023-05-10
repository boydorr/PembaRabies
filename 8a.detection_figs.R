# Estimated detection ----

# pkgs & scripts ----
library(treerabid)
library(data.table)
library(tidyverse)
library(lubridate)
library(magrittr)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggridges)
library(here)
library(ggbeeswarm)
select <- dplyr::select

# data ----
recover_detect <- fread("Output/trees/detection_sims.csv")
est_detect <- fread("Output/trees/detection_ests.csv")
scenarios <- fread("Output/trees/scenarios.csv")

# join up with scenarios ----
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
est_detect <- est_detect[scenarios, on = "scenario"]


recov_long <- melt(recover_detect, id.vars = c("true", "sim", "nobs"))
recov_summ <- 
  recov_long[, .(mean = mean(value), max = max(value), min = min(value)), 
             by = c("true", "sim", "nobs", "variable")]
recov_summ <- filter(recov_summ, nobs != 1000, variable != "est_unsorted")
recov <-
  ggplot(recov_summ) +
  geom_pointrange(aes(x = true, y = mean, ymin = min, ymax = max, 
                      fill = factor(nobs), color = factor(nobs)), 
                  shape = 21, 
                  alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_abline(slope = 1, intercept = -0.1, linetype = 2, color = "grey50") +
  scale_fill_brewer(palette = "Accent",
                    name = "Number of observations") +  
  scale_color_manual(values = c("black", NA), name = "Number of observations") +
  labs(x = "True detection prob", y = "Estimated detection prob") +
  theme_minimal_hgrid(font_size = 12) +
  guides(fill = guide_legend(override.aes = list(shape = 21, alpha = 1)))

# estimated detects ----
est_summ <-
  est_detect[, .(mean_detect = mean(detection)), by = c("era", "sim", "data_used", "prune_type")]
time_labs <- c("2010 - 2014", "2016 - 2019", "Overall")
names(time_labs) <- c("First", "Second", "Overall")

est_qs <- # But the 95%CIs are by sim - NOT across sims!
  est_detect[, .(mean_detect = mean(detection),
                 LCI = quantile(detection, 0.025),
                 UCI = quantile(detection, 0.975)), 
             by = c("era", "sim", "data_used", "prune_type")]
est_qs[prune_type == "Unpruned" & data_used == "Epi & genetic data" & era == "First"]$sim

# estimate total number of rabid dogs in each outbreak
first_detect_U <- est_summ[prune_type == "Unpruned" & 
                         data_used == "Epi & genetic data" & 
                         era == "First"]$mean_detect
quantile(first_detect_U, c(0.025, .5, .975))
mean(first_detect_U)

second_detect_U <- est_summ[prune_type == "Unpruned" & 
                             data_used == "Epi & genetic data" & 
                             era == "Second"]$mean_detect
quantile(second_detect_U, c(0.025, .5, .975))
mean(second_detect_U)

first_detect_TD <- est_summ[prune_type == "Pruned by time & distance" & 
                             data_used == "Epi & genetic data" & 
                             era == "First"]$mean_detect
quantile(first_detect_TD, c(0.025, .5, .975))
mean(first_detect_TD)
92/mean(first_detect_TD) # estimated cases in first outbreak

second_detect_TD <- est_summ[prune_type == "Pruned by time & distance" & 
                              data_used == "Epi & genetic data" & 
                              era == "Second"]$mean_detect
quantile(second_detect_TD, c(0.025, .5, .975))
mean(second_detect_TD)
97/mean(second_detect_TD)


ests_all <- 
  ggplot(data = est_summ) +
  geom_density(aes(x = mean_detect, 
                   fill = factor(era), group = factor(era)), 
               alpha = 0.75) +
  scale_fill_brewer(palette = "Dark2", name = "Time period", 
                    labels = time_labs) +
  theme_minimal_hgrid(font_size = 12) +
  labs(x = "Estimated detection probability", y = "Density") +
  facet_grid(data_used ~ prune_type, labeller = label_wrap_gen(20))
  
ggsave("figures/supplement/sfig_detection_se.jpeg", ests_all, height = 8, 
       width = 10, bg = "white")

ests_main <-
  ggplot(data = est_summ[prune_type == "Unpruned" & 
                         data_used == "Epi & genetic data"]) +
  geom_density(aes(x = mean_detect, 
                   fill = factor(era), group = factor(era)), 
               alpha = 0.75) +
  scale_fill_brewer(palette = "Dark2", name = "Time period",
                    labels = time_labs) +
  theme_minimal_hgrid(font_size = 12) +
  labs(x = "Estimated detection probability", y = "Density")

ests_first <- ggplot(data = est_summ[prune_type == "Unpruned" & 
                                        data_used == "Epi & genetic data" &
                                         era == "First"]) +
  geom_density(aes(x = mean_detect), alpha = 0.75) +
  theme_minimal_hgrid(font_size = 12) +
  labs(x = "Detection probability", y = "Density") + scale_x_continuous(lim = c(0,1))
ests_second <- ggplot(data = est_summ[prune_type == "Unpruned" & 
                                       data_used == "Epi & genetic data" &
                                       era == "Second"]) +
  geom_density(aes(x = mean_detect), alpha = 0.75) +
  theme_minimal_hgrid(font_size = 12) +
  labs(x = "Detection probability", y = "Density") + scale_x_continuous(lim = c(0,1))
ests_by_time <- (ests_first  / ests_second) 
ggsave("figures/ests_by_time.pdf", ests_by_time, height = 6, width = 4)


# probability of detecting at least one case given a chain size of x
chain_dets <- expand_grid(est_summ[prune_type == "Unpruned" & 
                            data_used == "Epi & genetic data"], 
                          size = seq_len(20))
chain_dets <- mutate(chain_dets, prob =1 - (1 - mean_detect)^size)

chain_dets_summary <- 
  chain_dets %>%
  group_by(size, era) %>%
  summarize(mean_prob = mean(prob))

chains <-
  ggplot() +
  geom_segment(data = chain_dets_summary, 
               aes(x = size - 1, xend = size + 1, 
                   y = mean_prob, yend = mean_prob, color = factor(era)), 
               size = 1.25) + 
  geom_quasirandom(data = chain_dets,
                   aes(x = size, y = prob, color = factor(era), 
                       group = era), alpha = 0.25) +
  labs(aes(x = "Probability of detecting chain of size x")) +
  scale_color_brewer(palette = "Dark2", name = "Time period", 
                     labels = time_labs, guide = "none") +
  theme_minimal_hgrid(font_size = 12) +
  labs(y = "Probability of detecting \nat least one case", x = "Chain size")

sfig_detection <- 
  (recov  / ests_main / chains) + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") 
ggsave("figures/supplement/sfig_detection_main.jpeg", sfig_detection, height = 10, width = 8)

