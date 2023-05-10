# tree sensitivity ----

# pkgs ----
library(data.table)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(lubridate)
library(glue)
library(igraph)
library(ggraph)
library(treerabid)
source("R/plot_lineage_ts.R")

# load in outputs ----
ttrees <- fread("Output/trees/trees_all.gz")
links_consensus <- fread("Output/trees/links_consensus_raw.csv") 
tree_consensus <- fread("Output/trees/links_consensus_consistent.csv")
scenarios <- fread("Output/trees/scenarios.csv")
case_dt <- fread("Output/trees/case_dt_cleaned.csv")
case_dates <- case_dt[, .(id_case = id, symptoms_started)]
lineages <- fread("Output/trees/lineage_dt_cleaned.csv")
comp_times_dist <- fread("Output/trees/comp_diffs.gz")

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

# topological uncertainty w/in trees ----
# join with case dates
ttrees <- ttrees[case_dates, on = "id_case"]

# also get scores in terms of most frequently chosen
links_all <- ttrees[, .(links = .N, 
                        prob = .N/202), by = c("scenario", "id_case", "id_progen")]
ttrees <- links_all[, c("id_case", "id_progen", "prob", "scenario")][ttrees, on = c("id_case", "id_progen", "scenario")]
links_consensus[, chosen := TRUE]
ttrees <- links_consensus[, c("id_case", "id_progen", "chosen", "scenario")][ttrees, on = c("id_case", "id_progen", "scenario")]
ttrees[, chosen := ifelse(is.na(chosen), FALSE, TRUE)]

# score the trees filtering out the known progens
known_progens <- unique(ttrees[type == "traced"]$id_case)
topo <- ttrees[!(id_case %in% known_progens)][,
                           .(score_prob = prod(prob), 
                           score_freq = sum(chosen)/.N), by = c("scenario", "sim")]
topo <- scenarios[topo, on = "scenario"]
mcc_trees <- ttrees[, .(mcc = sum(mcc)), by = c("scenario", "sim")][mcc > 0]
mcc_trees <- topo[mcc_trees, on = c("scenario", "sim")]

sfig_topo_uncertainty <-
  ggplot(data = topo) + 
  ggridges::geom_density_ridges(aes(y = as.numeric(prune_type), x = score_freq, fill = data_used,
                                    color = data_used, 
                                    group = interaction(as.numeric(topo$prune_type), topo$data_used)), 
                                alpha = 0.5) +
  geom_segment(data = mcc_trees, aes(x = score_freq, xend = score_freq, 
                                     y = as.numeric(prune_type) - 0.15,
                                     yend = as.numeric(prune_type) + 0.15, 
                                     color = data_used), 
               linetype = 2) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = levels(topo$prune_type)) +
  labs(x = "Proportion of consensus links in tree \n (N = 1000 trees)", y = "") +
  scale_fill_manual(values = data_used_col, aesthetics = c("color", "fill"), 
                    name = "Data used") +
  theme_minimal_hgrid(font_size = 12)

# agreement between trees ----
links_consensus <- scenarios[links_consensus, on = "scenario"]
cons_comp <- links_consensus[links_consensus, on = "id_case", allow.cartesian = TRUE]
setnafill(cons_comp, cols = c("id_progen", "i.id_progen"), fill = 0)
cons_comp[, match := id_progen == i.id_progen]
cons_comp <-
  cons_comp[, .(prop_matching = sum(match)/.N), 
            by = c("data_used", "prune_type", 
                   "i.data_used", "i.prune_type")]
cons_comp[, c("same_data", "same_prune") := 
            .(data_used == i.data_used, prune_type == i.prune_type)]
fwrite(cons_comp, "Output/trees/tree_agreement.csv")

# progenitor vs. lineage probabilities ----
links_consensus[, type := "Highest probability progenitor"]
links_prob <- links_consensus[use_gen == TRUE]
trees_gen <- scenarios[ttrees, on = "scenario"][use_gen == TRUE & lineage == 0]

# also filter out known progens
gen_summary <- trees_gen[, .(prob = .N/1000, type = "Highest probability lineage"), 
                         by = c("lineage_chain", "scenario", "id_case")]
gen_summary %>% 
  group_by(scenario, id_case) %>% 
  filter(prob == max(prob)) %>%
  select(id_case, prob, type, scenario) %>%
  bind_rows(select(links_prob, id_case, prob, type, scenario)) %>%
  left_join(scenarios) %>%
  filter(!id_case %in% known_progens) %>%
  as.data.table() -> prob_comp

prob_comp_cols <- c("#3D3131", "#F44242")
sfig_prob_comps <-
  ggplot(data = prob_comp) + 
  ggridges::geom_density_ridges2(aes(y = as.numeric(prune_type), x = prob, fill = type,
                                    color = type,
                                    group = interaction(as.numeric(prune_type), type)), 
                                alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = levels(prob_comp$prune_type)) +
  labs(x = "Probability of link or lineage \n selected most frequently", y = "") +
  scale_fill_manual(values = prob_comp_cols, aesthetics = c("color", "fill"), 
                    name = "Type of link") +
  theme_minimal_hgrid(font_size = 12) +
  theme(axis.text.y = element_blank())

# lineage heat map ----
lin_dates <- lineages[, .(earliest_date = min(date_sampled), 
                          lineage_chain = lineage), by = "lineage"] # and zero should come first
all_cases <- expand.grid(id_case = unique(gen_summary$id_case), 
                         lineage_chain = lin_dates$lineage_chain, 
                         scenario = scenarios$scenario)
gen_summary <- gen_summary[all_cases, on = c("id_case", "lineage_chain", "scenario")]
gen_summary <- scenarios[gen_summary, on = "scenario"][use_gen == TRUE]
gen_summary <- case_dates[gen_summary, on = "id_case"]
gen_summary <- lin_dates[gen_summary, on = "lineage_chain"]

setnafill(gen_summary, cols = "prob", fill = 0)
lin_dates %>%
  mutate(labs = ifelse(lineage_chain == 0, "Unsampled",
                       glue("L{lineage_chain} ({earliest_date})")), 
         name = interaction(earliest_date, lineage)) -> lin_dates
lin_labs <- lin_dates$labs
names(lin_labs) <- lin_dates$name
gen_summary$prune_type <- factor(gen_summary$prune_type, 
                                 levels = c("Pruned by time & distance",
                                            "Pruned by time",  
                                            "Unpruned"))

period_sepdate <- as.character(max(gen_summary$symptoms_started[gen_summary$symptoms_started < "2015-01-01"]))
sfig_lin_probs <-
  ggplot(gen_summary) +
  geom_tile(aes(y = reorder(interaction(earliest_date, lineage), desc(earliest_date)), 
                x = reorder(symptoms_started, symptoms_started), 
                fill = prob)) +
  scale_y_discrete(labels = lin_labs) +
  geom_vline(xintercept = period_sepdate, ) +
  scale_fill_distiller(direction = 1) +
  facet_wrap(~prune_type, ncol = 1, strip.position = "right", 
             labeller = label_wrap_gen()) + 
  labs(y = "Lineages ordered by date", x = "Unsampled cases ordered by date", 
       fill = "Probability") + 
  theme_minimal_hgrid(font_size = 12) +
  theme(axis.text.x = element_blank())

sfig_tree_se <- 
  ((sfig_topo_uncertainty | sfig_prob_comps) + plot_layout(guides = "collect")) / sfig_lin_probs + 
  plot_annotation(tag_levels = "A") 

ggsave("figures/supplement/sfig_tree_se.jpeg", height = 10, width = 10)

# time & distance distributions between cases (split by era) ----
comp_diffs <- scenarios[comp_times_dist, on = "scenario"]
comp_diffs[, era := ifelse(t < "2015-01-01", "2010 - 2014", 
                           "2016 - 2020")]
comp_diffs <- comp_diffs[type != "traced"]
era_cols <- alpha(scales::brewer_pal(palette = "Dark2")(2), 0.1)

# Add in the underlying distributions
t_max <- max(comp_diffs$t_diff, na.rm = TRUE) + 30
dist_max <- max(comp_diffs$dist_diff, na.rm = TRUE) + 1000

# times and distances
ref_d <-  data.frame(dist_diff = rweibull(10000, 
                                          shape = params_treerabid$DK_shape_weibull, 
                                          scale = params_treerabid$DK_scale_weibull))
ref_t <- data.frame(t_diff = rlnorm(10000, 
                                    meanlog = params_treerabid$SI_meanlog, 
                                    sdlog = params_treerabid$SI_sdlog))

sfig_tdiffs <-
  ggplot(comp_diffs) +
  geom_density(aes(x = t_diff/7, group = interaction(sim, era, 
                                                   scenario),
                   color = era)) +
  scale_x_continuous(trans = "log", breaks = scales::breaks_log()) +
  labs(x = "Time between linked cases (weeks)", y = "Density") + 
  scale_color_manual(values = era_cols, name = "Time period") +
  geom_density(data = ref_t, aes(x = t_diff/7), linetype = 2) +
  facet_grid(data_used ~ prune_type) +
  theme_half_open(font_size = 12) +
  guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                  fill = era_cols)))

ggsave("figures/supplement/sfig_tdiffs.jpeg", sfig_tdiffs, height = 8, width = 8, 
       bg = "white")

sfig_ddiffs <- 
  ggplot(comp_diffs) +
  geom_density(aes(x = scales::oob_squish(dist_diff, c(100, Inf)), 
                   group = interaction(sim, era,  scenario),
                   color = era)) +
  scale_x_continuous(trans = "log", breaks = scales::breaks_log()) +
  labs(x = "Distance between linked cases (meters)", y = "Density") + 
  scale_color_manual(values = era_cols, name = "Time period") +
  geom_density(data = ref_d, 
               aes(x = scales::oob_censor(dist_diff, c(100, Inf))), 
               linetype = 2) +
  facet_grid(data_used ~ prune_type) +
  theme_half_open(font_size = 12) +
  guides(color = guide_legend(override.aes = list(alpha = 1, 
                                                  fill = era_cols)))

ggsave("figures/supplement/sfig_ddiffs.jpeg", sfig_ddiffs, height = 8, width = 8, 
       bg = "white")

# maj | mcc | consensus tree comparison ----
tree_consensus <- tree_consensus[scenarios, on = "scenario"]
con_plots <- 
  lapply(split(tree_consensus,
               interaction(tree_consensus$prune_type, 
                           tree_consensus$data_used)), 
       function(x) {
         if(x$data_used[1] == "Epi data only") {
           x <- x[, -c("lineage_chain", "lineage")]
           x <- lineages[x, on = "id_case"]
           x$lineage_chain <- x$lineage
         }
         out <- plot_lin_ts(x, case_dates, unsampled = 21)
         out <- out + labs(subtitle = paste(x$data_used[1], "\n", 
                                            x$prune_type))
       }
)

consensus_plots_ts <- 
  wrap_plots(con_plots) + 
  plot_annotation(title = "Consensus tree comparison") &
  panel_border(size = 0.1, color = "black")

ggsave("figures/supplement/consensus_tree_check.jpeg", 
       consensus_plots_ts, height = 12, width = 12)

# mcc trees
tree_mcc <- ttrees[mcc_trees, on = c("scenario", "sim")][scenarios, on = "scenario"]

mcc_plots <- 
  lapply(split(tree_mcc,
               interaction(tree_mcc$prune_type, 
                           tree_mcc$data_used)), 
         function(x) {
           if(x$data_used[1] == "Epi data only") {
             x <- x[, -c("lineage_chain", "lineage")]
             x <- lineages[x, on = "id_case"]
             x$lineage_chain <- x$lineage
           }
           out <- plot_lin_ts(x, case_dates, unsampled = 21)
           out <- out + labs(subtitle = paste(x$data_used[1], "\n", 
                                              x$prune_type))
         }
  )


mcc_plots_ts <- 
  wrap_plots(mcc_plots) + 
  plot_annotation(title = "MCC tree comparison") &
  panel_border(size = 0.1, color = "black")

ggsave("figures/supplement/mcc_tree_check.jpeg", 
       mcc_plots_ts, height = 12, width = 12)

# maj trees
maj_trees <- ttrees[, .(maj = sum(majority)), by = c("scenario", "sim")][maj > 0]
tree_maj <- ttrees[maj_trees, on = c("scenario", "sim")][scenarios, on = "scenario"]

maj_plots <- 
  lapply(split(tree_maj,
               interaction(tree_maj$prune_type, tree_maj$data_used)), 
         function(x) {
           if(x$data_used[1] == "Epi data only") {
             x <- x[, -c("lineage_chain", "lineage")]
             x <- lineages[x, on = "id_case"]
             x$lineage_chain <- x$lineage
           }
           out <- plot_lin_ts(x, case_dates, unsampled = 21)
           out <- out + labs(subtitle = paste(x$data_used[1], "\n", 
                                              x$prune_type))
         }
  )


maj_plots_ts <- 
  wrap_plots(maj_plots) + 
  plot_annotation(title = "Majority tree comparison") &
  panel_border(size = 0.1, color = "black")

ggsave("figures/supplement/majority_tree_check.jpeg", 
       maj_plots_ts, height = 12, width = 12)

