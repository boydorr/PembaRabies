# Run trees with treerabid ----
# Takes about 13 minutes with three cores

library(here)
library(treerabid)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(sf)
library(lubridate)
library(raster)
library(foreach)
library(iterators)
library(doRNG)
library(doParallel)
library(igraph)
library(cowplot)
library(ggplot2)
library(ggraph)
library(patchwork)
select <- dplyr::select

# data
case_dt <- read_csv("Output/PembaAnimalCases.csv")
pemba_shp <- st_read(here("Output/GIS/PembaVill_NBS2012/PembaVill_NBS2012_cleaned.shp"))

case_dt %>%
  janitor::clean_names() %>%
  mutate(symptoms_started = ymd(symptoms_started)) %>%
  dplyr::filter(!is.na(symptoms_started),
                !is.na(utm_easting),
                !is.na(utm_northing),
                symptoms_started >= ymd("2010-01-01"),
                suspect %in% "Yes") %>%
  # get uncertainty in days
  mutate(days_uncertain = case_when(symptoms_started_accuracy == "+/- 14 days" ~ 14L,
                                    symptoms_started_accuracy == "+/- 7 days" ~ 7L,
                                    symptoms_started_accuracy == "+/- 28 days" ~ 28L,
                                    symptoms_started_accuracy == "0" ~ 0L,
                                    TRUE ~ 0L),
         owned = ifelse(owner %in% "Known", TRUE, FALSE),
         exclude_progen = ifelse(species %in% "Domestic dog", FALSE, TRUE)) -> case_dt

lineages <- read_csv("data/pembaSeq_geneticClusters.csv")
lin_meta <- read_csv("data/pembaSeq_epiLab.csv")

lineages %>%
  filter(!is.na(Sample_sequenceID)) %>%
  mutate(cluster2 = as.numeric(factor(coalesce(alt, cluster)))) %>%
  dplyr::select(Sample.ID, Sample_sequenceID, cluster, cluster2, Notes) %>%
  left_join(dplyr::select(lin_meta, Sample.ID, ID)) %>%
  janitor::clean_names() %>%
  right_join(case_dt, by = "id") %>%
  tidyr::replace_na(list(cluster = 0, cluster2 = 0)) -> case_dt
case_dt$nocluster <- 0 # null lineages

# duplicate check: max(tabulate(case_dt$id)) should = 1
fwrite(case_dt, "Output/trees/case_dt_cleaned.csv")

# ct data with dates
case_dates <- data.table(id_case = case_dt$id,
                         symptoms_started = case_dt$symptoms_started)

lineages <- data.table(select(case_dt, id_case = id, lineage = cluster2, 
                              date_sampled = symptoms_started))
fwrite(lineages, "Output/trees/lineage_dt_cleaned.csv")

# Run trees using epi & gen data with & w/out pruning ----
# weibull for distance & lognormal for serial
tree_pars <- tidyr::expand_grid(si_pdist = "lnorm", 
                                dist_pdist = "weibull",
                                convolve = "baseline",
                                prune = TRUE, 
                                time_cutoff = c(1, 0.99),
                                dist_cutoff = c(1, 0.99),
                                use_known = TRUE, 
                                use_gen = c(TRUE, FALSE),
                                nsim = 1000)
tree_pars %<>%
  # don't actually want to cutoff only by distance
  filter(!(dist_cutoff == 0.99 & time_cutoff == 1)) %>%
  mutate(prune = ifelse(time_cutoff == 1 & dist_cutoff == 1, FALSE, TRUE),
         scenario = row_number())
tree_pars$seed <- 42 * 1:nrow(tree_pars) + 42


# Run in parallel -----
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY = FALSE, fill = TRUE)
}

system.time({
  trees_all <-
    foreach(i = iter(tree_pars, by = "row"), .combine = comb) %do% {
      
      
      if(i$use_gen) {
        lin_dt <- copy(lineages)
      } else {
        lin_dt <- data.table(id_case = case_dt$id, lineage = 0)      
      }
      
      ttrees <-
        boot_trees(id_case = case_dt$id,
                   id_biter = case_dt$biter_id,
                   x_coord = case_dt$utm_easting,
                   y_coord = case_dt$utm_northing,
                   owned = case_dt$owned,
                   date_symptoms = case_dt$symptoms_started, # needs to be in a date class
                   days_uncertain = case_dt$days_uncertain,
                   use_known_source = i$use_known,
                   exclude_progen = case_dt$exclude_progen,
                   lineages = lin_dt,
                   prune = i$prune,
                   si_fun = treerabid::si_lnorm1,
                   dist_fun = treerabid::dist_weibull1,
                   cutoff = c(time = i$time_cutoff, dist = i$dist_cutoff),
                   params = treerabid::params_treerabid,
                   N = 1000,
                   seed = i$seed)
      
      # way to join back up with scenarios
      id <- i[, "scenario"]
      
      # Summarize the trees
      known_progens <- unique(ttrees[type == "traced"]$id_case)
      links_all <- build_all_links(ttrees, N = i$nsim)
      links_consensus_consistent <- 
        build_consensus_links(links_all, case_dates = case_dates,
                              lineages = lin_dt, known_progens = known_progens, 
                              link_all = TRUE)
      
      # consensus without loops fixed
      links_consensus_raw <- links_all[links_all[, .I[which.max(links)], 
                                                 by = c("id_case")]$V1]
      
      # majority & mcc
      tree_ids <- c(mcc = 
                      build_consensus_tree(links_consensus_consistent, ttrees, links_all,
                                           type = "mcc", output = "sim"), 
                    majority = 
                      build_consensus_tree(links_consensus_consistent, ttrees, links_all,
                                           type = "majority", output = "sim"))
      ttrees$mcc <- ifelse(ttrees$sim %in% tree_ids["mcc"], 1, 0)
      ttrees$majority <- ifelse(ttrees$sim %in% tree_ids["majority"], 1, 0)
      
      # also get reff
      reff <- ttrees[, .(reff = .N), 
                     by = c("sim", "id_progen")][, 
                                          .(sim, id_case = id_progen, reff)]
      ttrees <- reff[ttrees, on = c("sim", "id_case")]      
      setnafill(ttrees, cols = "reff", fill = 0)
      
      list(ttrees = cbind(ttrees, id), 
           links_consensus_consistent = cbind(links_consensus_consistent, id), 
           links_consensus_raw = cbind(links_consensus_raw, id))
      
    }
})

parallel::stopCluster(cl)

# Output files ---
fwrite(trees_all$ttrees, "Output/trees/trees_all.gz")
fwrite(trees_all$links_consensus_consistent, "Output/trees/links_consensus_consistent.csv") 
fwrite(trees_all$links_consensus_raw, "Output/trees/links_consensus_raw.csv") # may not need this in the end

# for joining back up with metadata on parameters used (don't want to add cols to 
# the full large data.tables)
fwrite(tree_pars, "Output/trees/scenarios.csv")


