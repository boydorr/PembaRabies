# Estimating detection probabilities from trees ----

# Takes abt 20 minutes (particularly the sim kappa part)

# sub_cmd:=-t 2 -n 5 -jn test -wt 1m -sn

if(Sys.getenv("SLURM_JOB_ID") != "") {
  ncores <- as.numeric(Sys.getenv("SLURM_NTASKS")) 
} else {
  ncores <- parallel::detectCores() - 1
}

print(ncores)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# Packages 
library(treerabid)  
library(data.table)
library(foreach)
library(iterators)
library(doRNG)
library(tidyr)
library(dplyr)
library(lubridate)

# Estimate detection ----
ttrees <- fread("Output/trees/trees_all.gz")
scenarios <- fread("Output/trees/scenarios.csv")

# choose 100 random, mcc, & majority for each ----
tree_ids <- ttrees[, .(maj = sum(majority), mcc = sum(mcc)), 
                   by = c("scenario", "sim")][maj > 0 | mcc > 0]
mcc_maj_trees <- ttrees[tree_ids[, .(scenario, sim)], on = c("sim", "scenario")]

# sample trees
set.seed(1677)
selected <- ttrees[!tree_ids[, .(scenario, sim)], on = c("sim", "scenario")][, .(sim = sample(sim, 100)), by = "scenario"]
comp_times_dist <- rbind(mcc_maj_trees, ttrees[selected, on = c("sim", "scenario")])
fwrite(comp_times_dist, "Output/trees/comp_diffs.gz")

# known and select cols
known <- unique(ttrees[type == "traced"]$id_case)
comp_times_dist[, known_kappa := ifelse(id_case %in% known, 1, 0)]
comp_times_dist <- comp_times_dist[, c("t", "t_diff", "scenario", "sim", "known_kappa")]

# Split up and run estimates (overall & by time period)
comp_times_dist[, levs := interaction(sim, 
                                      scenario,
                                      drop = TRUE)]
comp_times_dist <- comp_times_dist[!is.na(t_diff)]
t_diffs <- split(comp_times_dist, comp_times_dist$levs)

system.time({
  est_detect <-
    foreach(k = iter(t_diffs), .combine = rbind) %do% {
      
      split_by_time <- list(overall = k, 
                            first = k[t < ymd("2015-01-01")], 
                            second = k[t > ymd("2015-01-01")])
      ests <- lapply(split_by_time, 
                      function(x) {
                        out <- fit_sims_pi(t_diff = x$t_diff, 
                                 nsims = 10, 
                                 candidate_pis = seq(0.01, 0.99, by = 0.01),
                                 si_fun = treerabid::si_fun_lnorm, 
                                 params = treerabid::params_treerabid, 
                                 alpha = 0.01, 
                                 known_kappas = x$known_kappa,
                                 seed = as.numeric(x$levs[1]))
                        
                         data.table(detection = out, 
                                    x[1, -c("known_kappa", "t_diff", "levs")], 
                                    nobs = nrow(x))
                      }
                    )
      ests[[1]]$era <- "Overall"
      ests[[2]]$era <- "First"
      ests[[3]]$era <- "Second"
      
      rbindlist(ests)
      
    }
})

# Also do a check against simulated kappa vals across nobs ----
set.seed(1456)
sim_kappa <- 
  rbindlist(
    lapply(
      seq_len(1000),  
      function(x) {
        z <- runif(1, min = 0.01, max = 0.99)
        rbindlist(
        lapply(c(100, 200), function(n) {
          t_diff <- sim_times_pi(si_fun_lnorm, nobs = n, 
                                 params = treerabid::params_treerabid, alpha = 0.01,
                                 pi = z)
          ests_sorted <- fit_sims_pi(t_diff, nsims = 10, 
                                     candidate_pis = seq(0.01, 0.99, by = 0.01),
                                     si_fun_lnorm, params = treerabid::params_treerabid, 
                                     alpha = 0.01, seed = round(z) * 1000)
          data.table(true = z, est_sorted = ests_sorted, 
                     sim = x, nobs = n)
        })
      )
    }
  )
)

# Each tree lookup ----
parallel::stopCluster(cl)

# Write out files of the trees and the links (consensus & all)
fwrite(est_detect, "Output/trees/detection_ests.csv")
fwrite(sim_kappa, "Output/trees/detection_sims.csv")

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/Serengeti_Rabies/output/"
syncfrom <- "mrajeev@della.princeton.edu:Serengeti_Rabies/output/trees"
