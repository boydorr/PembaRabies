# Make figures for supplement ----

# pkgs ----
library(data.table)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(lubridate)
library(sf)
library(ggforce)
library(animation)
library(ggraph)
library(igraph)
library(treerabid)
source("R/plot_lineage_ts.R")

# load in outputs ----
ttrees <- fread("Output/trees/trees_all.gz")
links_consensus <- fread("Output/trees/links_consensus_consistent.csv") 
scenarios <- fread("Output/trees/scenarios.csv")
case_dt <- fread("Output/trees/case_dt_cleaned.csv")
case_dates <- case_dt[, .(id_case = id, symptoms_started)]
pemba_shp <- st_read("Output/GIS/PembaVill_NBS2012/PembaVill_NBS2012_cleaned.shp")

# filter to the main one to present ---
# unpruned using both epi & genetic data
links_consensus <- links_consensus[scenario == 1]
ttrees <- ttrees[scenario == 1]
mcc_tree <- ttrees[mcc == 1]
maj_tree <- ttrees[majority == 1]

out <- plot_lin_ts_vertical(links_consensus, case_dates)
out[[1]] <- 
  out[[1]] + 
  labs(tag = "B") + 
  theme(axis.line.y = element_line(color = "grey", size = 0.75), 
        plot.tag = element_text(face = "bold"))

out_h <- plot_lin_ts(links_consensus, case_dates)
out_h[[1]] <- 
  out_h[[1]] + 
  labs(tag = "B") + 
  theme(axis.line.y = element_line(color = "grey", size = 0.75), 
        plot.tag = element_text(face = "bold"))


# map & time series
case_dt <- data.table(case_dt)
links_consensus <- links_consensus[case_dt[, .(id_case = id, utm_easting, utm_northing, symptoms_started)], on = "id_case"]
links_consensus[, period := ifelse(year(symptoms_started) < 2015, "2010 - 2016", "2016 - 2020")]
links_consensus$lineage_chain <- factor(links_consensus$lineage_chain , levels = 1:7)
shps = c(NA, rep(22, 9))
names(shps) <- c(0, 1:9)
cols_chains <- c("grey50",
                 "#ebac23", "#b80058", "#006e00", "#00bbad",
                 "#ff9287", "#5954d6", "#d163e6", "#b24502", "#00c6f8",
                 "#878500", "#008cf9")
names(cols_chains) <- c(0:9)
cols_chains <- cols_chains[2:(max(as.numeric(links_consensus$lineage_chain)) + 1)]
shps <- shps[2:(max(as.numeric(links_consensus$lineage_chain)) + 1)]

pemba_maps <-
  ggplot() +
  geom_sf(data = pemba_shp, fill = "light grey", color = "light grey") + # Change colour here!
  geom_point(data = links_consensus,
             aes(x = utm_easting, y = utm_northing,
                 color = factor(lineage_chain),
                 shape = ifelse(lineage == 0, 16, 15),
                 size = ifelse(lineage == 0, 1.5, 2.5), 
                 alpha = ifelse(lineage == 0, 0.8, 1)),
             stroke = 1.5) +
  geom_point(data = filter(links_consensus, lineage != 0),
             aes(x = utm_easting, y = utm_northing),
             color = "black", size = 2.5, shape = 0,
             stroke = 1) +
  scale_color_manual(values = cols_chains, labels = names(cols_chains),
                     name = "Lineage", drop = TRUE) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  facet_wrap(~ period, ncol = 2) +
  theme_map() +
  theme(legend.position = "bottom", plot.margin = margin(0, 0, 0, 0)) +
  guides(color = guide_legend(override.aes = list(shape = 15))) +
  labs(tag = "C")

links_consensus %>%
  mutate(month = floor_date(symptoms_started, unit = "months")) %>%
  group_by(month, lineage_chain) %>%
  summarise(cases = n(),
            start_date = min(month)) -> top_chain_hist

# chain_ts <-
#   ggplot(top_chain_hist) +
#   geom_col(aes(x = as.numeric(month),
#                y = cases, fill = lineage_chain),
#            position = position_stack(reverse = FALSE)) +
#   scale_y_reverse() +
#   scale_x_reverse() +
#   scale_fill_manual(values = cols_chains, guide = "none") +
#   labs(x = "", y = "Number of cases") +
#   cowplot::theme_minimal_vgrid() +
#   coord_flip() +
#   theme(axis.text.y = element_blank()) +
#   labs(tag = "A")

chain_ts <-
  ggplot(top_chain_hist) +
  geom_col(aes(x = month, y = cases, fill = lineage_chain)) +
  scale_x_date(date_breaks = "15 months", date_labels = "%b %Y") + 
  scale_fill_manual(values = cols_chains, guide = "none") +
  labs(x = "", y = "Number of cases") +
  cowplot::theme_minimal_hgrid() +
  labs(tag = "A")

chains_main <- 
  (chain_ts + 
   out_h +
   pemba_maps) + 
  plot_layout(ncol = 3, widths = c(0.75, 1.75, 1.5), guides = "collect") &
  theme(legend.position = "bottom")

chains_main <- 
  (chain_ts + 
     out_h +
     pemba_maps) + 
  plot_layout(ncol = 1, heights = c(1.75, 1, 2.5), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/main/chains_map_ts.jpeg", chains_main, height = 8, width = 8)


# animations ----
links_consensus %>%
  select(id_case, x_coord_to = utm_easting, y_coord_to = utm_northing, id_progen, 
         chain_id = lineage_chain, date_infectious = symptoms_started) %>%
  mutate(month = floor_date(date_infectious, unit = "month")) %>%
  left_join(select(case_dt, id_progen = id, x_coord_from = utm_easting,  
                   y_coord_from = utm_northing, date_exposed = symptoms_started)) %>%
  ungroup() %>%
  mutate(group = 1:nrow(.), 
         date_exposed = case_when(is.na(date_exposed) ~ ymd(date_infectious), 
                                  TRUE ~ ymd(date_exposed)), 
         x_coord_from = ifelse(is.na(x_coord_from), x_coord_to, x_coord_from), 
         y_coord_from = ifelse(is.na(y_coord_from), y_coord_to, y_coord_from)) -> chain_segs

bez_pts <- get_bezier_pts(from = as.matrix(select(chain_segs, lat = y_coord_from, long = x_coord_from)),
                          to = as.matrix(select(chain_segs, lat = y_coord_to, long = x_coord_to)),
                          frac = 0.25,
                          transform = function(x) log(1 / x) * 500)

bez_pts <- left_join(chain_segs, bez_pts, by = "group")

# trying it here
bez_pts %>% 
  mutate(date_infectious = floor_date(date_infectious, unit = "month"), 
         date_exposed = floor_date(date_exposed, unit = "month")) -> bez_pts
top_chain_hist$chain_id <- top_chain_hist$lineage_chain

# sampled
bez_pts$sampled <- bez_pts$id_case %in% links_consensus$id_case[links_consensus$lineage != 0]

# save as a movie (more manageable file size and res than gif/pdf)
system.time(
  saveVideo({
    make_gif(bez_pts, top_chain_hist, cols_chains,
             pemba_shp, 
             date_seqs = seq.Date(from = ymd("2010-01-01"), 
                                  to = ymd("2018-10-01"), 
                                  by = "month"), 
             fade_out = 1)  
  }, 
  img.name = 'trans_plot', title = 'Reconstructed transmission trees', 
  description = 
    c('Reconstructed transmission trees using lognormal SI and weibull dk'),
  interval = 0.5, nmax = 50, ani.dev = "png", ani.type = "png",
  ani.res = 800, ani.width = 800, ani.height = 800))
