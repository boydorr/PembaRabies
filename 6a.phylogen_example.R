
# Example phylogenetic reconcile
examp <- data.frame(from = c(1, 1, 3, 4, 5, 5, 6, 6, 7),
                    to = c(2, 3, 4, 5, 6, 7, 8, 9, 10),
                    prob = c(0.25, 0.75, 0.4, 0.9, 0.1, 0.95, 0.5, 0.9, 0.55))

lineages <- data.frame(id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       lineage = c("A", "", "", "B", "", "", "", "A", "", "B"))
lineages$lineage[lineages$lineage == ""] <- "Unsampled"
gr <- graph_from_data_frame(examp, vertices = lineages)

ggraph(gr) +
  geom_edge_link(aes(alpha = prob), edge_width = 1.2) +
  geom_node_point(aes(shape = lineage, color = lineage,
                      fill = ifelse(lineage == "Unsampled", TRUE, FALSE)), size = 4) +
  scale_shape_manual(values = c(22, 25, 21)) +
  scale_fill_manual(values = c("white", "grey50")) +
  scale_color_manual(values = c("blue", "red", "white"), name = "Lineage") +
  scale_edge_alpha(name = "Scaled probability") +
  theme_graph() +
  guides(color = guide_legend(override.aes = list(fill = c("white", "white", "grey50"), 
                                                  shape = c(22, 25, 21))))

ggsave("figs/example.tiff", dpi = 600, height = 7, width = 6)
