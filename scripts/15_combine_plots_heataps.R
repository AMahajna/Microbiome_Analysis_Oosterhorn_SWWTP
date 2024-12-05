# Convert the heatmap to a ggplot object
heatmap1_ggplot <- as.ggplot(heatmap_func)
heatmap2_ggplot <- as.ggplot(p)

# Combine plots using plot_grid
combined_plot_heatmap <- plot_grid(
  heatmap1_ggplot, heatmap2_ggplot,
  labels = c("A", "B"),  # Labels for the plots
  nrow = 2,              # Correct argument: nrow (not nrows)
  label_size = 14        # Size of the labels
)

png("figures/combined_plot_heatmap.png", width = 8, height = 12, units = "in", res = 1000)
print(combined_plot_heatmap)
dev.off()
