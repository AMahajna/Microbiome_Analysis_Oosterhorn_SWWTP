prevalence_plot_labels <- prevalence_plot_labels +
  theme(
    axis.text.x = element_text(
      angle = 90,   # Rotate x-axis labels
      vjust = 1,    # Vertical justification
      hjust = 1,    # Horizontal justification
      size = 6      # Adjust font size (smaller)
    )
  )

removal_plot <- removal_plot +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)  # Rotate x-axis labels
  )

# Combine plots with cowplot
combined_plot_RE_prevalence <- plot_grid(
  removal_plot, prevalence_plot_labels, 
  labels = c("A", "B"), # Labels for the plots
  nrow = 2,             # Correct argument: nrow (not nrows)
  label_size = 14       # Size of the labels
)

png("figures/combined_plot_RE_prevalence.png", width = 8, height = 6, units = "in", res = 1000)

tiff("figures/combined_plot_RE_prevalence.tiff", width = 8, height = 6, units = "in", res = 1000)
print(combined_plot_RE_prevalence)
dev.off()

