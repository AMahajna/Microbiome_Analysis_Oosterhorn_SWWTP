# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

color_vector <- rainbow(ncol(assay(tse)))

# Open the PNG device
png(filename="figures/rarefaction_curve_species.png", units = 'in', width = 9, height = 6, res = 1000)

# Set up the plot margins to create space for the legend outside the plot
par(mar = c(5, 5, 5, 10))  # Increase the right margin to make space for the legend

# Plot the rarefaction curve
rarecurve(
  t(assay(tse)), 
  step = 500, 
  col = color_vector, 
  label = TRUE, 
  xlim = c(0, 35000), 
  ylim = c(0, 700), 
  lwd = 2
)

# Add the legend outside of the plot area
legend(
  "topright", 
  legend = colnames(assay(tse)), 
  col = color_vector, 
  lty = 1, 
  lwd = 2, 
  cex = 0.5,        # Reduce text size
  ncol = 1,         # Use a single column for the legend
  inset = c(-0.1, 0),  # Adjust legend position outside the plot
  text.width = 10,   # Increase the text width to ensure it fits
  box.lwd = 0,     # Adjust border width (optional)
  box.col = "white",
  text.font = 2,     # Bold text
  x.intersp = 0.1, 
  y.intersp = 1.25, 
  xpd = TRUE         # Allow legend to be drawn outside the plot area
)

# Close the PNG device
dev.off()

################################################################################
ranks = c('Kingdom', 'Phylum', 'Class', 'Order')
for (r in ranks) {
  altExp(tse_pathway,r) <- agglomerateByRank(tse_pathway, r, agglomerate.tree = TRUE)
}

color_vector <- rainbow(ncol(assay(altExp(tse_pathway,"Class"))))

# Open the PNG device
png(filename="figures/rarefaction_curve_functions.png", units = 'in', width = 9, height = 6, res = 1000)

# Set up the plot margins to create space for the legend outside the plot
par(mar = c(5, 5, 5, 10))  # Increase the right margin to make space for the legend

# Plot the rarefaction curve
rarecurve(
  t(assay(altExp(tse_pathway,"Class"))), 
  step = 500, 
  col = color_vector, 
  label = FALSE, 
  xlim = c(0, 25000), 
  ylim = c(0, 600), 
  lwd = 2, 
  ylab = "Cumulative Metabolic Pathways Count"
)

# Add the legend outside of the plot area
legend(
  "topright", 
  legend = colnames(assay(altExp(tse_pathway,"Class"))), 
  col = color_vector, 
  lwd = 2, 
  lty = 1, 
  cex = 0.5,        # Reduce text size
  ncol = 1,         # Use a single column for the legend
  inset = c(-0.1, 0),  # Adjust legend position outside the plot
  text.width = 10,   # Increase the text width to ensure it fits
  box.lwd = 0,     # Adjust border width (optional)
  box.col = "white",
  text.font = 2,     # Bold text
  x.intersp = 0.1, 
  y.intersp = 1.25, 
  xpd = TRUE         # Allow legend to be drawn outside the plot area
)

# Close the PNG device
dev.off()

################################################################################
