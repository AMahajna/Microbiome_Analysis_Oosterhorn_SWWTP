## Plot species accumulation curve
ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
for (r in ranks) {
  altExp(tse,r) <- agglomerateByRank(tse, r, agglomerate.tree = TRUE)
}

#Calculate species accumulation
species_accum <- specaccum(t(assay(tse)), method = 'random')
genera_accum <- specaccum(t(assay(altExp(tse,"Genus"))), method = 'random')
family_accum <- specaccum(t(assay(altExp(tse,"Family"))), method = 'random')

bacterial_species_accum <- specaccum(t(assay(tse_bacteria)), method = 'random')

#active_species_accum <- specaccum(t(assay(tse_active)), method = 'random')

ranks = c('Kingdom', 'Phylum', 'Class', 'Order')
for (r in ranks) {
  altExp(tse_pathway,r) <- agglomerateByRank(tse_pathway, r, agglomerate.tree = TRUE)
}
#enzyme_accum <- specaccum(t(assay(tse_pathway)), method = 'random')
pathway_accum <- specaccum(t(assay(altExp(tse_pathway,"Class"))), method = 'random')

png(filename="figures/accum_curve.png" ,units = 'in',width=9, height=6, res=1000)
# Plot species accumulation curve
plot(species_accum, xlab = "Number of Samples", ylab = "Number of Features",
     col = "skyblue", lwd = 3,  xlim = c(0, 35), ci = 0, ylim = c(0,1500))

# Add the second curve using the lines() function
lines(
  genera_accum,
  col = "darkgreen",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  family_accum,
  col = "darkorange",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  bacterial_species_accum,
  col = "slategrey",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add the second curve using the lines() function
lines(
  pathway_accum,
  col = "firebrick",      # Different color for the second curve
  lwd = 3, # Line width
  ci = 0
)

# Add a legend to the plot
legend(
  "topleft",
  legend = c("Metabolic Pathways","Bacterial species","Species","Genera","Families"),  # Labels for each curve
  col = c("firebrick","slategrey","skyblue","darkgreen","darkorange"),         # Corresponding colors
  lwd = 3,
  cex = 1,# Line width in the legend,
  bty = "n"
)
dev.off()

