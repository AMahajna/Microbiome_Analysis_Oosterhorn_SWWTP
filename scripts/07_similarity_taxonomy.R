################################################################################
tse_bacteria_rda = altExp(tse_bacteria, "Species")

#Rarefy dat 
#otu_table <- assay(tse_bacteria_rda, "counts")
#min_depth <- min(rowSums(otu_table))
#rarefied_otu_table <- rrarefy(otu_table , min_depth)
#assay(tse_bacteria_rda, "counts") <- rarefied_otu_table

# Calculate the list of sequencing depths across samples
sequencing_depths <- colSums(assay(tse_bacteria_rda))
# Calculate variation between highest and lowest sequencing depth
depth_variation <- max(sequencing_depths)/min(sequencing_depths)
depth_variation

# Perform RDA
#PERMANOVA assumes that the dispersion (variance) of the response variables is 
#similar across groups. Non-homogeneous dispersion can violate this assumption, 
#potentially leading to misleading results.

tse_bacteria_rda <- runRDA(tse_bacteria_rda,
                           assay.type = "relabundance",
                           formula = assay ~Capacity_blowers_.+INF_COD_mg_O2_per_l + T_avg_C,
                           distance = "bray",
                           na.action = na.exclude)

#Full formula
#formula = assay ~ INF_Cl_mg_per_l + INF_COD_mg_O2_per_l + INF_Nkj_mg_N_per_l + INF_PO4o_mg_P_per_l+  INF_SO4_µg_per_l + INF_TSS_mg_per_l + Glycerol_kg +Return_sludge_m3_per_h+ Inf_Flow_m3_per_h + Capacity_blowers_. + DW_AT_g_per_l + SVI_10 + T_avg_C,

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse_bacteria_rda, "RDA"), "significance")

#Print out rda_info
rda_info$permanova |>
  knitr::kable()
#The resulting warning is not critical as we don't aim to generalize the model

#rda_info$homogeneity |>
#  knitr::kable()

# Generate RDA plot coloured by season
rda_plot= plotRDA(tse_bacteria_rda, "RDA", colour_by = "Season")
#png(filename="figures/RDA_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
rda_plot
#dev.off()

#Statistically significant (P < 0.05)
################################################################################
#Visualize dbRDA loadings

# Extract loadings for first eigenvector
rda <- reducedDim(tse_bacteria_rda, "RDA")
rda <- attr(rda, "rda")
coef <- rda$CCA$v
coef <- coef[, 1, drop = FALSE]

# Get the taxa with biggest weights
top_coef <- head(coef[rev(order(abs(coef))), , drop = FALSE], 20)
# Sort weights in increasing order
top_coef <- top_coef[order(top_coef), ]

# Create data.frame
df_rda <- data.frame(
  x = top_coef,
  y = factor(names(top_coef), unique(names(top_coef))))

# Create a plot
#png(filename="figures/RDA_loadings_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
ggplot(df_rda, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  labs(x = "", y= "", title = "Top Taxa") +
  theme_bw()
#dev.off()
################################################################################



