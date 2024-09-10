source(file = "scripts/01_install_load_packages.R")
source(file = "scripts/02_import_clean_data.R")
#source(file = "scripts/03_feature_accumulation_curve.R")
source(file = "scripts/04_exploration_quality_control_taxonomy.R")
source(file = "scripts/05_exploration_quality_control_function.R")




################################################################################
#Diversity  
#Visualizing significance in group-wise comparisons

index <- "shannon"
group_var <- "Season"

tse_bacteria <- estimateDiversity(tse_bacteria, assay.type = "counts", index = "shannon")


# Calculate p values
pvals <- pairwise.wilcox.test(
  tse_bacteria[[index]], tse_bacteria[[group_var]], p.adjust.method = "fdr")
# Put them to data.frame format
pvals <- pvals[["p.value"]] |>
  as.data.frame()
varname <- "group1"
pvals[[varname]] <- rownames(pvals)
# To long format
pvals <- reshape(
  pvals,
  direction = "long",
  varying = colnames(pvals)[ !colnames(pvals) %in% varname ],
  times = colnames(pvals)[ !colnames(pvals) %in% varname ],
  v.names = "p",
  timevar = "group2",
  idvar = "group1"
) |>
  na.omit()
# Add y-axis position
pvals[["y.position"]] <- apply(pvals, 1, function(x){
  temp1 <- tse[[index]][ tse[[group_var]] == x[["group1"]] ]
  temp2 <- tse[[index]][ tse[[group_var]] == x[["group2"]] ]
  temp <- max( c(temp1, temp2) )
  return(temp)
})
pvals[["y.position"]] <- max(pvals[["y.position"]]) +
  order(pvals[["y.position"]]) * 0.2
# Round values
pvals[["p"]] <- round(pvals[["p"]], 3)

# Create a boxplot
p <- plotColData(
  tse_bacteria, x = group_var, y = index, show_violin = TRUE, shape_by = "Season",
  point_size = 5, color_by = "Season") +
  theme(text = element_text(size = 10)) +
  stat_pvalue_manual(pvals) 

#png(filename="figures/diversity_plot.png" ,units = 'in',width=9, height=6, res=1000)
p
#dev.off()

################################################################################

#ls("package:mia")
#co-abundant groups as CAGs, which are clusters of taxa that co-vary across samples 
#getUniqueTaxa(tse_bacteria, rank = "Phylum")
