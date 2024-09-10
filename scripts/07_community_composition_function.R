##Diversity heatmap 

#set taxonomy ranks 
colnames(rowData(tse_pathway))<- c("SEED_subsystem_functional_category", 
                                   "description", 
                                   "SEED_subsystem",
                                   "enzyme")
setTaxonomyRanks(colnames(rowData(tse_pathway)))
ranks =getTaxonomyRanks()

for (r in ranks) {
  altExp(tse_pathway,r) <- agglomerateByRank(tse_pathway, r, agglomerate.tree = TRUE)
}
#length(unique(rowData(tse_pathway)[ ,1]))

tse_functional_category <- agglomerateByRank(tse_pathway,rank = "SEED_subsystem_functional_category")

# Add clr-transformation on samples
tse_functional_category <- transformAssay(tse_functional_category, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)

# Add standardize-transformation on features (taxa)
tse_functional_category <- transformAssay(tse_functional_category, assay.type = "clr",
                                          MARGIN = "features", 
                                          method = "standardize", name = "clr_z")

# Gets the assay table
mat <- assay(tse_functional_category, "clr_z")


#png(filename="figures/heatmap_functional_category.png" ,units = 'in',width=9, height=6, res=1000)
# Creates the heatmap
pheatmap(mat)
#dev.off()



core_description = getPrevalentFeatures(altExp(tse_pathway,"description"), detection = 0, prevalence = 99/100,
                                        rank = "description", sort = TRUE)

#check
sum((getPrevalence(altExp(tse_pathway,"description"), detection = 1/100, prevalence =1, 
                   sort = TRUE, assay.type = "counts",
                   as.relative = TRUE))== 1)


core_subsystem = getPrevalentFeatures(altExp(tse_pathway,"SEED_subsystem"), detection = 0, prevalence = 99/100,
                                      rank = "SEED_subsystem", sort = TRUE)


core_enzyme = getPrevalentFeatures(tse_pathway, detection = 0, prevalence = 0.999,
                                   rank = "enzyme", sort = TRUE)

sum((getPrevalence(tse_pathway, tank ="enzyme", detection = 1/100, prevalence =1,  
                   sort = TRUE, assay.type = "counts",
                   as.relative = TRUE))== 1)

tse_enzyme <- agglomerateByRank(tse_pathway,rank = "enzyme")

tse_enzyme  <- transformAssay(tse_enzyme , method = "relabundance")

# Add clr-transformation on samples
tse_enzyme <- transformAssay(tse_enzyme, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)

# Add standardize-transformation on features (taxa)
tse_enzyme <- transformAssay(tse_enzyme, assay.type = "clr",
                             MARGIN = "features", 
                             method = "standardize", name = "clr_z")



tse_enzyme_subset <- tse_enzyme[core_enzyme, ]

# Add clr-transformation
tse_enzyme_subset <- transformAssay(tse_enzyme_subset, method = "clr",
                                    MARGIN="samples",
                                    assay.type = "counts")
# Does standardize-transformation
tse_enzyme_subset  <- transformAssay(tse_enzyme_subset , assay.type = "clr",
                                     MARGIN = "features", 
                                     method = "standardize", name = "clr_z")

# Gets the assay table

mat <- assay(tse_enzyme_subset, "clr_z")

png(filename="figures/heatmap_core_enzyme.png" ,units = 'in',width=9, height=6, res=1000)
# Creates the heatmap
pheatmap(mat)
dev.off()

# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

# to view the tree, run
taxa_tree
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters

# Adds information to rowData
rowData(tse_enzyme_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_enzyme_subset))), ]

# Prints taxa and their clusters
rowData(tse_enzyme_subset)$clusters

# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 



# Order data based on 
tse_enzyme_subset <- tse_enzyme_subset[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- unfactor(colData(tse_enzyme_subset)$Season)

############
# Use left_join to merge df1 and df2 by the shared column 'Date'
sample_data$Date = rownames(sample_data)
df2 = as.data.frame(colData(tse_enzyme_subset)[ , c("Date","Season")])
sample_data$Date = as.Date(sample_data$Date)
sample_data <- sample_data %>%
  left_join(df2, by = "Date")
rownames(sample_data) = sample_data$Date
sample_data$Date = NULL
############
png(filename="figures/heatmap_core_enzyme_clusters.png" ,units = 'in',width=9, height=6, res=1000)

breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors, border_color = grey)
dev.off()
