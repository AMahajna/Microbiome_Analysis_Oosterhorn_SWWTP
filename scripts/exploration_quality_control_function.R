################################################################################
##Abundance 

#relative abundance for the top-10 phylum over a log-scaled axis
#png(filename="figures/density_plot_bacterial_phylum.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(altExp(tse_pathway, "Order"), layout = "jitter", 
                     assay.type = "relabundance",
                     n = 10, point_size=2, point_shape=19, colour_by="Season",
                     point_alpha=0.5) + 
  scale_x_log10(label=scales::percent)
#dev.off()

#relative abundance for the top-10 species over a log-scaled axis
#png(filename="figures/density_plot_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
plotAbundanceDensity(altExp(tse_bacteria, "Species"), layout = "density", 
                     assay.type = "relabundance",
                     n = 10, colour_by="Season", 
                     point_alpha=1/10) +
  scale_x_log10(label=scales::percent)
#dev.off()

################################################################################
##Core community 

core_phylum = getPrevalence(
  tse, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(8)

core_species =getPrevalence(
  tse, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(14)

core_phylum_bacteria = getPrevalence(
  tse_bacteria, rank = "Phylum",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(5)

core_species_bacteria = getPrevalence(
  tse_bacteria, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(13)

core_class_bacteria = getPrevalence(
  tse_bacteria, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(12)

#mapTaxonomy(tse_bacteria, taxa = "s__flagellatus")
#Homonyms in Taxonomy 

################################################################################
##Prevalence  
#Creating alternative experiment for each taxonomic level

#Prevalence of Phylum in total community labeled by kingdom 
rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

#png(filename="figures/prevalence.png" ,units = 'in',width=9, height=6, res=1000)
plotRowData(altExp(tse,"Phylum"), "prevalence", point_size=5,
            colour_by = "Kingdom")
#dev.off()

#Core community is made of 8 Phyla (5 bacteria and 3 Eukaryota)
sum(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as.relative = TRUE, prevalence = 1) == 1)
################################################################################
#Prevalence taxonomic tree

#tse_bacteria <- tse[
#  rowData(tse)$Kingdom %in% c("k__Bacteria"), ]

altExps(tse_bacteria) <- lapply(
  altExps(tse_bacteria), function(y){
    rowData(y)$prevalence <- getPrevalence(
      y, detection = 1/100,
      sort = FALSE,
      assay.type = "counts",
      as.relative = TRUE)
    return(y)
  })

top_phyla <- getTopFeatures(
  altExp(tse_bacteria,"Phylum"),
  method="sum",
  top=5L,
  assay.type="counts")
top_phyla_mean <- getTopFeatures(
  altExp(tse_bacteria,"Phylum"),
  method="mean",
  top=5L,
  assay.type="counts")

#Number of ranks can also be decreased
x <- unsplitByRanks(tse_bacteria, ranks = taxonomyRanks(tse_bacteria)[1:7])
x <- addHierarchyTree(x)

#png(filename="figures/prevalence_tree_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
plotRowTree(
  x[rowData(x)$Phylum %in% top_phyla,],
  edge_colour_by = "Phylum",
  tip_colour_by = "prevalence",
  node_colour_by = "prevalence")
#dev.off()


################################################################################
#Quality Control 
tse_bacteria <- addPerCellQC(tse_bacteria)

#png(filename="figures/library_size.png" ,units = 'in',width=9, height=6, res=1000)
plotColData(tse_bacteria,"sum","Season", colour_by = "Season", point_size = 5, 
            shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()
################################################################################
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse_bacteria))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis


df_quality <- as.data.frame(colData(tse_bacteria)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df_quality, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#The distribution of calculated library sizes can be visualized as a histogram (left)
#or by sorting the samples by library size (right).

LibrarySize =p1 + p2

png(filename="figures/library_size_bacteria.png" ,units = 'in',width=9, height=6, res=1000)
#pdf("figures/LibrarySize.pdf")
print(LibrarySize)
dev.off() 

