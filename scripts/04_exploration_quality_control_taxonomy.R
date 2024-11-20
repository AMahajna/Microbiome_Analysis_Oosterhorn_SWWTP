################################################################################
##Abundance: SAMSA2 Data for active bacterial community 

#relative abundance for the top-20 phylum
#png(filename="figures/abundance_density_plot_active_20.png" ,units = 'in',width=9, height=6, res=1000)
top_20 = plotAbundanceDensity(tse_active, layout = "jitter", 
                     assay.type = "relabundance",
                     n = 20, point_size=2, point_shape=19, colour_by="Season",
                     point_alpha=0.5)
#dev.off()

#relative abundance for the top-5 species
#png(filename="figures/abundance_density_plot_active_5.png" ,units = 'in',width=9, height=6, res=1000)
top_5 = plotAbundanceDensity(tse_active, layout = "density", 
                     assay.type = "relabundance",
                     n = 5, colour_by="Season", 
                     point_alpha=1/10) 
#dev.off()

combined_plot_density <- plot_grid(
  top_20, 
  top_5, 
  ncol = 1, 
  labels = c("A", "B"),  # Labels for the plots
  label_size = 12        # Customize label size
)

#png(filename="figures/combined_abundance_density_plot.png" ,units = 'in',width=9, height=6, res=1000)
print(combined_plot_density)
#dev.off()
################################################################################
##Core community 
core_active_species = getPrevalence(
  tse_active, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(20)

cat("There are 19 core active species in the microbiome")
print(core_active_species)
################################################################################
################################################################################
################################################################################

##Core community from tse (generated using Kraken2 includes all Kingdoms)

################################################################################

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
#Create tse_phylum 
tse_phylum <- agglomerateByRank(tse, rank = "Phylum", update.tree = TRUE)
altExp(tse, "Phylum") <- tse_phylum

#Prevalence of Phylum in total community labeled by kingdom 
rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

prevalence_plot = plotRowData(altExp(tse,"Phylum"), "prevalence", point_size=5,
                              colour_by = "Kingdom")
#png(filename="figures/prevalence.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_plot)
#dev.off()

#Core community is made of 8 Phyla (5 bacteria and 3 Eukaryota)
sum(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as.relative = TRUE, prevalence = 1) == 1)


################################################################################
#Prevalence taxonomic tree

altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")

rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 1/100, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)


top_phyla <- getTopFeatures(
  altExp(tse,"Phylum"),
  method="sum",
  top=5L,
  assay.type="counts")
top_phyla_mean <- getTopFeatures(
  altExp(tse,"Phylum"),
  method="mean",
  top=5L,
  assay.type="counts")

#Number of ranks can also be decreased
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:7])
x <- addHierarchyTree(x)

prevalence_tree = plotRowTree(
  x[rowData(x)$Phylum %in% top_phyla,],
  edge_colour_by = "Phylum",
  tip_colour_by = "prevalence",
  node_colour_by = "prevalence")

#png(filename="figures/prevalence_tree.png" ,units = 'in',width=9, height=6, res=1000)
#print(prevalence_tree)
#dev.off()

################################################################################
################################################################################
#Quality Control 
tse <- addPerCellQC(tse)

library_size_season = plotColData(tse,"sum","Season", colour_by = "Season", point_size = 5, 
            shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#png(filename="figures/library_size_season_kraken2.png" ,units = 'in',width=9, height=6, res=1000)
#print(library_size_season)
#dev.off()
################################################################################
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse))) +
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


df_quality <- as.data.frame(colData(tse)) %>%
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

#png(filename="figures/library_size_kraken2.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize)
#dev.off() 
################################
#Quality Control SAMSA2 
tse_active <- addPerCellQC(tse_active)

library_size_season_active = plotColData(tse_active,"sum","Season", colour_by = "Season", point_size = 5, 
                                           shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

png(filename="figures/library_size_season_active.png" ,units = 'in',width=9, height=6, res=1000)
print(library_size_season_active)
dev.off()
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse_active))) +
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


df_quality <- as.data.frame(colData(tse_active)) %>%
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

LibrarySize_Active =p1 + p2

#png(filename="figures/library_size_samsa2.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize_Active)
#dev.off() 

##################################################
#Quality Control SAMSA2 functions
tse_pathway <- addPerCellQC(tse_pathway)

library_size_season_pathway = plotColData(tse_pathway,"sum","Season", colour_by = "Season", point_size = 5, 
                                         shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#png(filename="figures/library_size_season_pathway.png" ,units = 'in',width=9, height=6, res=1000)
#print(library_size_season_pathway)
#dev.off()
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse_pathway))) +
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


df_quality <- as.data.frame(colData(tse_pathway)) %>%
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

LibrarySize_pathway =p1 + p2

#png(filename="figures/library_size_pathway.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize_pathway)
#dev.off() 

##################################################
