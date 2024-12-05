################################################################################
# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

################################################################################
##Prevalence  
#Create tse_phylum 
tse_phylum <- agglomerateByRank(tse, rank = "Phylum", update.tree = TRUE)
tse_class <- agglomerateByRank(tse, rank = "Class", update.tree = TRUE)
altExp(tse, "Phylum") <- tse_phylum
altExp(tse, "Class") <- tse_class

#Prevalence of Phylum in total community labeled by kingdom 
rowData(altExp(tse,"Phylum"))$prevalence <- 
  getPrevalence(altExp(tse,"Phylum"), detection = 0, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

prevalence_plot = plotRowData(altExp(tse,"Phylum"), "prevalence", point_size=5,
                              colour_by = "Kingdom")



rowData(altExp(tse,"Class"))$prevalence <- 
  getPrevalence(altExp(tse,"Class"), detection = 0, 
                sort = FALSE,
                assay.type = "counts", as.relative = TRUE)

df_prvalence_class = as.data.frame(rowData(altExp(tse,"Class"))$prevalence)

df_prevalence = as.data.frame(rowData(altExp(tse,"Phylum"))$prevalence)

df_core_phyla = df_prevalence[df_prevalence == 1]

#png(filename="figures/prevalence.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_plot)
#dev.off()

df_prevalence <- rowData(altExp(tse, "Phylum"))[ ,c("Kingdom","Phylum","prevalence")]

# Assuming df_prevalence is your dataset
prevalence_plot_labels = ggplot(df_prevalence, aes(x = Phylum, y = prevalence, color = Kingdom)) +
  geom_point(size = 4, alpha = 0.7) +  # Add alpha for transparency to avoid overplotting
  scale_color_brewer(palette = "Set1") +  # Use a nice color palette for distinction
  labs(x = "Phylum", y = "Prevalence") +  # Add a title
  theme_minimal() +  # Use minimal theme for a clean look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
    axis.title = element_text(size = 12),  # Adjust axis title size
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Style the title
  )

png(filename="figures/prevalence_phyla.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_plot)
dev.off()

png(filename="figures/prevalence_phyla_labels.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_plot_labels)
dev.off()

#Core community is made of 8 Phyla (5 bacteria and 3 Eukaryota)
sum(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as.relative = TRUE, prevalence = 1) == 1)


################################################################################

core_species = getPrevalence(
  tse, rank = "Species",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(15)

cat("There are 14 core species in the microbiome")
print(core_species)

df_core_species <- data.frame(
  Species = names(core_species),
  Value = as.numeric(core_species)
)

# Subset by feature
selected <- rowData(tse)$Species %in% df_core_species$Species &
  !is.na(rowData(tse)$Species)
tse_sub <- tse[selected, ]

#write_xlsx(as.data.frame(rowData(tse_sub)), "output_data/core_species.xlsx")
################################################################################


