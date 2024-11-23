#SAMSA2 data 
#tse_active 
#tse_pathway 

# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]

#accumulation curve 
species_acc = specaccum(t(assay(tse_active)), method = 'exact')
species_function_curve = as.data.frame(species_acc[[4]])

enzyme_acc = specaccum(t(assay(tse_enzymes)), method = 'exact')
species_function_curve = cbind(species_function_curve, as.data.frame(enzyme_acc[[4]]))

tse_subsystem <- agglomerateByRank(tse_pathway, rank = "Class")
subsystem_acc = specaccum(t(assay(tse_subsystem)), method = 'exact')
species_function_curve = cbind(species_function_curve, as.data.frame(subsystem_acc[[4]]))

tse_functional_category <- agglomerateByRank(tse_pathway, rank = "Phylum")
functional_category_acc = specaccum(t(assay(tse_functional_category)), method = 'exact')
species_function_curve = cbind(species_function_curve, as.data.frame(functional_category_acc[[4]]))

tse_description <- agglomerateByRank(tse_pathway, rank = "Kingdom")
description_acc = specaccum(t(assay(tse_description)), method = 'exact')
species_function_curve = cbind(species_function_curve, as.data.frame(description_acc[[4]]))


colnames(species_function_curve)[1] <- "Species_Accumulation"
colnames(species_function_curve)[2] <- "Enzyme_Accumulation"
colnames(species_function_curve)[3] <- "SEED_Subsystem_Accumulation"
colnames(species_function_curve)[4] <- "Description_Accumulation"
colnames(species_function_curve)[5] <- "SEED_subsystem_functional_category_Accumulation"
# Create a new row to add
new_row <- data.frame(Species_Accumulation = 0, Enzyme_Accumulation= 0,SEED_Subsystem_Accumulation =0,
                      Description_Accumulation= 0, 
                      SEED_subsystem_functional_category_Accumulation=0)  # New values for the first row

# Add the new row at the top of the data frame and shift the rest of the rows down
species_function_curve <- rbind(new_row, species_function_curve[1:(nrow(species_function_curve)), ])

###################################################
function_acc_curve = ggplot() +
  # Scatter plot for species vs. SEED subsystem functional category
  geom_point(data = species_function_curve, aes(x = Species_Accumulation, y = SEED_subsystem_functional_category_Accumulation, color = SEED_subsystem_functional_category_Accumulation), color = "darkgreen") +
  
  # Regression line for species vs. SEED subsystem functional category
  geom_smooth(data = species_function_curve, aes(x = Species_Accumulation, y = SEED_subsystem_functional_category_Accumulation), method = "auto", se = FALSE, color = "palegreen") +
  
  # Regression line for genera accumulation
  # geom_smooth(data = species_function_curve, aes(x = Species_Accumulation, y = Enzyme_Accumulation), method = "auto", se = FALSE, color = "darkgreen") +
  
  # Scatter plot for species vs. SEED subsystem functional category
  geom_point(data = species_function_curve, aes(x = Species_Accumulation, y = SEED_Subsystem_Accumulation), color = "firebrick") +
  
  # Regression line for family accumulation
  geom_smooth(data = species_function_curve, aes(x = Species_Accumulation, y = SEED_Subsystem_Accumulation), method = "auto", se = FALSE, color = "lightcoral") +
  
  # Scatter plot for species vs. SEED subsystem functional category
  geom_point(data = species_function_curve, aes(x = Species_Accumulation, y = Description_Accumulation), color = "blue") +
  
  # Regression line for bacterial species accumulation
  geom_smooth(data = species_function_curve, aes(x = Species_Accumulation, y = Description_Accumulation), method = "auto", se = FALSE, color = "skyblue") +
  
  # Theme and labels
  theme_minimal() +
  labs(x = "Observed Species Richness", y = "Observed Function Richness")

########################################
enzyme_acc_curve = ggplot() +
  # Scatter plot for species vs. SEED subsystem functional category
  geom_point(data = species_function_curve, aes(x = Species_Accumulation, y = Enzyme_Accumulation, color = SEED_subsystem_functional_category_Accumulation), color = "darkorange") +
  
  # Regression line for species vs. SEED subsystem functional category
  geom_smooth(data = species_function_curve, aes(x = Species_Accumulation, y = Enzyme_Accumulation), method = "auto", se = FALSE, color = "lightsalmon") +
  
  # Theme and labels
  theme_minimal() +
  labs(x = "Observed Species Richness", y = "Observed Enzyme Richness") 

###########


# Create a dummy dataframe for the legend
legend_data <- data.frame(
  x = rep(1, 5),
  y = 1:5,
  color = c("black","darkorange", "firebrick", "blue", "darkgreen"),
  label = c("Legend","Enzymes", "SEED Subsystems", "Description", "SEED Subsystems Functional Category")
)

# Create a plot for the manual legend
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = color)) +
  geom_line(aes(group = label), size = 1) +  # Lines for each label
  scale_color_manual(values = c("black","darkorange", "firebrick", "blue", "darkgreen")) +  # Colors for each line
  scale_x_continuous(limits = c(0, 1)) +  # Adjust x axis to make space for lines
  scale_y_continuous(limits = c(0, 1)) +  # Adjust y axis to make space for lines
  theme_void() +  # Remove axis labels
  theme(legend.position = "none") +  # Remove default legend
  labs(color = "Legend")  # Add legend title


# Manually adjust the labels to be placed next to each other horizontally
legend_plot <- legend_plot +
  annotation_custom(grob = textGrob("Legend", gp = gpar(col = "black", fontsize = 10)), xmin = -0.95, ymin = 0.5) +
  annotation_custom(grob = textGrob("Enzymes", gp = gpar(col = "darkorange", fontsize = 10)), xmin = -0.65, ymin = 0.5) +
  annotation_custom(grob = textGrob("SEED Subsystems", gp = gpar(col = "firebrick", fontsize = 10)), xmin = -0.2, ymin = 0.5) +
  annotation_custom(grob = textGrob("Description", gp = gpar(col = "blue", fontsize = 10)), xmin = 0.15, ymin = 0.5) +
  annotation_custom(grob = textGrob("SEED Subsystems Functional Category", gp = gpar(col = "darkgreen", fontsize = 10)), xmin = 0.65, ymin = 0.5)

# Display the legend plot
legend_plot
###################
#Print all 

combined_plot_acc_curve <- plot_grid(
  enzyme_acc_curve, 
  function_acc_curve, 
  legend_plot,
  ncol = 1, 
  labels = c("A", "B"),  # Labels for the plots
  label_size = 10,       # Customize label size
  rel_heights = c(4, 4, 1)  # Adjust relative heights (first two = full size, last = 1/4 size)
)

#png(filename="figures/species_function_acc_curve.png" ,units = 'in',width=9, height=6, res=1000)
#print(combined_plot_acc_curve)
#dev.off()



