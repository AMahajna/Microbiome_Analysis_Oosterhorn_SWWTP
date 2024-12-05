################################################################################
# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

core_pathways = getPrevalence(
  tse_pathway, rank = "Class",detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) %>% head(32)

selected <- rowData(tse_pathway)$Class %in% rownames(as.data.frame(core_pathways)) &
  !is.na(rowData(tse_pathway)$Class)
tse_core_pathways <- tse_pathway[selected, ]

tse_core_pathways = agglomerateByRank(tse_core_pathways, rank = "Class")
tse_core_pathways = transformAssay(
  x = tse_core_pathways, assay.type = "relabundance", method = "clr", name = "clr")

df_core_pathways = t(assay(tse_core_pathways, "clr"))

df = as.data.frame(cbind(df_core_pathways, RE))

control_LOOCV = trainControl( method="LOOCV", returnResamp = 'all')
control_CV = trainControl( method = "cv", number = 5)

################################################################################

set.seed(1100)

N_df = df[, c(1:31, which(names(df) == "N_removal"))]

N_rf_model <- train(
  N_removal ~ ., data = N_df,
  method = "rf",
  trControl = control_LOOCV,
  importance = TRUE
)

# Get feature importance
N_importance <- varImp(N_rf_model, scale = TRUE)
N_importance_plot = plot(N_importance, main = "A) Variable Importance Plot for Predicting Nitrogenous Pollutants Removal Efficiency")

# Plot feature importance
png(filename="figures/feature_importance_N.png" ,units = 'in',width=9, height=6, res=1000)
print(N_importance_plot)
dev.off()
################################################################################

BOD_df = df[, c(1:31, which(names(df) == "BOD_removal"))]

BOD_rf_model <- train(
  BOD_removal ~ ., data = BOD_df,
  method = "rf",
  trControl = control_LOOCV,
  importance = TRUE
)

# Get feature importance
BOD_importance <- varImp(BOD_rf_model, scale = TRUE)
BOD_importance_plot = plot(BOD_importance, main = "B) Variable Importance Plot for Predicting BOD Removal Efficiency")

# Plot feature importance
png(filename="figures/feature_importance_BOD.png" ,units = 'in',width=9, height=6, res=1000)
print(BOD_importance_plot)
dev.off()
################################################################################
P_df = df[, c(1:31, which(names(df) == "P_removal"))]

P_rf_model <- train(
  P_removal ~ ., data = P_df,
  method = "rf",
  trControl = control_LOOCV,
  importance = TRUE
)

# Get feature importance
P_importance <- varImp(P_rf_model, scale = TRUE)
P_importance_plot = plot(P_importance, main = "C) Variable Importance Plot for Predicting Phosphorus Pollutants Removal Efficiency")

# Plot feature importance
png(filename="figures/feature_importance_P.png" ,units = 'in',width=9, height=6, res=1000)
print(P_importance_plot)
dev.off()
################################################################################
png(filename="figures/importance.png" ,units = 'in',width=11, height=16, res=1000)
grid.arrange(N_importance_plot, BOD_importance_plot, P_importance_plot)
dev.off()

################################################################################
################################################################################

RE_long <- pivot_longer(RE, 
                        cols = c("N_removal","BOD_removal","P_removal"), 
                        names_to = "quality_type", 
                        values_to = "value")

#Create the plot with colorblind-friendly colors

removal_plot= ggplot(RE_long, aes(x = Date, y = value, color = quality_type, group = quality_type)) +
  geom_line() +                  # Add lines
  geom_point() +                 # Add dots
  theme_minimal() +              # Use a minimal theme
  labs(x = "Date", y = "Removal Efficiency [%]", color = "Legend") + # Labels
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + # Tilt x-axis labels
  scale_color_brewer(palette = "Set1") # Colorblind-friendly palette

png(filename="figures/removal_efficiency.png" ,units = 'in',width=9, height=6, res=1000)
print(removal_plot)
dev.off()

################################################################################

