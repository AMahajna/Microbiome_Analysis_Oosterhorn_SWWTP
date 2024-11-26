# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]

################################################################################
#All subsystems classes of stress 
# Stress Response 
#selected_stress <- rowData(tse_pathway)$Kingdom %in% c("Stress Response") &
#  !is.na(rowData(tse_pathway)$Kingdom)
#tse_stress <- tse_pathway[selected_stress, ]

#tse_stress = agglomerateByRank(tse_stress, rank = "Phylum")

#tse_stress = transformAssay(
#  x = tse_stress, assay.type = "relabundance", method = "clr", pseudocount = TRUE,
#  name = "clr")

#environmental data is standardized and abundance data is CLR transformed
#df_env_stress = cbind(as.data.frame(colData(tse_stress))[ ,10:22],t(as.data.frame(assay(tse_stress, "relabundance"))))
################################################################################

tse_func_cat = agglomerateByRank(tse_pathway, rank = "Kingdom")
tse_func_cat = transformAssay(
  x = tse_func_cat, assay.type = "relabundance", method = "clr", pseudocount = TRUE,
  name = "clr")

df_env_func_cat = cbind(as.data.frame(colData(tse_func_cat))[ ,10:22],t(as.data.frame(assay(tse_func_cat, "relabundance"))))

################################################################################
set.seed(102020023)

control_LOOCV = trainControl( method="LOOCV", returnResamp = 'all')
#control_CV = trainControl( method = "cv", number = 5)

y <- df_env_func_cat[ ,"Stress Response"]
X = df_env_func_cat[ ,1:13]
data_combined = cbind(X,y)

rf_model <- train(
  y ~ ., data = data_combined,
  method = "rf",
  trControl = control_LOOCV,
  importance = TRUE
)

# Get feature importance
importance <- varImp(rf_model)
print(importance)

# Plot feature importance
#png(filename="figures/feature_importance_rf_loocv_relabundance_stress_response.png" ,units = 'in',width=9, height=6, res=1000)
plot(importance)
#dev.off()





