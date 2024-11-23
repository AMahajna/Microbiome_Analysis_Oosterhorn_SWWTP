#Which function from pathway most contributes to removal efficiency 
#consider each one alone (BOD, Sulphates, and TN alone and TRE all together)
#Use RFE with all possible base models in their, and see if you can average the
# featuer importance to see which function most contributes to the removal efficiency 


# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]

################################################################################
#Adjust this:
#Define the core functions which you want to relate to the removal efficiency 


###############################################################################
###############################################################################
##Recursive Feature Elimination  

#rfe1 is fitted to Random Forest model with with 10 k-fold cross correlation  
#ref2 is fitted to bagged trees model with with leave one out cross correlation
#Weekly_Data <- weekly_data

#13 input and diversity column 14 is output 
#normalization of the variable values and splitting of input and target values  
x <-Weekly_Data[,1:13]
normalization <- preProcess(x)
x <- predict(normalization, x)
x <- as.data.frame(x)
y<- as.data.frame(diversity)

#training scheme: setting up the controls for each recursive feature elimination 
control_RF_CV = rfeControl(functions=rfFuncs, method="cv", repeats = 5, number = 10, returnResamp = 'all')
control_RF_LOOCV = rfeControl(functions=rfFuncs, method="LOOCV", returnResamp = 'all')
control_TB_LOOCV = rfeControl(functions=treebagFuncs, method="LOOCV", returnResamp = 'all')

#reproducible data 
set.seed(121321)

#split data- 80% for training and 20% for testing 
inTrain <- createDataPartition(Weekly_Data$diversity, p= .80, list = FALSE)[,1]

x_train <- x[ inTrain, ]
x_test <- x[-inTrain, ]

y_train <- y[ inTrain,1]
y_test <- y[ -inTrain,1]

#run RFE
results_rfe1 <- rfe(x =x_train , y= y_train , sizes=c(1:13),rfeControl=control_RF_CV)
sprintf("The optimal number of variables is: %s", results_rfe1$bestSubset)
sprintf("Optimal Variable: %s", results_rfe1$optVariables)

# RFE LOOCV
results_rfe_rf_loocv <- rfe(x =x_train , y= y_train , sizes=c(1:13),rfeControl=control_RF_LOOCV)
sprintf("The optimal number of variables is: %s", results_rfe_rf_loocv$bestSubset)
sprintf("Optimal Variable: %s", results_rfe_rf_loocv$optVariables)

results_rfe2 <- rfe(x =x_train , y= y_train , sizes=c(1:13),rfeControl=control_TB_LOOCV)
sprintf("The optimal number of variables is: %s", results_rfe2$bestSubset)
sprintf("Optimal Variable: %s", results_rfe2$optVariables)
#plot <- ggplot(data = results_rfe1, metric = "RMSE")

png(filename="figures/Variables_RMSE.png", units ='in', height=8, width=8, res = 1000)
trellis.par.set(caretTheme())
plot1 <- plot(results_rfe1, type = c("g", "o"))
plot2 <- xyplot(results_rfe1, 
                type = c("g", "p"), 
                ylab = "RMSE CV Estimates")
print(plot1, split=c(1,1,1,2), more=TRUE)
print(plot2, split=c(1,2,1,2))
dev.off()

#plot for RMSE of rfe2
plot(results_rfe2, type = c("g", "o"))

#plot variable importance for random forest with cross validation 
png(filename="figures/Var_Imp.png", units ='in', height=5, width=5, res = 1000)
varimp_data <- data.frame(feature = row.names(varImp(results_rfe1))[1:3],
                          importance = varImp(results_rfe1)[1:3, 1])
plot_var_imp = ggplot(data = varimp_data, 
                      aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.25, hjust=1.25,color="white", size=4) + 
  theme_bw() + theme(legend.position = "none", 
                     axis.text.y = element_text(angle = 45, hjust = 1, size =11)) +
  coord_flip()
print(plot_var_imp)
dev.off()

# Combine plots and label them
#combined_plot <- plot_grid(plot1, plot_var_imp, labels = c("a", "b"), ncol = 1)
# Create text labels for the plots
label_a <- textGrob("a", gp = gpar(fontsize = 16), x = 0, y = 1, just = c("left", "top"))  # Position in top-left corner
label_b <- textGrob("b", gp = gpar(fontsize = 16), x = 0, y = 1, just = c("left", "top"))  # Position in top-left corner


# Combine plots and add labels
combined_plot <- grid.arrange(
  arrangeGrob(plot1, top = label_a),  # Add label 'a' above the first plot
  arrangeGrob(plot_var_imp, top = label_b),  # Add label 'b' above the second plot
  ncol = 1
)
# Save the combined plot as a PNG
ggsave("figures/combined_plot_RMSE_var_imp.png", combined_plot, height=8, width=8)

#plot variable importance for bagged trees with "leave-one-out cross-validation" 
png(filename="figures/Var_Imp_LOOCV.png", units ='in', height=4, width=4, res = 1000)
varimp_data <- data.frame(feature = row.names(varImp(results_rfe2))[1],
                          importance = varImp(results_rfe2)[1, 1])
ggplot(data = varimp_data, 
       aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.25, hjust=1.25,color="white", size=4) + 
  theme_bw() + theme(legend.position = "none") +
  coord_flip()
dev.off()

#plot the density distribution of the RMSE for number of variables included in model
png(filename="figures/Density.png", units ='in', height=5, width=8, res = 1000)
#plot5 = stripplot(results_rfe1)
#plot6 = histogram(results_rfe1)
plot7 = densityplot(results_rfe1)
print(plot7)
dev.off()




