# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]


################################################################################

ranks = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species')
for (r in ranks) {
  altExp(tse,r) <- agglomerateByRank(tse, r, agglomerate.tree = TRUE)
  altExp(tse_bacteria,r) <- agglomerateByRank(tse_bacteria, r, agglomerate.tree = TRUE)
}

altExps(tse) <-
  lapply(altExps(tse),
         function(y){
           rowData(y)$prevalence <- 
             getPrevalence(y, detection = 1/100, 
                           sort = FALSE,
                           assay.type = "counts", 
                           as.relative = TRUE)
           y
         })


top_phyla <- getTopFeatures(
  altExp(tse,"Phylum"),
  method="sum",
  top=10L,
  assay.type="counts")


top_phyla_mean <-  getTopFeatures(altExp(tse,"Phylum"),
                         method="mean",
                         top=5L,
                         assay.type="counts")

x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addHierarchyTree(x)


prevalence_tree = plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")


png(filename="figures/prevalence_tree.png" ,units = 'in',width=9, height=6, res=1000)
print(prevalence_tree)
dev.off()

################################################################################
