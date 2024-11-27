# Directory creation
directories <- c("input_data", "input_data/func_result", "input_data/org_result", "input_data/subsystems_result", 
                 "output_data", "figures", "scripts")

# Create directories if they don't exist
lapply(directories, function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
})

# List of required packages
required_packages <- c("data.table", "vegan", "BiocManager", "biomformat", "mia", "miaTime", "miaViz", "ggtree", 
                       "scuttle", "ggplot2", "readr", "readxl", "plyr", "dplyr", "tidyverse", "remotes", "lubridate", 
                       "pheatmap", "ape", "forcats", "scater", "phyloseq", "ggpubr", "patchwork", "stringr", 
                       "ComplexHeatmap", "shadowtext", "car", "bluster", "kableExtra","grid","gridExtra","cowplot",
                       "RColorBrewer","writexl", "caret","FeatureTerminatoR")

# Install and load packages
lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("mia", "miaTime", "miaViz", "ggtree", "scuttle", "phyloseq", "ComplexHeatmap")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
})