################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("input_data/func_result")){dir.create("input_data/func_result")}
if(!dir.exists("input_data/org_result")){dir.create("input_data/org_result")}
if(!dir.exists("input_data/subsystems_result")){dir.create("input_data/subsystems_result")}

if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 
# Check if a package is installed and load it

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}
library(vegan)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

if (!requireNamespace("biomformat", quietly = TRUE)) {
  install.packages("biomformat")
}
library(biomformat)

if (!requireNamespace("mia", quietly = TRUE)) {
  BiocManager::install("mia")
}
library(mia)

if (!requireNamespace("miaTime", quietly = TRUE)) {
  BiocManager::install("miaTime")
}
library(miaTime)

if (!requireNamespace("miaViz", quietly = TRUE)) {
  BiocManager::install("miaViz")
}
library(miaViz)

if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}
library(ggtree)

if (!requireNamespace("scuttle", quietly = TRUE)) {
  BiocManager::install("scuttle")
}
library(scuttle)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  BiocManager::install("ggplot2")
}
library(ggplot2)

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}
library(plyr)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)

if (!requireNamespace("lubridate", quietly = TRUE)) {
  install.packages("lubridate")
}
library(lubridate)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
library(ape)

if (!requireNamespace("forcats", quietly = TRUE)) {
  install.packages("forcats")
}
library(forcats)

if (!requireNamespace("scater", quietly = TRUE)) {
  install.packages("scater")
}
library(scater)

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}
library(phyloseq)


if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install("ggpubr")
}
library(ggpubr)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install("patchwork")
}
library(patchwork)

if (!requireNamespace("stringr", quietly = TRUE)) {
  install("stringr")
}
library(stringr)

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

if (!requireNamespace("shadowtext", quietly = TRUE)) {
  BiocManager::install("shadowtext")
}
library(shadowtext)
