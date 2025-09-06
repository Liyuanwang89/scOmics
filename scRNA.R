################################
########-snRNA pipeline-########
################################
# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggsci)

# Set color palettes
colpalettes <- unique(c(
  pal_npg("nrc")(10), pal_aaas("default")(10), pal_nejm("default")(8),
  pal_lancet("lanonc")(9), pal_jama("default")(7), pal_jco("default")(10),
  pal_ucscgb("default")(26), pal_d3("category10")(10), pal_locuszoom("default")(7),
  pal_igv("default")(51), pal_uchicago("default")(9), pal_startrek("uniform")(7),
  pal_tron("legacy")(7), pal_futurama("planetexpress")(12), pal_rickandmorty("schwifty")(12),
  pal_simpsons("springfield")(16), pal_gsea("default")(12)
))

# Set working directory
setwd("/home/wangliyuan/scRNA_datas/combine_merge")

# Define a function to process each dataset
process_obj <- function(data.dir, project.name) {
  # Read data and create Seurat object
  data <- Read10X(data.dir = data.dir)
  obj <- CreateSeuratObject(
    counts = data,
    project = project.name,
    min.cells = 3,
    min.features = 200
  )
  
  # Calculate mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  # Subset data (fixed the subset issue from original code)
  obj <- subset(obj, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 20)
  obj <- runDecontX(obj)  
  obj <- doubletFinder_v3(obj, 
                                   PCs = 1:50, 
                                   pN = 0.25, 
                                   pK = 0.09, 
                                   nExp = round(ncol(obj) * 0.05))
  # Normalization and feature selection
  obj <- NormalizeData(object = obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:10)
  
  return(obj)
}

# Define all tissues and time points
tissues <- c(
  "bladder", "bone_marrow", "cerebellum", "cerebrum", "dorsal_muscle",
  "duodenum", "fat", "heart", "ileum", "jejunum", "kidney", "leg_muscle",
  "liver", "lung", "lymph", "ovary", "pancreas", "PBMC", "spleen",
  "stomach", "trachea"
)

time_points <- c("E65", "D1", "Y1", "Y3")

# Initialize list to store all objects
all_objects <- list()

# Process all datasets using nested loops
for (tissue in tissues) {
  for (time in time_points) {
    # Skip combinations that don't exist (like ovary_E65 which wasn't in original code)
    if (tissue == "lymph" && time == "E65") next
    if (tissue == "ovary" && time == "E65") next
    if (tissue == "PBMC" && time == "E65") next
    
    # Construct path and project name
    path <- paste0(
      "/home/wangliyuan/scRNA_datas/scRNA_integ_all_202204/sc_cellranger_data/",
      tissue, "_", time
    )
    proj_name <- paste0(tissue, "_", time)
    
    # Check if directory exists before processing
    if (dir.exists(path)) {
      cat("Processing:", proj_name, "\n")
      obj <- process_obj(path, proj_name)
      all_objects[[proj_name]] <- obj
    }
  }
}

# Integration pipeline
features <- SelectIntegrationFeatures(object.list = all_objects)
obj.anchors <- FindIntegrationAnchors(object.list = all_objects, anchor.features = features)
obj.integrated <- IntegrateData(anchorset = obj.anchors)

# Set default assay and save
DefaultAssay(object = obj.integrated) <- "integrated"
saveRDS(obj.integrated, file = "IntegrateData.rds")

# Scale, PCA, and UMAP
obj.integrated <- ScaleData(object = obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(object = obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", dims = 1:30)

# Clustering
obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 40, dims = 1:30)
obj.integrated <- FindClusters(object = obj.integrated, resolution = 0.5)
obj.integrated <- RunTSNE(object = obj.integrated, dims = 1:30, seed.use = 2)

# Save final object
DefaultAssay(obj.integrated) <- "RNA"
saveRDS(obj.integrated, file = "obj.integrated_clusters.rds")

# Find markers and save results
inputmarkers <- FindAllMarkers(obj.integrated)
write.table(inputmarkers, "markers_results.txt", quote = FALSE, sep = "\t")

# Save final integrated object
saveRDS(obj.integrated, file = "obj.integrated.rds")

