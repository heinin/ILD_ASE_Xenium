# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 12/14/2023
# Description: Utilities
# ==============================================================================

# ======================================
# Import libraries
# ======================================

suppressPackageStartupMessages({library(tidyverse)
                                library(googlesheets4)})

# ======================================
# Helper functions
# ======================================

get_pcs <- function(seurat_object, reduction_name="pca") {
  
  # Determine percent of variation associated with each PC
  pct <- seurat_object[[reduction_name]]@stdev / sum(seurat_object[[reduction_name]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % 
  # variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co1
  
  # Determine the difference between variation of PC and subsequent PC and
  # selecting last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  #pcs <- min(co1, co2)
  
  c(co1, co2)
  
}

recluster <- function(seurat_object){
  # Normalize and scale data and run PCA
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- NormalizeData(seurat_object,
                                 normalization.method = "LogNormalize",
                                 verbose = F)
  seurat_object <- ScaleData(seurat_object,
                             features = rownames(seurat_object),
                             verbose = F)
  seurat_object <- RunPCA(seurat_object,
                          reduction.name = "pca",
                          features = rownames(seurat_object),
                          verbose = F)
  # Get PCs for UMAP
  npcs <- get_pcs(seurat_object)[2]
  message(paste0("# PCs for UMAP: ", npcs))
  # Find neighbors, cluster and UMAP
  seurat_object <- RunUMAP(seurat_object,
                           reduction = "pca",
                           reduction.name = "umap",
                           dims = 1:npcs,
                           return.model = TRUE)
  seurat_object <- FindNeighbors(seurat_object,
                                 reduction = "pca",
                                 dims = 1:npcs,
                                 graph.name = c("nn",
                                                "snn"))
  seurat_object <- FindClusters(seurat_object,
                                resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                graph.name = "snn")
  
  
  seurat_object
}
