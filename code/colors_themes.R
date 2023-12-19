# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 12/14/2023
# Description: Colors and themes for the lung spatial ASE plots
# ==============================================================================

# ======================================
# Import libraries
# ======================================

suppressPackageStartupMessages({library(tidyverse)
                                library(googlesheets4)
                                library(ggthemes)
                                library(ggplot2)
                                library(ggrepel)})

# ======================================
# Cell type markers
# ======================================

epithelial_features <- c("EGFR", "DUOX1", "NKX2-1", "AGER", "RTKN2", "NAPSA",
                         "PGC", "SFTA2", "SFTPC", "SFTPD", "KRT14", "KRT15",
                         "KRT5", "KRT6A", "S100A2", "TP63", "KRT17", "AGR3",
                         "C20orf85", "FOXJ1", "GCLM", "DMBT1", "EPCAM", "KRT18",
                         "MGST1", "MMP7", "FOXI1", "MUC5B", "SCGB1A1",
                         "SCGB3A2", "WFDC2", "ATF3", "KRT8", "SOX9", "SPINK1",
                         "GKN2", "MMP10", "SOX2", "CDH26", "TP73", "CFTR",
                         "HES1", "PKM", "SOX4", "NUCB2", "RNASE1", "SAA2",
                         "AKR1C1", "AKR1C2", "BPIFA1", "CEACAM5", "ERN2",
                         "FCGBP", "GSR", "LTF", "CCNA1", "ICAM1", "ITGB1")
endothelial_features <- c("APLN", "CA4", "HEY1", "BMPR2", "CD34", "EPAS1",
                          "FCN3", "GNG11", "PECAM1", "APLNR", "COL15A1",
                          "PLVAP", "ACKR1", "POSTN", "CLDN5", "RAMP2", "ZEB1",
                          "HAS1", "KDR", "CDKN2A")
immune_features <- c("PPARG", "BANK1", "CD19", "CD79A", "LTB", "MS4A1",
                     "TNFRSF13C", "CD86", "GZMB", "HLA-DRA", "CCR7", "CXCR4",
                     "PTPRC", "TCL1A", "CD69", "CD4", "CD8A", "CD8B", "CD2",
                     "CD28", "CD3D", "CD3E", "CD3G", "FOXP3", "GZMK", "TRAC",
                     "ITM2C", "CD27", "CCL5", "LCK", "FABP4", "MARCO", "MCEMP1",
                     "SPP1", "FCN1", "S100A12", "S100A8","S100A9", "CCL22",
                     "ITGAM", "NFKB1",  "IFIT2", "FGFBP2", "GNLY", "KLRB1",
                     "KLRC1", "NKG7", "LILRA4", "BCL2", "CD79B", "CXCR5", 
                     "CXCL9", "GPR183", "HLA-DQA1", "KLRG1", "BCL2L11", "CD52", 
                     "SLC1A3", "TNFRSF9", "CTLA4",  "IL2RA", "LAG3", "PDCD1",
                     "PDIA6", "PIM2", "IL7R", "LEF1", "FASLG", "HAVCR2",
                     "ISG20", "CPA3","KIT", "TPSAB1", "C1QC", "CD68", "MS4A7",
                     "AIF1", "CD14", "FCGR3A", "FCER1G", "SLC25A37", "CD247",
                     "GZMA", "IRF7")
mesenchymal_features <- c("MFAP5", "PI16", "SFRP2", "ELN", "FAP", "AXL", "LGR6", 
                          "COL1A1", "COL1A2", "COL3A1", "DCN", "FN1", "HAS2", 
                          "LUM", "MEG3", "SPARCL1", "CTHRC1")

# ======================================
# Colors
# ======================================

samples <- c("TILD126", "TILD010", "VUHD115", "TILD074", "TILD028", "TILD055",
             "TILD041", "THD0026", "TILD001", "TILD062", "TILD093", "TILD103",
             "TILD123", "VUILD104", "VUILD48", "VUILD96", "TILD049", "VUHD113",
             "TILD080", "TILD111", "TILD113", "TILD059", "THD0016", "TILD109",
             "VUILD102", "TILD102", "VUILD78", "TILD136", "TILD006", "TILD030",
             "TILD037", "VUILD59", "VUILD91")
sample_col <- colorRampPalette(tableau_color_pal(palette = "Tableau 20")(20))(length(samples))
names(sample_col) <- samples

clusters <- c(0, seq(1, 20))
cluster_col <- colorRampPalette(gdocs_pal()(10))(length(clusters))
names(cluster_col) <- clusters

sample_type <- c("ILD", "HD")
sample_type_col <- c("violetred", "cadetblue")
names(sample_type_col) <- sample_type




