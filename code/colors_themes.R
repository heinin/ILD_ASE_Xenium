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
                                library(dplyr)
                                library(ggplot2)
                                library(ggrepel)
                                library(viridis)})

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
endothelial_features <- c("APLN", "HEY1", "BMPR2", "CD34", "EPAS1",
                          "FCN3", "GNG11", "PECAM1", "APLNR", "COL15A1",
                          "PLVAP", "ACKR1", "POSTN", "CLDN5", "RAMP2", "ZEB1",
                          "HAS1", "KDR", "CDKN2A", "CCL21", "PROX1", "PDPN",
                          "DKK2", "GJA5", "BMX", "SPRY1", "VWA1", "MYC",
                          "HBEGF", "PLA1A", "CPE", "PTGDS", "CA4",
                          "VIPR1", "RGCC", "CYB5A", "ADGRL2")
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
                          "LUM", "MEG3", "SPARCL1", "CTHRC1", "CLU", "APOD",
                          "FBLN1", "CST3", "IGFBP6", "SCARA5", "CD34",
                          "ACTA2", "PDGFRB", "MYH11", "TAGLN", "DES", "ACTG2",
                          "MYLK", "WNT2", "A2M", "GPC3", "MACF1", "CES1",
                          "LIMCH1", "MSLN", "UPK3B", "HP", "WT1", "RGS5",
                          "HIGD1B", "GJA4")

# ======================================
# ggplot themes
# ======================================

my_theme <- theme(axis.text = element_text(size = 9),
                  axis.title = element_text(size = 9),
                  plot.title = element_text(size = 10))

manuscript_theme <- theme(axis.text = element_text(size = 6),
                          axis.title = element_text(size = 7),
                          plot.title = element_text(size = 7),
                          legend.text = element_text(size = 6),
                          legend.title = element_text(size = 7),
                          strip.text.x = element_text(size = 7))
                          
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
sample_type_col <- c("#a8344e", "#1E5B89")
names(sample_type_col) <- sample_type

# Cell lineage colors
lineage_col <- c("#1E5B89", "#466b53", "#a8344e", "#85497c")
names(lineage_col) <- c("Epithelial", "Immune", "Endothelial", "Mesenchymal")

# Cell type annotations for color names
gs4_deauth()
ct_annot  <- gs4_get("https://docs.google.com/spreadsheets/d/1SDfhxf6SjllxXEtNPf32ZKTEqHC9QJW3BpRYZFhpqFE/edit?usp=sharing")

# Epithelial cell type colors
epi_annot <- read_sheet(ct_annot, sheet = "Epithelial")
#epi_celltype_col <- colorRampPalette(c("lightgreen", "darkgreen"))(length(epi_annot$annotation_1))
#epi_celltype_col <- epi_annot$color
epi_celltype_col <- viridis(length(epi_annot$snn_res.0.8),
                            alpha = 1,
                            begin = 0.4,
                            end = 1,
                            direction = 1,
                            option = "D")
names(epi_celltype_col) <- epi_annot$annotation_1

# Immune cell type colors
immune_annot <- read_sheet(ct_annot, sheet = "Immune")
immune_cts <- setdiff(immune_annot$annotation_1, c("Endothelial",
                                                   "Epithelial",
                                                   "Endo_Mesen",
                                                   "Fibroblast"))
immune_celltype_col <- viridis(length(immune_cts),
                               alpha = 1,
                               begin = 0.4,
                               end = 1,
                               direction = 1,
                               option = "B")
names(immune_celltype_col) <- immune_cts

# Mesenchymal cell type colors
mesen_annot <- read_sheet(ct_annot, sheet = "Mesenchymal")
mesen_cts <- setdiff(mesen_annot$annotation_1, c("Lymphatic_C1"))

mesen_celltype_col <- viridis(length(mesen_cts),
                              alpha = 1,
                              begin = 0.4,
                              end = 1,
                              direction = 1,
                              option = "A")
names(mesen_celltype_col) <- mesen_cts

# Endothelial cell type colors
endo_annot <- read_sheet(ct_annot, sheet = "Endothelial")
endo_cts <- endo_annot$annotation_1
endo_celltype_col <- viridis(length(endo_cts),
                              alpha = 1,
                              begin = 0.4,
                              end = 1,
                              direction = 1,
                              option = "G")
names(endo_celltype_col) <- endo_cts

# All celltypes
all_celltypes_annot <- read_sheet(ct_annot, sheet = "All celltypes, annotated, merged")

epi_cts_annotation_2 <- all_celltypes_annot %>%
  filter(lineage == "Epithelial") %>%
  dplyr::select(annotation_2) %>% unlist() %>% as.character()

imm_cts_annotation_2 <- all_celltypes_annot %>%
  filter(lineage == "Immune") %>%
  dplyr::select(annotation_2) %>% unlist() %>% as.character()

endo_cts_annotation_2 <- all_celltypes_annot %>%
  filter(lineage == "Endothelial") %>%
  dplyr::select(annotation_2) %>% unlist() %>% as.character()

mesen_cts_annotation_2 <- all_celltypes_annot %>%
  filter(lineage == "Mesenchymal") %>%
  dplyr::select(annotation_2) %>% unlist() %>% as.character()

epi_col <- viridis(length(epi_cts_annotation_2),
                   alpha = 1,
                   begin = 0.4,
                   end = 1,
                   direction = 1,
                   option = "D")
names(epi_col) <- epi_cts_annotation_2

imm_col <- viridis(length(imm_cts_annotation_2),
                   alpha = 1,
                   begin = 0.4,
                   end = 1,
                   direction = 1,
                   option = "B")
names(imm_col) <- imm_cts_annotation_2

endo_col <- viridis(length(endo_cts_annotation_2),
                    alpha = 1,
                    begin = 0.4,
                    end = 1,
                    direction = 1,
                    option = "A")
names(endo_col) <- endo_cts_annotation_2

mesen_col <- viridis(length(mesen_cts_annotation_2),
                     alpha = 1,
                     begin = 0.4,
                     end = 1,
                     direction = 1,
                     option = "G")
names(mesen_col) <- mesen_cts_annotation_2

epi_col <- colorRampPalette(c("#B9EAF2", "#1E5B89"))(length(epi_cts_annotation_2))
names(epi_col) <- epi_cts_annotation_2[order(epi_cts_annotation_2)]
imm_col <- colorRampPalette(c("#F3EBA2", "#466b53"))(length(imm_cts_annotation_2))
names(imm_col) <- imm_cts_annotation_2[order(imm_cts_annotation_2)]
endo_col <- colorRampPalette(c("#F4ABB8", "#a8344e"))(length(endo_cts_annotation_2))
names(endo_col) <- endo_cts_annotation_2[order(endo_cts_annotation_2)]
mesen_col <- colorRampPalette(c("#D7D0CD", "#85497c"))(length(mesen_cts_annotation_2))
names(mesen_col) <- mesen_cts_annotation_2[order(mesen_cts_annotation_2)]

annotation_2_col <- c(epi_col, imm_col, endo_col, mesen_col)

epi_cts_annotation_3 <- all_celltypes_annot %>%
  filter(lineage == "Epithelial") %>%
  dplyr::select(annotation_3) %>% unlist() %>% as.character()

imm_cts_annotation_3 <- all_celltypes_annot %>%
  filter(lineage == "Immune") %>%
  dplyr::select(annotation_3) %>% unlist() %>% as.character()

endo_cts_annotation_3 <- all_celltypes_annot %>%
  filter(lineage == "Endothelial") %>%
  dplyr::select(annotation_3) %>% unlist() %>% as.character()

mesen_cts_annotation_3 <- all_celltypes_annot %>%
  filter(lineage == "Mesenchymal") %>%
  dplyr::select(annotation_3) %>% unlist() %>% as.character()

epi_annotation_3_col <- colorRampPalette(c("#B9EAF2", "#1E5B89"))(length(epi_cts_annotation_3))
names(epi_annotation_3_col) <- epi_cts_annotation_3[order(epi_cts_annotation_3)]
imm_annotation_3_col <- colorRampPalette(c("#F3EBA2", "#466b53"))(length(imm_cts_annotation_3))
names(imm_annotation_3_col) <- imm_cts_annotation_3[order(imm_cts_annotation_3)]
endo_annotation_3_col <- colorRampPalette(c("#F4ABB8", "#a8344e"))(length(endo_cts_annotation_3))
names(endo_annotation_3_col) <- endo_cts_annotation_3[order(endo_cts_annotation_3)]
mesen_annotation_3_col <- colorRampPalette(c("#D7D0CD", "#85497c"))(length(mesen_cts_annotation_3))
names(mesen_annotation_3_col) <- mesen_cts_annotation_3[order(mesen_cts_annotation_3)]

annotation_3_col <- c(epi_annotation_3_col, imm_annotation_3_col,
                      endo_annotation_3_col, mesen_annotation_3_col)
