# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 12/14/2023
# Description: Colors and themes for the lung spatial ASE plots
# ==============================================================================

# ======================================
# Import libraries
# ======================================

library(tidyverse)
library(googlesheets4)
library(ggthemes)
library(ggplot2)
library(ggrepel)

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




