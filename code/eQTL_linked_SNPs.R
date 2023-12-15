#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 3/22/2023
# Description: Finding eSNPs in LD with SNPs overlapping exons
#==============================================================================

#==============================================================================
# Load libraries
#==============================================================================

suppressMessages({library(lifecycle)
                  library(readJDX)
                  library(VariantAnnotation)
                  library(vcfR)
                  library(data.table)
                  library(snpStats)
                  library(tidyverse)
                  library(LDheatmap)
                  library(plyr)
                  library(pheatmap)
                  library(RColorBrewer)
                  library(googlesheets4)
                  library(AnnotationDbi)
                  library(org.Hs.eg.db)
                  library(BSgenome.Hsapiens.NCBI.GRCh38)
                  library(GenomicFeatures)
                  library(ComplexHeatmap)
                  library(viridis)})

#==============================================================================
# Environment variables
#==============================================================================

sinkall(filename = paste0("/scratch/hnatri/ILD/eQTL_linked_snps.Rout"))
set.seed(1234)

#==============================================================================
# Import data
#==============================================================================

# Genotype data
vcf <- readVcf("/labs/banovich/IPF/eQTL/qtl_mapping_vcf/filtered_MAF-HWE-INDPW.vcf", verbose = FALSE)
var_ranges <- rowRanges(vcf)
gt <- geno(vcf)
gt <- gt@listData$GT
snp_matrix <- genotypeToSnpMatrix(vcf)

rm(vcf)

# FFPE blocks in the IPF sheet
# https://docs.google.com/spreadsheets/d/1rCrE8KAelXHX_MHxi48119ZLmJhs29QZb1eGcc11mg0/edit?usp=sharing
gs4_deauth()
ipf_metadata <- gs4_get("https://docs.google.com/spreadsheets/d/1rCrE8KAelXHX_MHxi48119ZLmJhs29QZb1eGcc11mg0/edit?usp=sharing")
sheet_names(ipf_metadata)
ffpe_info <- read_sheet(ipf_metadata, sheet = "FFPE Blocks")

length(unique(ffpe_info$`Sample Name`))
length(intersect(ffpe_info$`Sample Name`, rownames(snp_matrix$genotypes)))
intersect(ffpe_info$`Sample Name`, rownames(snp_matrix$genotypes))
setdiff(ffpe_info$`Sample Name`, rownames(snp_matrix$genotypes))
sort(setdiff(rownames(snp_matrix$genotypes), ffpe_info$`Sample Name`))

rownames(snp_matrix$genotypes) <- gsub("VUHD71", "VUHD071", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD72", "VUHD072", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD73", "VUHD073", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD74", "VUHD074", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD78", "VUHD078", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD75", "VUHD075", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD94", "VUHD094", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD92", "VUHD092", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD95", "VUHD095", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD65", "VUHD065", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD67", "VUHD067", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD66", "VUHD066", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD68", "VUHD068", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD70", "VUHD070", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD84", "VUHD084", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD92", "VUHD092", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD69", "VUHD069", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD85", "VUHD085", rownames(snp_matrix$genotypes))
rownames(snp_matrix$genotypes) <- gsub("VUHD98", "VUHD098", rownames(snp_matrix$genotypes))

ffpe_info$`Sample Name` <- gsub("VUHD71", "VUHD071", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD72", "VUHD072", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD73", "VUHD073", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD74", "VUHD074", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD78", "VUHD078", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD75", "VUHD075", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD94", "VUHD094", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD92", "VUHD092", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD95", "VUHD095", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD65", "VUHD065", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD67", "VUHD067", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD66", "VUHD066", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD68", "VUHD068", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD70", "VUHD070", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD84", "VUHD084", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD92", "VUHD092", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD69", "VUHD069", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD85", "VUHD085", ffpe_info$`Sample Name`)
ffpe_info$`Sample Name` <- gsub("VUHD98", "VUHD098", ffpe_info$`Sample Name`)

ffpe_info[which(ffpe_info$`Sample Name` %in% setdiff(ffpe_info$`Sample Name`, rownames(snp_matrix$genotypes))),]

gt_ffpe_samples <- intersect(ffpe_info$`Sample Name`, rownames(snp_matrix$genotypes))

#var_info <- as.data.frame(var_info)
#gt_data <- as.data.frame(gt_data)

#==============================================================================
# eQTL and gene info
#==============================================================================

# All significant eQTL
mashr_snp_info <- fread("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/filtered_MAF-HWE-INDPW_snp-info.txt", header=T)
mashr_sighits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-significant-eQTL.rds")
mashr_tophits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-top-eQTL.rds")

# Interaction eQTL
inteqtl_topsnp_sig <- read.table("/labs/banovich/IPF/eQTL/inteQTL_topSNPs.tsv", sep = "\t", header = T)

esnps <- unique(c(mashr_sighits$rowData$variant_id, inteqtl_topsnp_sig$rsid))
egenes <- unique(mashr_sighits$rowData$feature_id, inteqtl_topsnp_sig$gene)

# Gene annotations
genes_gtf <- read.table("/labs/banovich/SingleCell/CellRanger/3_1_0/Ensemble_93/PipelineData/Projects/IPF/References/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf", sep = "\t")
head(genes_gtf)
gene_ids <- sapply(strsplit(genes_gtf$V9," "), `[`, 2)
gene_ids <- gsub(";", "", gene_ids)
gene_names <- sapply(strsplit(genes_gtf$V9," "), `[`, 6)
gene_names <- gsub(";", "", gene_names)
genes_gtf$feature_id <- gene_ids
genes_gtf$feature_name <- gene_names
genes_gtf <- dplyr::filter(genes_gtf, V3 == "gene")
genes_gtf <- dplyr::select(genes_gtf, V1, V4, V5, feature_id, feature_name)
colnames(genes_gtf) <- c("chr", "start", "stop", "ensembl", "gene")

# Panel genes
panel_genes <- read.csv("/home/hnatri/ILD_spatial_ASE/hLung_XeniumV1.csv")

#==============================================================================
# Calculating LD
#==============================================================================

# R squared between all pairs of SNPs
genotypes <- snp_matrix$genotypes

# For every eGene, finding eSNPs in LD with SNPs overlapping the
# gene body
result_list <- lapply(egenes, function(gene){
  message(gene)
  # gene <- egenes[2]
  gene_info <- genes_gtf[which(genes_gtf$ensembl == gene),]
  snps <- mashr_snp_info[which(mashr_snp_info$snp_chr==as.numeric(gene_info$chr) & mashr_snp_info$snp_loc < gene_info$stop & mashr_snp_info$snp_loc > gene_info$start),]$rsid
  
  if(length(snps)==0){return(NULL)}

  # LD between these SNPs and eSNPs
  esnps_subset <- mashr_snp_info[which(mashr_snp_info$rsid %in% esnps & mashr_snp_info$snp_chr == as.numeric(gene_info$chr)),]$rsid
  genotypes_subset <- genotypes
  genotypes_subset@.Data <- genotypes_subset@.Data[,which(colnames(genotypes_subset@.Data) %in% unique(c(esnps_subset, snps)))]
  ld <- ld(genotypes_subset, depth = 20000, stats = "R.squared")
  
  # Which eSNPs had gene body SNPs with R^2 >0.8
  length(snps[which(snps %in% colnames(ld))])
  length(snps)
  
  ld <- as.matrix(ld)
  ld <- ld[esnps_subset, snps]
  ld <- as.data.frame(ld)
  ld$esnp <- rownames(ld)
  ld <- pivot_longer(ld, cols = setdiff(colnames(ld), "esnp"), names_to = "gene_body_snp", values_to = "R2")
  ld_sig <- ld[which(ld$R2>0.8),]

  ld_sig$gene <- unique(gene_info$gene)
  ld_sig$ensembl <- unique(gene_info$ensembl)
  ld_sig$chr <- unique(gene_info$chr)
  ld_sig$start <- min(gene_info$start)
  ld_sig$stop <- max(gene_info$stop)
  ld_sig$esnp_loc <- plyr::mapvalues(x = ld_sig$esnp,
                                     from = mashr_snp_info$rsid,
                                     to = mashr_snp_info$snp_loc)

  return(ld_sig)
})

names(result_list) <- egenes
result_list <- plyr::compact(result_list)

#saveRDS(result_list, "/scratch/hnatri/ILD/eSNP_LD_spatial.rds")
result_df <- rbindlist(result_list, use.names=TRUE, fill=TRUE)
write.table(result_df, "/scratch/hnatri/ILD/eSNP_LD_spatial.tsv", quote = F, row.names = F, sep = "\t")

sinkall()
q(save = "no")

result_df <- read.table("/scratch/hnatri/ILD/eSNP_LD_spatial.tsv", sep = "\t", header = T)
result_df$gene_snp <- paste0(result_df$gene, "|", result_df$esnp)

#==============================================================================
# Visualizing basic metrics
#==============================================================================

# Distribution of # linked loci
hist(table(result_df$gene_snp), breaks = 100, xlab = "Gene body SNPs with LD>0.8", main = "")

# How many eQTL had at least two linked loci overlapping the gene?
as.data.frame(table(result_df$gene_snp)) %>% filter(Freq>=2) %>% dim()

as.data.frame(table(result_df$gene_snp)) %>%
  filter(Freq>10) %>%
  ggplot(aes(x = reorder(Var1, (-Freq)), y = Freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

as.data.frame(table(as.data.frame(table(result_df$gene_snp))$Freq)) %>%
  ggplot(aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Matrix of genes with >2 linked SNPs, heterozygosity across all individuals.
# For each sample, calculating mean dosage across all SNPs overlapping each of
# the 3677 genes, or just the eSNP
keep_genes <- as.data.frame(table(result_df$gene_snp)) %>% filter(Freq>=1)
dim(keep_genes)

keep_genes <- separate(data = keep_genes, col = Var1, into = c("gene", "snp"), sep = "\\|")
#keep_snps <- c()

head(gt)
gt_data <- gt[which(rownames(gt) %in% result_df[which(result_df$gene %in% keep_genes$gene),]$esnp),]
#rm(gt)
#gt_data <- gt[which(rownames(gt) %in% keep_snps),]
gt_data <- gsub("1/1", 2, gt_data)
gt_data <- gsub("0/1", 1, gt_data)
gt_data <- gsub("1/0", 1, gt_data)
gt_data <- gsub("0/0", 0, gt_data)
gt_data_num <- matrix(as.numeric(gt_data),
                  ncol = ncol(gt_data))
colnames(gt_data_num) <- colnames(gt_data)
rownames(gt_data_num) <- rownames(gt_data)

# Metadata
eqtl_info <- as.data.frame(mashr_sighits$rowData)
eqtl_info$gene_symbol <- mapIds(org.Hs.eg.db, keys = eqtl_info$feature_id, keytype = "ENSEMBL", column = "SYMBOL")
eqtl_info$gene_snp <- paste0(eqtl_info$gene_symbol, "|", eqtl_info$variant_id)

result_df$type <- mapvalues(x = result_df$gene_snp,
                                 from = eqtl_info$gene_snp,
                                 to = eqtl_info$type)
result_df$type[-which(result_df$gene_snp %in% eqtl_info$gene_snp)] <- NA

# Plot data
result_df_plot <- result_df[which(result_df$esnp %in% rownames(gt_data_num)),]
gt_data_num_plot <- gt_data_num[which(rownames(gt_data_num) %in% result_df_plot[,"esnp"]),]
result_df_plot <- result_df[which(result_df$esnp %in% rownames(gt_data_num_plot)),]
result_df_plot <- distinct(dplyr::select(result_df_plot, -c("gene_body_snp", "R2", "ensembl", "start", "stop")))
result_df_plot <- result_df_plot[!duplicated(result_df_plot$esnp),]
result_df_plot <- result_df_plot[match(rownames(gt_data_num_plot), result_df_plot$esnp),]

annotation <- data.frame(esnp=as.character(result_df_plot$esnp), eQTL_type="type") %>%
  column_to_rownames("esnp")
annotation$eQTL_type <- mapvalues(x = rownames(annotation),
                                  from = result_df$esnp,
                                  to = result_df$type)

p <- pheatmap(gt_data_num_plot, scale="row", annotation_row = annotation,
              show_rownames = T,
              annotation_colors = list(eQTL_type=c("global"="darkgreen",
                                                   "multi-state"="orange",
                                                   "unique"="tomato3",
                                                   "int-eQTL"="darkblue")),
              color = colorRampPalette(c("navy", "white", "red"))(50),
              cutree_cols=3, cutree_rows=3)

filename <- "/home/hnatri/ILD_spatial_ASE/dosage_prettyheatmap.pdf"
pdf(file = filename,
    width = 12,
    height = 10)
p
dev.off()

png(file="/home/hnatri/ILD_spatial_ASE/dosage_prettyheatmap.png",
    width=1200, height=1000)
pheatmap(gt_data_num_plot, scale="row", annotation_row = annotation,
         show_rownames = F,
         color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()

# Same with ComplexHeatmap

# Color
type_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(setdiff(unique(result_df$type), NA)))
names(type_col) <- setdiff(unique(result_df$type), NA)

row_ha <- rowAnnotation(
  df = data.frame(type = result_df_plot$type),
  annotation_height = unit(4, "mm"),
  col = list("type" = type_col)
)

plot_func <- function(){
  hm <- Heatmap(gt_data_num_plot, 
                name = "Dosage", # Title of legend
                col = viridis(100),
                use_raster = T,
                column_title = NULL, row_title = NULL,
                show_column_names = TRUE,
                cluster_columns = TRUE,
                show_row_names = FALSE,
                cluster_rows = TRUE,
                row_km = 5,
                right_annotation = row_ha,
                height = nrow(gt_data_num_plot)*unit(0.007, "mm"),
                column_names_gp = gpar(fontsize = 6))  # Text size for row names
  
  heatmap <- draw(hm)
}

p <- plot_func()

# Saving to a file
filename <- "/home/hnatri/ILD_spatial_ASE/dosage_complexheatmap.pdf"
pdf(file = filename,
    width = 15,
    height = 10)
p
dev.off()

#png(file="/home/hnatri/ILD_spatial_ASE/dosage_complexheatmap.png",
#    width=1200, height=1000)
#p
#dev.off()

# Extracting cluster info
r.dend <- row_dend(p)  # Extract row dendrogram
rcl.list <- row_order(p)  #E xtract clusters (output is a list)

# Loop to extract SNPs for each cluster
for (i in 1:length(row_order(p))){
  if (i == 1) {
  clu <- t(t(row.names(gt_data_num_plot[row_order(p)[[i]],])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("esnp", "Cluster")
  } else {
    clu <- t(t(row.names(gt_data_num_plot[row_order(p)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
    }
}

head(out)

result_df_plot <- merge(result_df_plot, out, by = "esnp")
result_df_info <- merge(result_df_plot, eqtl_info, by = "gene_snp")

unique(result_df_plot[which(result_df_plot$type=="global"),]$nSignificant)
unique(result_df_info[which(result_df_info$type.y=="global"),]$nSignificant)
unique(result_df_info[which(result_df_info$eQTL_type=="global"),]$nSignificant)

result_df_info$type <- result_df_info$type.y

result_df_info <- distinct(dplyr::select(result_df_info, -c("type.y", "type.x", "sc_eQTL")))
result_df_info$int_eQTL <- ifelse(result_df_info$gene_snp %in% inteqtl_topsnp_sig$gene_rsid, "TRUE", "FALSE")
result_df_info$base_panel <- ifelse(result_df_info$gene %in% panel_genes$gene_name, "TRUE", "FALSE")

#write.table(result_df_info, "/home/hnatri/ILD_spatial_ASE/SNPs_selected_genes.tsv", sep = "\t", quote = F, row.names = F)
#result_df_info <- read.table("/home/hnatri/ILD_spatial_ASE/SNPs_selected_genes.tsv", sep = "\t", header = T)

#==============================================================================
# Adding manually selected genes
#==============================================================================

# Manually selected genes
# https://docs.google.com/spreadsheets/d/1f7HtQC07PWa1ij9UxTmgrEMX2OnHDVcrpn7Oa_VNIms/edit?usp=sharing
candidates <- gs4_get("https://docs.google.com/spreadsheets/d/1f7HtQC07PWa1ij9UxTmgrEMX2OnHDVcrpn7Oa_VNIms/edit?usp=sharing")
sheet_names(candidates)
candidates <- read_sheet(candidates, sheet = "Updated candidates")
result_df_candidates <- result_df[which(result_df$gene_snp %in% candidates$gene_snp),]

# Adding linked SNP location
varinfo_df <- as.data.frame(var_ranges@ranges)
varinfo_df$chr <- as.data.frame(as.data.frame(var_ranges@seqnames))[,"value"]
varinfo_df$ref <- as.character(var_ranges$REF)
varinfo_df$alt <- as.data.frame(var_ranges$ALT)[,"value"]

colnames(varinfo_df) <- gsub("names", "gene_body_snp", colnames(varinfo_df))
result_df_candidates <- merge(result_df_candidates, varinfo_df, by = "gene_body_snp")
result_df_candidates$panel <- ifelse(result_df_candidates$gene %in% panel_genes$gene_name, "TRUE", "FALSE")
table(result_df_candidates$panel)

#write.table(result_df_candidates, "/home/hnatri/ILD_spatial_ASE/result_df_candidates.tsv", sep = "\t", quote = F, row.names = F)
#result_df_candidates <- read.table("/home/hnatri/ILD_spatial_ASE/result_df_candidates.tsv", sep = "\t", header = T)

prop_table <- as.data.frame(table(result_df_info[,"Cluster"], as.character(result_df_info[,"type"])))
colnames(prop_table) <- c("Cluster", "Type", "Freq")
prop_table <- spread(prop_table, Cluster, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, Cluster, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)

ggplot(prop_table, aes(x=Cluster, y=Freq, fill = Type)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw()

# Adding effect sizes for each eQTL in each cell type
effects <- as.data.frame(mashr_sighits$posterior_means)
effects$gene_snp <- rownames(mashr_sighits$posterior_means)
result_df_info$gene_snp <- paste0(result_df_info$feature_id, "|", result_df_info$rsid)

# # of linked SNPs
keep_genes$gene_snp <- paste0(keep_genes$gene, "|", keep_genes$snp)
colnames(keep_genes) <- c("gene", "snp", "n_linked_snps", "gene_snp")
result_df_info$gene_snp <- paste0(result_df_info$gene_symbol, "|", result_df_info$esnp)
result_df_info <- merge(result_df_info, keep_genes, by = "gene_snp", all = F)

#write.table(result_df_info, "/home/hnatri/ILD_spatial_ASE/SNPs_selected.tsv", sep = "\t", quote = F, row.names = F)
#result_df_info <- read.table("/home/hnatri/ILD_spatial_ASE/SNPs_selected.tsv", sep = "\t", header = T)
  
# Summary of sample genotypes for selected loci
result_df_candidates <- read.table("/home/hnatri/ILD_spatial_ASE/result_df_candidates.tsv", sep = "\t", header = T)
ffpe_info_tgen <- ffpe_info[which(ffpe_info$Location=="tgen"),]

gt_data_ffpe <- gt_data[,which(colnames(gt_data) %in% ffpe_info_tgen$`Sample Name`)]
gt_data_ffpe <- gt_data_ffpe[which(rownames(gt_data_ffpe) %in% c(result_df_candidates$esnp, result_df_candidates$gene_body_snp)),]

#write.table(result_df_candidates, "/scratch/hnatri/ILD/ASE_candidate_SNP_gts_FFPE_samples.tsv", quote = F, sep = "\t", row.names = F)
#result_df_candidates <- read.table("/scratch/hnatri/ILD/ASE_candidate_SNP_gts_FFPE_samples.tsv", header = T, sep = "\t")

# Sequences 50bp upstream and downstream
#gr <- GRanges(seqnames = paste0('chr',seq(1,3)),IRanges(start=c(1,100000,100000),end=c(2,100001,100001)))
result_df_candidates$region_start <- result_df_candidates$start.y-50
result_df_candidates$region_stop <- result_df_candidates$start.y+50

result_df_candidates$gene_chr_linkedsnploc_ref_alt <- paste0(result_df_candidates$gene, "_", result_df_candidates$chr.x, "_", result_df_candidates$start.y, "_", result_df_candidates$ref, "_", result_df_candidates$alt)

gr <- makeGRangesFromDataFrame(result_df_candidates[,-which(colnames(result_df_candidates) %in% c("start", "end", "stop"))],
                               keep.extra.columns=T,
                               ignore.strand=T,
                               seqinfo=NULL,
                               seqnames.field="chr.x",
                               start.field="region_start",
                               end.field="region_stop",
                               #strand.field="strand",
                               starts.in.df.are.0based=F)
tmp <- Views(Hsapiens, gr)
sequence <- DNAStringSet(as.character(tmp))

names(sequence) <- result_df_candidates$gene_chr_linkedsnploc_ref_alt

writeXStringSet(sequence, format="fasta", file="/scratch/hnatri/ILD/ASE_candidate_SNP_sequences.fasta")

