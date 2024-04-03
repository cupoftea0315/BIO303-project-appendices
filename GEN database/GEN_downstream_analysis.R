library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)

# Read data
rm(list=ls())
load("GEN_after.Rdata")

analyze_DESeq2 <- function(data_frame) {
  
  # Gets the name of the data frame
  data_name <- deparse(substitute(data_frame))
  
  # Preprocessing: Remove genes with negative values and round non-integer values
  data_frame[data_frame < 0] <- NA  # Set negative values to NA
  data_frame <- data_frame[complete.cases(data_frame), ]  # Remove rows with NA values
  data_frame <- round(data_frame)  # Round non-integer values to the nearest integer
  
  # Volcano plot
  dds <- DESeqDataSetFromMatrix(countData = data_frame, colData = data.frame(group=group), design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Extract differentiation genes
  diff_genes <- rownames(res)[which(res$padj < 0.05)]
  
  # Create a data frame containing the logFC and p-value for each gene
  df <- data.frame(gene_name = rownames(res), logFC = res$log2FoldChange, P.Value = res$pvalue)
  
  # Mark differentially expressed genes and calculate -log10(p-value)
  df$group <- ifelse(df$logFC>=3&df$P.Value<=0.01, "Up", 
                     ifelse(df$logFC<=-3&df$P.Value<=0.01, "Down", "Not sig"))
  df$pvalue_log10 <- (-log10(df$P.Value))
  volcano_plot <- ggplot(df, aes(x=logFC, y=pvalue_log10)) + 
    geom_point(aes(color=group)) +
    scale_color_manual(values=c("dodgerblue", "gray", "firebrick")) +
    geom_vline(xintercept = 3, linetype = "dashed") +
    geom_vline(xintercept = -3, linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_label_repel(data=df[df$group!="Not sig",],
                     aes(label=gene_name), 
                     box.padding = unit(0.35, "lines"),
                     size = 2) +
    labs(title = "Volcano plot", y=expression(-log10), x="log(Fold Change)") + 
    theme_test()
  ggsave(paste0("GEN_Volcano_plot(", data_name, ").png"), plot = volcano_plot)
  
  # Check the differential expression genes function for up-regulated genes
  up_genes <- subset(df$gene_name, df$group == 'Up')
  keys_up <- NULL
  if (length(up_genes) > 0) {
    up_genes <- sapply(strsplit(up_genes, "_"), "[", 1)
    keys_up <- bitr(up_genes, fromType = "GENEID", toType = "SYMBOL", OrgDb = EnsDb.Hsapiens.v79)
    ego_up <- enrichGO(gene = keys_up$SYMBOL,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  }
  
  # Check the differential expression genes function for down-regulated genes
  down_genes <- subset(df$gene_name, df$group == 'Down')
  keys_down <- NULL
  if (length(down_genes) > 0) {
    down_genes <- sapply(strsplit(down_genes, "_"), "[", 1)
    keys_down <- bitr(down_genes, fromType = "GENEID", toType = "SYMBOL", OrgDb = EnsDb.Hsapiens.v79)
    ego_down <- enrichGO(gene = keys_down$SYMBOL,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
  }
  
  # Export up-regulated differential expression genes names and functions
  if (!is.null(keys_up)) {
    up_genes_names <- as.data.frame(keys_up)
    output_path_up_names <- paste0("GEN_differential_express_genes_up(", data_name, ").csv")
    write.csv(up_genes_names, file = output_path_up_names, row.names = FALSE)
    sorted_ego_up <- ego_up[order(ego_up$p.adjust), ]
    top_up_functions <- head(sorted_ego_up$Description, 20)
  } else {
    top_up_functions <- character(0)
  }
  
  # Export down-regulated differential expression genes names and functions
  if (!is.null(keys_down)) {
    down_genes_names <- as.data.frame(keys_down)
    output_path_down_names <- paste0("GEN_differential_express_genes_down(", data_name, ").csv")
    write.csv(down_genes_names, file = output_path_down_names, row.names = FALSE)
    sorted_ego_down <- ego_down[order(ego_down$p.adjust), ]
    top_down_functions <- head(sorted_ego_down$Description, 20)
  } else {
    top_down_functions <- character(0)
  }
  
  # Heatmap
  top_genes <- head(rownames(res)[order(res$pvalue)], 50)
  pheatmap(data_frame[top_genes, ], fontsize_row = 8, filename = paste0("GEN_Heatmap(", data_name, ").png"))
  
  # GO enrichment analysis
  gene_names <- sapply(strsplit(rownames(data_frame), "_"), "[", 1)
  keys <- bitr(gene_names, fromType = "GENEID", toType = "SYMBOL", OrgDb = EnsDb.Hsapiens.v79)
  ego <- enrichGO(gene = keys$SYMBOL,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  plot <- dotplot(ego, showCategory=30, font.size = 5)
  plot <- plot + theme(axis.text.y = element_text(hjust = 1))
  ggsave(paste0("GEN_GO_enrichment_analysis_dotplot(", data_name, ").png"), plot = plot, width = 10, height = 6, units = "in")
}

analyze_DESeq2(rawdata_combined)
analyze_DESeq2(rawdata_combined_ComBat)
analyze_DESeq2(rawdata_combined_ComBatseq)