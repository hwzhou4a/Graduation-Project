# ==============================================================================
# RNA-seq Analysis: Symbol-First Aggregation Strategy
# Workflow: Import -> Annotate -> Aggregate to Symbol -> Filter -> DEA/PCA -> Export
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(DESeq2)
  library(biomaRt)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggtext)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(pheatmap)
})

# 设置工作目录
rsem_dir <- "D:/Graduation Project/RNASeq/genome_alignment/03_rsem_quantification"
setwd(rsem_dir)

# ==============================================================================
# Part 1: 数据导入与注释准备
# ==============================================================================
cat("[Part 1] Importing and Annotating...\n")

# 1. 文件列表与元数据
files <- list.files(pattern = "\\.genes\\.results$", recursive = TRUE, full.names = TRUE)
names(files) <- basename(files) %>% str_replace(".genes.results", "")

sample_names <- names(files)
conditions <- str_extract(sample_names, "ZHW-[CK]") %>% str_sub(-1)
short_names <- str_sub(sample_names, -2)

colData <- data.frame(
  row.names = sample_names,
  condition = factor(conditions, levels = c("C", "K")),
  short_name = short_names
)

# 2. 导入数据 (Ensembl ID Level)
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# 3. 清理 ID 版本号
clean_ids <- sub("\\..*", "", rownames(txi$counts))
rownames(txi$counts) <- clean_ids
rownames(txi$abundance) <- clean_ids
rownames(txi$length) <- clean_ids

# 4. 获取注释 (Protein Coding 筛选用)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = rownames(txi$counts),
                    mart = mart)

# ==============================================================================
# Part 2: 核心步骤 - 聚合到 Gene Symbol 层面
# ==============================================================================
cat("\n[Part 2] Aggregating Ensembl IDs to Gene Symbols...\n")

# 1. 筛选 Protein Coding 并关联 Symbol
# 我们先在 ID 层面筛选出 Protein Coding，避免非编码基因干扰 Symbol 聚合
pc_ids <- gene_annot %>% 
  filter(gene_biotype == "protein_coding") %>%
  filter(!is.na(external_gene_name) & external_gene_name != "")

# 仅保留数据中存在的 ID
common_ids <- intersect(rownames(txi$counts), pc_ids$ensembl_gene_id)
pc_ids <- pc_ids %>% filter(ensembl_gene_id %in% common_ids)

# 提取子集矩阵
mat_counts <- txi$counts[common_ids, ]
mat_tpm    <- txi$abundance[common_ids, ]
mat_len    <- txi$length[common_ids, ]

# 创建映射向量 (ID -> Symbol)
id_to_symbol <- pc_ids$external_gene_name
names(id_to_symbol) <- pc_ids$ensembl_gene_id
# 确保顺序一致
match_idx <- match(rownames(mat_counts), names(id_to_symbol))
group_symbols <- id_to_symbol[match_idx]

# 2. 执行聚合 (Aggregation)
# Counts 和 TPM: 直接求和
agg_counts <- rowsum(mat_counts, group = group_symbols)
agg_tpm    <- rowsum(mat_tpm, group = group_symbols)

# Length: 加权平均 (Weighted Average by TPM)
# 公式: sum(Length * TPM) / sum(TPM)
# 这是一个更科学的合并基因长度的方法，用于 DESeq2 offset
weighted_len_num <- rowsum(mat_len * mat_tpm, group = group_symbols)
# 避免除以0 (如果某基因在所有样本TPM都是0，直接用简单的平均长度)
agg_tpm_safe <- agg_tpm
agg_tpm_safe[agg_tpm_safe == 0] <- 1 
agg_len <- weighted_len_num / agg_tpm_safe

# 对于全0表达的基因，恢复其简单平均长度
zero_expr_genes <- rownames(agg_tpm)[rowSums(agg_tpm) == 0]
if(length(zero_expr_genes) > 0) {
  # 重新计算简单平均
  simple_len <- rowsum(mat_len, group = group_symbols)
  # 这里只是近似处理，因为全0基因在后续会被过滤掉，所以影响不大，但为了程序健壮性加上
  counts_per_symbol <- table(group_symbols)
  agg_len[zero_expr_genes, ] <- simple_len[zero_expr_genes, ] / as.vector(counts_per_symbol[zero_expr_genes])
}

cat(sprintf("   > Aggregation complete. Unique Symbols: %d\n", nrow(agg_counts)))

# 3. 重建 txi 对象 (Symbol Level)
txi_symbol <- list(
  counts = agg_counts,
  abundance = agg_tpm,
  length = agg_len,
  countsFromAbundance = txi$countsFromAbundance
)

# [FIX] Fix length=0 issue for DESeq2
if (any(txi_symbol$length == 0)) {
  txi_symbol$length[txi_symbol$length == 0] <- 1
}

# ==============================================================================
# Part 3: 构建 DESeq2 对象并应用统一过滤
# ==============================================================================
cat("\n[Part 3] Filtering Non-expressed Genes...\n")

# 1. 构建 dds
dds <- DESeqDataSetFromTximport(txi_symbol, colData = colData, design = ~ condition)

# 2. 统一过滤: 剔除所有样本中 Count 和为 10 的基因
# (注意：这是你要求的 "filter those non-expressed genes by dds")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 这个 filtered_genes 列表将用于后续所有步骤
filtered_genes <- rownames(dds)

cat(sprintf("   > Final Genes remaining (Protein Coding + expressed): %d\n", length(filtered_genes)))

# ==============================================================================
# Part 4: 差异分析 (DEA)
# ==============================================================================
cat("\n[Part 4] Running DESeq2...\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "K", "C"))

setwd("D:/Graduation Project/RNASeq/DEG")
# 保存结果 (现在行名直接就是 Symbol，不需要再 join 注释了！)
write.csv(as.data.frame(res), "DEA_Results_Symbol.csv")

# 筛选显著差异基因
sig_genes_df <- as.data.frame(res) %>% filter(padj < 0.001 & abs(log2FoldChange) > 2)
cat(sprintf("   > Significant DEGs: %d\n", nrow(sig_genes_df)))

# ==============================================================================
# Part 5: 可视化 (PCA, Heatmap & Volcano)
# ==============================================================================
cat("\n[Part 5] Visualization (Symbol Level)...\n")

setwd("D:/Graduation Project/RNASeq/PCA")

# --- 0. Update Metadata for Plotting ---
# Create specific labels for plotting
# "Control" for C
# "*LMNA* knockout" for K (The asterisks tell ggtext to italicize LMNA)
colData$plot_condition <- ifelse(colData$condition == "C", "Control", "*LMNA* knockout")
# Ensure the factor level order (Control first)
colData$plot_condition <- factor(colData$plot_condition, levels = c("Control", "*LMNA* knockout"))

# Define a manual color palette to ensure consistency
my_colors <- c("Control" = "#F8766D", "*LMNA* knockout" = "#00BFC4")

# --- 1. PCA Plots (with Italics & Short Names) ---

# Function to draw PCA
draw_custom_pca <- function(mat, title_text, filename) {
  # Calculate PCA
  pca <- prcomp(t(mat), scale. = TRUE)
  pca_df <- as.data.frame(pca$x)
  
  # Calculate variance
  pca_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  # Merge with metadata
  # Ensure row names match
  pca_df <- merge(pca_df, colData, by="row.names")
  
  p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=plot_condition)) +
    geom_point(size=5) +
    # Use short_name (C1, K1) for labels
    geom_text_repel(aes(label=short_name), size=5, box.padding=0.5, show.legend=FALSE) +
    scale_color_manual(values=my_colors, name="Condition") +
    xlab(paste0("PC1: ", pca_var[1], "% variance")) +
    ylab(paste0("PC2: ", pca_var[2], "% variance")) +
    ggtitle(title_text) +
    theme_bw() +
    theme(
      # This enables markdown parsing for the legend (Italics!)
      legend.text = element_markdown(size=12),
      legend.title = element_text(size=13, face="bold"),
      plot.title = element_text(size=16, face="bold", hjust=0.5),
      axis.title = element_text(size=14),
      axis.text = element_text(size=12)
    )
  
  ggsave(filename, p, width=7, height=6)
  return(p)
}

# A. PCA: VST (Symbol Level)
# Generate vst_data
vst_data <- vst(dds, blind = FALSE)

# We need to extract the matrix from the DESeqTransform object
vst_mat <- assay(vst_data)
# Filter for variable genes (optional but recommended for PCA noise reduction)
# Here we use the filtered_genes list from previous steps
vst_mat_subset <- vst_mat[rownames(vst_mat) %in% filtered_genes, ]

cat("   > Drawing VST PCA...\n")
draw_custom_pca(vst_mat_subset, "PCA of Gene Expression (VST)", "Refined_PCA_VST.pdf")

# B. PCA: Log2(TPM+1)
# Use txi_symbol from previous steps
log_tpm_mat <- log2(txi_symbol$abundance[filtered_genes, ] + 1)

cat("   > Drawing LogTPM PCA...\n")
draw_custom_pca(log_tpm_mat, "PCA: Symbol Level (Log2 TPM)", "Refined_PCA_LogTPM.pdf")


# --- 2. Heatmaps (Z-score, Reordered Columns) ---
# Q: Are we using Z-score?
# A: YES. The argument scale="row" in pheatmap performs Z-score standardization across rows.

setwd("D:/Graduation Project/RNASeq/Heatmap")

# Define Column Order: C1, C2, C3, K1, K2, K3
# Find sample names that correspond to these short names
ordered_indices <- order(colData$condition, colData$short_name)
ordered_samples <- rownames(colData)[ordered_indices]

# Prepare Annotation for Heatmap
# Note: pheatmap legends don't support markdown italics easily.
# We will use plain text "LMNA knockout" here.
annot_df <- data.frame(Condition = ifelse(colData$condition == "C", "Control", "LMNA knockout"))
rownames(annot_df) <- rownames(colData)
annot_colors <- list(Condition = c("Control" = "#F8766D", "LMNA knockout" = "#00BFC4"))

# Select Top 50 DEGs
top_genes <- rownames(sig_genes_df)[order(sig_genes_df$padj)][1:min(50, nrow(sig_genes_df))]

draw_custom_heatmap <- function(mat, title_text, filename) {
  # 1. Subset matrix by genes
  mat_sub <- mat[top_genes, ]
  
  # 2. Reorder columns (Samples) explicitly
  mat_sub <- mat_sub[, ordered_samples]
  
  # 3. Update column names to Short Names (C1, K1...) for the plot
  # We need a mapping vector
  new_colnames <- colData[ordered_samples, "short_name"]
  
  # Draw
  pdf(filename, width=8, height=9)
  pheatmap(mat_sub,
           scale = "row",          # THIS IS THE Z-SCORE CALCULATION
           cluster_rows = TRUE,    # Cluster genes
           cluster_cols = FALSE,   # DO NOT cluster samples (Keep manual order)
           annotation_col = annot_df,
           annotation_colors = annot_colors,
           labels_col = new_colnames, # Use short names on x-axis
           show_rownames = TRUE,
           main = title_text,
           color = colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100))
  dev.off()
}

cat("   > Drawing VST Heatmap...\n")
draw_custom_heatmap(vst_mat, "Heatmap (VST)", "Refined_Heatmap_VST.pdf")

cat("   > Drawing LogTPM Heatmap...\n")
draw_custom_heatmap(log_tpm_mat, "Heatmap of Top 50 DEGs (LogTPM)", "Refined_Heatmap_LogTPM.pdf")


setwd("D:/Graduation Project/RNASeq/Volcano")
# --- 3. Volcano Plot ---

cat("   > Drawing Volcano Plot...\n")

# Prepare data
volc_data <- as.data.frame(res)
volc_data$symbol <- rownames(volc_data)
# Add labels for Top 10 significant genes
volc_data$label <- NA
top10_idx <- order(volc_data$padj)[1:10]
volc_data$label[top10_idx] <- volc_data$symbol[top10_idx]

# Define status
volc_data$diff <- "NO"
volc_data$diff[volc_data$log2FoldChange > 2 & volc_data$padj < 0.001] <- "UP"
volc_data$diff[volc_data$log2FoldChange < -2 & volc_data$padj < 0.001] <- "DOWN"

# Plot
p_volc <- ggplot(volc_data, aes(x=log2FoldChange, y=-log10(padj), col=diff, label=label)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("DOWN"="#336699", "NO"="grey", "UP"="#CC3333"),
                     labels=c("Down-regulated", "Not significant", "Up-regulated")) +
  geom_vline(xintercept=c(-2, 2), linetype="dashed", color="black", alpha=0.5) +
  geom_hline(yintercept=-log10(0.001), linetype="dashed", color="black", alpha=0.5) +
  geom_text_repel(max.overlaps=20, box.padding=0.5) +
  labs(title="Volcano Plot: *LMNA* knockout vs Control",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value",
       color="Status") +
  theme_bw() +
  theme(plot.title = element_markdown(hjust = 0.5, face="bold")) # Allow italics in title

ggsave("Refined_Volcano.pdf", p_volc, width=8, height=7)


setwd("D:/Graduation Project/RNASeq/Enrichment")
# ==============================================================================
# Refined Part 6: Enrichment Analysis (Separated UP and DOWN DEGs)
# ==============================================================================
cat("\n[Part 6] Refined Enrichment Analysis (UP vs DOWN separated)...\n")

library(enrichplot)

# 1. Separate significant DEGs into UP and DOWN
# Because sig_genes_df is already filtered for abs(LFC) > 2, we just split by >0 and <0
up_symbols <- rownames(sig_genes_df %>% filter(log2FoldChange > 0))
down_symbols <- rownames(sig_genes_df %>% filter(log2FoldChange < 0))

cat(sprintf("   > Translating %d Up-regulated and %d Down-regulated genes...\n", 
            length(up_symbols), length(down_symbols)))

# 2. Convert to Entrez IDs
up_ids <- bitr(up_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_ids <- bitr(down_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Store in a list to loop through easily
gene_lists <- list()
if(nrow(up_ids) > 0) gene_lists[["UP"]] <- up_ids$ENTREZID
if(nrow(down_ids) > 0) gene_lists[["DOWN"]] <- down_ids$ENTREZID

if(length(gene_lists) > 0) {
  
  # Loop through UP and DOWN lists separately
  for (direction in names(gene_lists)) {
    cat(sprintf("\n   === Processing %s DEGs ===\n", direction))
    current_genes <- gene_lists[[direction]]
    
    # ----------------------------------------------------------------------------
    # A. GO Enrichment 
    # ----------------------------------------------------------------------------
    cat(sprintf("     > Running GO Analysis (p<0.05, q<0.05) for %s...\n", direction))
    
    ego <- enrichGO(gene          = current_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",      
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,       
                    readable      = TRUE)       
    
    if(!is.null(ego) && nrow(ego) > 0) {
      df_go <- as.data.frame(ego) %>% filter(pvalue < 0.05 & qvalue < 0.05)
      
      if(nrow(df_go) > 0) {
        write.csv(df_go, paste0("Enrichment_Table_GO_", direction, ".csv"), row.names = FALSE)
        
        p_go_dot <- dotplot(ego, showCategory=5, split="ONTOLOGY") + 
          facet_grid(ONTOLOGY~., scale="free") +
          ggtitle(paste0("GO Enrichment (", direction, " DEGs)")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        ggsave(paste0("Enrichment_BubblePlot_GO_", direction, ".pdf"), p_go_dot, width = 10, height = 12)
        cat("       [OK] GO results saved.\n")
      } else {
        cat("       [!] No GO terms left after strict filtering.\n")
      }
    } else {
      cat("       [!] No significant GO terms found.\n")
    }
    
    # ----------------------------------------------------------------------------
    # B. KEGG Enrichment
    # ----------------------------------------------------------------------------
    cat(sprintf("     > Running KEGG Analysis (p<0.05, q<0.2) for %s...\n", direction))
    
    ekegg <- enrichKEGG(gene         = current_genes,
                        organism     = 'hsa',    
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)      
    
    if(!is.null(ekegg) && nrow(ekegg) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      df_kegg <- as.data.frame(ekegg) %>% filter(pvalue < 0.05 & qvalue < 0.05)
      
      if(nrow(df_kegg) > 0) {
        write.csv(df_kegg, paste0("Enrichment_Table_KEGG_", direction, ".csv"), row.names = FALSE)
        
        p_kegg_dot <- dotplot(ekegg, showCategory=20) + 
          ggtitle(paste0("KEGG Pathway Enrichment (", direction, " DEGs)")) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        ggsave(paste0("Enrichment_BubblePlot_KEGG_", direction, ".pdf"), p_kegg_dot, width = 10, height = 12)
        cat("       [OK] KEGG results saved.\n")
      } else {
        cat("       [!] No KEGG terms left after strict filtering.\n")
      }
    } else {
      cat("       [!] No significant KEGG pathways found.\n")
    }
  }
} else {
  cat("   [!] No valid Entrez IDs found for enrichment analysis.\n")
}

cat("\n[Part 6] Finished.\n")

# ==============================================================================
# Part 7: 生成最终矩阵 (应用了 dds 过滤标准)
# ==============================================================================
cat("\n[Part 7] Generating Final Clean Matrices...\n")

# 我们使用 txi_symbol (已聚合) 并取 filtered_genes (dds 的非零子集)
# 这确保了导出的矩阵与 PCA/DEA 使用的数据完全一致

final_counts <- txi_symbol$counts[filtered_genes, ]
final_tpm    <- txi_symbol$abundance[filtered_genes, ]
final_logtpm <- log_tpm_mat # 已经在 Part 5 计算并过滤过

setwd("D:/Graduation Project/RNASeq/expression_matrices")
# 保存
write.csv(final_counts, "Final_Clean_Counts_Symbol.csv")
write.csv(final_tpm, "Final_Clean_TPM_Symbol.csv")
write.csv(final_logtpm, "Final_Clean_LogTPM_Symbol.csv")

cat("   > Pipeline Finished Successfully!\n")
