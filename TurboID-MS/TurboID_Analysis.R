# ==============================================================================
# TURBOID-MS COMPREHENSIVE ANALYSIS PIPELINE
# ==============================================================================

required_CRAN <- c("readr", "dplyr", "tidyr", "ggplot2", "ggVennDiagram", "factoextra", "FactoMineR")
required_Bioc <- c("limma", "clusterProfiler", "org.Hs.eg.db", "STRINGdb", "enrichplot", "pathview")

# Install missing packages
for(pkg in required_CRAN) { if(!require(pkg, character.only = TRUE)) install.packages(pkg, ask = FALSE) }
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for(pkg in required_Bioc) { if(!require(pkg, character.only = TRUE)) BiocManager::install(pkg, ask = FALSE) }

library(dplyr)
library(ggplot2)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(ggVennDiagram)
library(enrichplot)

# Go to working directory
setwd("D:/Graduation Project/TurboID-MS")

dir.create("Output_Plots", showWarnings = FALSE)
dir.create("Output_Tables", showWarnings = FALSE)

# ==============================================================================
# <<< SCIENTIFIC COLOR SCHEME >>>
# Define a professional, color-blind friendly palette for conditions
# ==============================================================================
cond_colors <- c(
  "D" = "#ED0000", # Red (Emerin_d95-99)
  "Q" = "#42B540", # Green (Emerin_Q133H)
  "P" = "#0099B4", # Teal (Emerin_P183T)
  "S" = "#925E9F", # Purple (Emerin_S54F)
  "W" = "#00468B", # Dark Blue (Wildtype Emerin)
  "N" = "#B0B0B0"  # Grey/Silver (3xNLS Control)
)

# ==============================================================================
# 1. Read Data & Extract Identifiers
# ==============================================================================
message("Reading data...")
file_name <- "20260423_084550_20260422_ZHW_Report_No_Normalized.tsv"
df_raw <- read.delim(file_name, stringsAsFactors = FALSE, check.names = FALSE)

# <<< ID PARSING FOR DISPLAY AND ENRICHMENT >>>
# Extract the first ID from semicolon-separated lists
df_raw$Primary_Accession <- sapply(strsplit(df_raw$PG.ProteinAccessions, ";"), `[`, 1)
raw_names <- sapply(strsplit(df_raw$PG.ProteinNames, ";"), `[`, 1)
df_raw$Primary_Name <- gsub("_.*$", "", raw_names) # e.g., "NUD4B_HUMAN" -> "NUD4B"

# Create a mapping dictionary for downstream functional analysis
id_map <- data.frame(Name = df_raw$Primary_Name, Accession = df_raw$Primary_Accession)

quant_cols <- grep("raw\\.PG\\.Quantity", colnames(df_raw), value = TRUE)
quant_matrix <- df_raw[, quant_cols]

clean_names <- gsub(".*\\[\\d+\\] ZHW_20260422_([A-Z0-9]+)\\.raw\\.PG\\.Quantity", "\\1", quant_cols)
colnames(quant_matrix) <- clean_names

# Initial Distribution
all_numeric <- as.numeric(as.matrix(quant_matrix))
log2_numeric <- log2(all_numeric[!is.na(all_numeric) & all_numeric > 0])
p_dist <- ggplot(data.frame(Value = log2_numeric), aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.5) + theme_minimal() +
  labs(title = "Overall Distribution of Log2(Raw Quantities)", x = "Log2(Quantity)", y = "Density")
ggsave("Output_Plots/1_Raw_Distribution.pdf", p_dist, width = 6, height = 4)

# ==============================================================================
# 2. Impute N/A with 1 & 3. Within-Group Normalization
# ==============================================================================
message("Imputing NAs with 1 and Normalizing...")
quant_matrix[is.na(quant_matrix)] <- 1

groups <- list(D=c("D1","D2","D3"), Q=c("Q1","Q2","Q3"), P=c("P1","P2","P3"), 
               S=c("S1","S2","S3"), W=c("W1","W2","W3"), N=c("N1","N2","N3"))

norm_matrix <- quant_matrix
for (grp in names(groups)) {
  cols <- groups[[grp]]
  medians <- apply(quant_matrix[, cols], 2, function(x) {
    val <- median(x[x > 1], na.rm = TRUE)
    if(is.na(val)) return(1) else return(val)
  })
  max_median <- max(medians)
  for (col in cols) {
    if (medians[col] > 1) norm_matrix[, col] <- norm_matrix[, col] * (max_median / medians[col])
  }
}

# ==============================================================================
# 4. Log2 Transformation & Pre-filtering Cutoff
# ==============================================================================
log_matrix <- log2(norm_matrix)

# <<< SET ROWNAMES TO PROTEIN NAMES FOR DISPLAY >>>
rownames(log_matrix) <- df_raw$Primary_Name

cutoff_value <- quantile(log2_numeric, 0.05, na.rm=TRUE)

# ==============================================================================
# 5. PCA Analysis
# ==============================================================================
message("Performing PCA...")
pca_res <- prcomp(t(log_matrix))
pca_data <- data.frame(Sample = rownames(pca_res$x),
                       PC1 = pca_res$x[,1], PC2 = pca_res$x[,2],
                       Condition = substr(rownames(pca_res$x), 1, 1))

# <<< APPLY PROFESSIONAL COLOR CODE TO PCA >>>
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, show.legend = FALSE) +
  scale_color_manual(values = cond_colors) +  # APPLIED PALETTE HERE
  theme_bw() +
  labs(title = "PCA of Intra-Condition Normalized TurboID Data (Log2)")
ggsave("Output_Plots/2_PCA_Plot.pdf", p_pca, width = 7, height = 6)

# ==============================================================================
# Helper Function for Limma, Volcano Plot, and GO
# ==============================================================================
run_comparison <- function(test_cond, ref_cond, prefix_name, do_GO = FALSE) {
  cols_test <- groups[[test_cond]]
  cols_ref  <- groups[[ref_cond]]
  all_cols  <- c(cols_test, cols_ref)
  
  keep <- apply(log_matrix[, all_cols], 1, max) >= cutoff_value
  sub_matrix <- log_matrix[keep, all_cols]
  
  condition <- factor(rep(c("Test", "Ref"), each = 3), levels = c("Ref", "Test"))
  design <- model.matrix(~ condition)
  fit <- lmFit(sub_matrix, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$ProteinName <- rownames(res)
  res$Status <- "Not_Sig"
  res$Status[res$logFC > 1 & res$adj.P.Val < 0.05] <- "Up"
  res$Status[res$logFC < -1 & res$adj.P.Val < 0.05] <- "Down"
  
  write.csv(res, paste0("Output_Tables/", prefix_name, "_Diff.csv"), row.names = FALSE)
  
  # <<< MATCH VOLCANO COLORS WITH PCA COLORS >>>
  # Up proteins take the Test color, Down proteins take the Ref color
  volc_palette <- c(
    "Up" = unname(cond_colors[test_cond]),
    "Down" = unname(cond_colors[ref_cond]),
    "Not_Sig" = "#DFDFDF" # Subtle Grey
  )
  
  p_volc <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = volc_palette) +
    theme_minimal() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = paste(prefix_name, "Volcano Plot"),
         x="Log2 Fold Change",
         y="-Log10 Adjusted P-value")
  ggsave(paste0("Output_Plots/", prefix_name, "_Volcano.pdf"), p_volc, width = 8, height = 7)
  
  # <<< USE CLUSTERPROFILER TO MAP UNIPROT TO ENTREZ >>>
  if(do_GO) {
    up_names <- res$ProteinName[res$Status == "Up"]
    up_uniprot <- id_map$Accession[id_map$Name %in% up_names]
    
    if(length(up_uniprot) > 5) {
      # Translate UniProt Accessions to Entrez IDs suitable for standard databases
      mapped_ids <- suppressMessages(bitr(up_uniprot, fromType = "UNIPROT", 
                                          toType = "ENTREZID", OrgDb = org.Hs.eg.db))
      
      if(nrow(mapped_ids) > 0) {
        # ont = "ALL" queries BP, CC, and MF simultaneously
        ego <- enrichGO(gene          = mapped_ids$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = "ENTREZID",
                        ont           = "ALL", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE) 
        
        if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          # Use dotplot faceted by ONTOLOGY to clearly separate BP, CC, and MF
          p_go <- dotplot(ego, split = "ONTOLOGY", showCategory = 5) + 
            facet_grid(ONTOLOGY ~ ., space = "free_y", scale = "free_y") +
            labs(title = paste(prefix_name, "GO Enrichment (BP, CC, MF)")) +
            theme(plot.title = element_text(size = 12, face = "bold"))
          
          ggsave(paste0("Output_Plots/", prefix_name, "_GO_All.pdf"), p_go, width = 10, height = 12)
          write.csv(as.data.frame(ego), paste0("Output_Tables/", prefix_name, "_GO_All.csv"))
        }
      }
    }
  }
  return(res)
}

# ==============================================================================
# 6 & 7. Execute Comparisons
# ==============================================================================
message("Running comparisons vs 3xNLS...")
res_WN <- run_comparison("W", "N", "1_Emerin_vs_3xNLS", do_GO = TRUE)
res_DN <- run_comparison("D", "N", "1_Emerin_d95-99_vs_3xNLS", do_GO = TRUE)
res_QN <- run_comparison("Q", "N", "1_Emerin_Q133H_vs_3xNLS", do_GO = TRUE)
res_PN <- run_comparison("P", "N", "1_Emerin_P183T_vs_3xNLS", do_GO = TRUE)
res_SN <- run_comparison("S", "N", "1_Emerin_S54F_vs_3xNLS", do_GO = TRUE)

message("Running comparisons vs WT Emerin...")
res_DW <- run_comparison("D", "W", "2_Emerin_d95-99_vs_Emerin", do_GO = FALSE)
res_QW <- run_comparison("Q", "W", "2_Emerin_Q133H_vs_Emerin", do_GO = FALSE)
res_PW <- run_comparison("P", "W", "2_Emerin_P183T_vs_Emerin", do_GO = FALSE)
res_SW <- run_comparison("S", "W", "2_Emerin_S54F_vs_Emerin", do_GO = FALSE)

# ==============================================================================
# 8 & 9. Venn Diagrams & Enrichment (GO, KEGG, PPI) of Common Genes
# ==============================================================================
message("Generating Venn Diagrams and Common Enrichment...")
list_Up <- list("D"=res_DW$ProteinName[res_DW$Status=="Up"], 
                "Q"=res_QW$ProteinName[res_QW$Status=="Up"],
                "P"=res_PW$ProteinName[res_PW$Status=="Up"], 
                "S"=res_SW$ProteinName[res_SW$Status=="Up"])

list_Dn <- list("D"=res_DW$ProteinName[res_DW$Status=="Down"], 
                "Q"=res_QW$ProteinName[res_QW$Status=="Down"],
                "P"=res_PW$ProteinName[res_PW$Status=="Down"], 
                "S"=res_SW$ProteinName[res_SW$Status=="Down"])

# <<< MODIFICATION: FUNCTION TO INJECT NAMES INTO VENN DIAGRAM >>>
modify_venn_labels <- function(p) {
  # Loop through ggplot layers to find the text/label layer
  for (i in seq_along(p$layers)) {
    if (!is.null(p$layers[[i]]$mapping$label)) {
      layer_data <- p$layers[[i]]$data
      # Ensure this layer contains region items (the actual protein names)
      if (!is.null(layer_data) && "item" %in% names(layer_data) && is.list(layer_data$item)) {
        
        # Create a new custom label: If < 11 items, stack their names. Otherwise, show count.
        layer_data$custom_label <- sapply(layer_data$item, function(x) {
          n <- length(x)
          if (n == 0) return("")
          if (n < 21) return(paste(x, collapse = "\n"))
          return(as.character(n))
        })
        
        # Update the layer data and point the label aesthetic to our new custom_label
        p$layers[[i]]$data <- layer_data
        p$layers[[i]]$mapping$label <- ggplot2::aes(label = custom_label)$label
        
        # Shrink the text slightly and tighten line height so lists fit nicely inside the bubbles
        p$layers[[i]]$aes_params$size <- 2.5
        p$layers[[i]]$aes_params$lineheight <- 0.8
      }
    }
  }
  return(p)
}

# Generate and intercept the UP Venn Diagram
p_venn_up <- ggVennDiagram(list_Up, label = "count") + 
  labs(title="Upregulated vs WT Emerin") + 
  scale_fill_gradient(low="white", high="red")
p_venn_up <- modify_venn_labels(p_venn_up)
ggsave("Output_Plots/3_Venn_Upregulated.pdf", p_venn_up, width=8, height=8) # Increased size to 8x8

# Generate and intercept the DOWN Venn Diagram
p_venn_dn <- ggVennDiagram(list_Dn, label = "count") + 
  labs(title="Downregulated vs WT Emerin") + 
  scale_fill_gradient(low="white", high="blue")
p_venn_dn <- modify_venn_labels(p_venn_dn)
ggsave("Output_Plots/3_Venn_Downregulated.pdf", p_venn_dn, width=8, height=8) # Increased size to 8x8


common_up <- Reduce(intersect, list_Up)
common_dn <- Reduce(intersect, list_Dn)
write.csv(data.frame(ProteinName = common_up), "Output_Tables/3_Common_Up_vs_Emerin.csv", row.names=F)
write.csv(data.frame(ProteinName = common_dn), "Output_Tables/3_Common_Down_vs_Emerin.csv", row.names=F)

# Initialize STRINGdb mapping object
message("Initializing STRINGdb for PPI Networks (this may take a moment to load the database)...")
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=150, input_directory="")

run_common_enrich <- function(prot_names, label) {
  if (length(prot_names) < 2) return(NULL)
  
  uniprots <- id_map$Accession[id_map$Name %in% prot_names]
  
  # --- 1. GO & KEGG ENRICHMENT ---
  mapped <- suppressMessages(bitr(uniprots, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  
  if (nrow(mapped) > 0) {
    ego <- enrichGO(gene = mapped$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      p_go <- dotplot(ego, split = "ONTOLOGY", showCategory = 5) + 
        facet_grid(ONTOLOGY ~ ., space = "free_y", scales = "free_y") +
        labs(title = paste(label, "GO Enrichment (BP, CC, MF)"))
      ggsave(paste0("Output_Plots/4_", label, "_GO_All.pdf"), p_go, width = 10, height = 12)
      write.csv(as.data.frame(ego), paste0("Output_Tables/4_", label, "_GO_All.csv"))
    }
    
    ekegg <- enrichKEGG(gene = mapped$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      p_kegg <- barplot(ekegg, showCategory = 15, title = paste(label, "KEGG Pathways"))
      ggsave(paste0("Output_Plots/4_", label, "_KEGG.pdf"), p_kegg, width = 10, height = 12)
      write.csv(as.data.frame(ekegg), paste0("Output_Tables/4_", label, "_KEGG.csv"))
    }
  }
  
  # --- 2. PROTEIN-PROTEIN INTERACTION (PPI) NETWORK ---
  df_to_map <- data.frame(UNIPROT = uniprots)
  mapped_string <- suppressWarnings(string_db$map(df_to_map, "UNIPROT", removeUnmappedRows = TRUE))
  
  if(nrow(mapped_string) > 1) {
    pdf(paste0("Output_Plots/5_", label, "_PPI_Network.pdf"), width = 8, height = 8)
    string_db$plot_network(mapped_string$STRING_id)
    dev.off()
    
    interactions <- string_db$get_interactions(mapped_string$STRING_id)
    write.csv(interactions, paste0("Output_Tables/5_", label, "_PPI_Interactions.csv"), row.names = FALSE)
  }
}

message("Running PPI and Enrichment for Common Upregulated...")
run_common_enrich(common_up, "Common_Upregulated")

message("Running PPI and Enrichment for Common Downregulated...")
run_common_enrich(common_dn, "Common_Downregulated")

message("Analysis complete! All files saved in Output_Plots and Output_Tables directories.")
