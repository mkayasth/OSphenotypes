library(limma)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(EnhancedVolcano)
library(GSVA)
library(combinat)
library(ComplexHeatmap)
library(grid)
library(reshape2)

##### Loading bulk-seq RNA & metadata file.

Expression <- readRDS("GeneExp_03192025.RDS")
metadata <- read_delim("Osteosarcoma_CellLines.txt")

# Arranging sample ID in the geneExpression table in the same order as metadata: ALT -> TA.
geneExpression <- Expression[, metadata$depmap_id]
saveRDS(geneExpression, "Output/geneExpression.rds")
# write.csv(geneExpression, "rna-seq.csv")


##### Filtering for lowly expressed genes. Only keeping genes with logTPM > 1 in at least 5% of samples.
keep <- rowSums(geneExpression > 1) >= ceiling(0.05 * ncol(geneExpression))
geneExpression <- geneExpression[keep, ]
rm(keep)

# saveRDS(geneExpression, "Output/geneExpression.rds")

##### Doing differential gene expression using Limma for our logTPM data.

# creating design matrix for limma.
metadata$TMMstatus <- as.factor(metadata$TMMstatus)
design <- model.matrix(~ 0 + TMMstatus, data = metadata)
colnames(design) <- levels(metadata$TMMstatus)


# Estimating array weights to account for sample-specific variability in library sizes.
weights <- arrayWeights(geneExpression, design = design)

# Fitting linear model using limma.
fit <- lmFit(geneExpression, design, weights = weights)
fit <- eBayes(fit, trend = TRUE)

# Contrasts: ALT vs TA && TA Vs. ALT.
contrast.matrix <- makeContrasts(
  altVSta = ALT - TA,
  taVSalt = TA - ALT,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# top results for ALT.
alt_results <- topTable(fit2, coef = "altVSta", number = Inf, adjust = "fdr")

# top results for TA.
ta_results <- topTable(fit2, coef = "taVSalt", number = Inf, adjust = "fdr")

# We only have 16 samples: 9 for ALT and 7 for TA. We, as a result, do not have significant adjusted p-value.
# Filtering for p-value < 0.05.
alt_candidates <- alt_results[alt_results$P.Value < 0.05 & alt_results$logFC > 1.5, ] # only taking positive values (genes upregulated in ALT.)
ta_candidates <- ta_results[ta_results$P.Value < 0.05 & ta_results$logFC > 1.5, ] # only taking positive values (genes upregulated in TA.)


# binding candidate genes together.
all_candidates <- rbind(alt_candidates, ta_candidates)

# genes from ALT will have gene status ALT and genes from TA will have gene status TA.
all_candidates <- all_candidates %>%
  mutate("Gene Status" = case_when(
    rownames(all_candidates) %in% rownames(alt_candidates) ~ "ALT",
    rownames(all_candidates) %in% rownames(ta_candidates) ~ "Telomerase"))

# expression table of just the candidates.
dge_gene <- geneExpression[rownames(geneExpression) %in% rownames(all_candidates), ]

##### Making volcano plot for TA vs. ALT.

ta_expression <- ta_results[, c("logFC", "P.Value")]
alt_expression <- alt_results[, c("logFC", "P.Value")]

ta_expression <- ta_expression %>%
  mutate("Gene Status" = case_when(
    ta_results$logFC > 1.5 & ta_results$P.Value < 0.05 ~ "Upregulated in TA samples.",
    ta_results$logFC < -1.5 & ta_results$P.Value < 0.05 ~ "Downregulated in TA samples.",
    TRUE ~"Not Significant"
  ))

# # Top 5 Genes in terms of fold change (to label in the volcano plot). Red denotes most upregulated in TA and blue denotes most downregulated in ALT phenotype.
# ta_top5 <- ta_expression %>%
#   arrange(desc(logFC)) %>%
#   slice_head(n = 5)
# 
# ta_bottom5 <- ta_expression %>%
#   arrange(logFC) %>%
#   slice_head(n = 5)
# 
# top_labels <- bind_rows(
#   ta_top5 %>% rownames_to_column("Gene"),
#   ta_bottom5 %>% rownames_to_column("Gene")
# )


# #Volcano plot ~ all samples.
# ggplot(data = ta_expression, aes(x = logFC, y = -log10(P.Value), color = `Gene Status`)) +
#   geom_point() + 
#   scale_color_manual(values = c("Upregulated in TA samples." = "red", "Downregulated in TA samples." = "blue")) +
#   theme_classic() + geom_text_repel(data = top_labels, 
#                                     aes(label = Gene),
#                                     vjust = 0.2, hjust = 0.2, size = 4,
#                                     color = "black", box.padding = 0.5,
#                                     point.padding = 0.5, max.overlaps = Inf) +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 15, face = "bold"), 
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 13, face = "bold"))
#   
# 
# rm(top_labels)

#Volcano plot ~ all samples.

# top labels -- signature genes obtained from iterative selection performed at the end of this file.
top_labels <- c("GASK1B", "SYCP2", "H2BC11", "MAB21L2", "CLSTN3", "ELOVL4", "MAP9", "COL24A1", "PGM2L1", "STK32A", "CRYAB", "ARMH4", "LFNG", "SYTL5", "CCL28",
                "TERT", "CCDC8", "ALDH1A3", "ERICH1", "TDRP", "SHFL", "DNPH1", "IRX2")

EnhancedVolcano(ta_expression,
                lab = rownames(ta_expression),
                x = 'logFC',
                y = 'P.Value',
                selectLab = top_labels,  # highlighting signature genes.
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P-value'),
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                arrowheads = FALSE,
                labFace = 'bold',
                boxedLabels = TRUE,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                colAlpha = 0.8,
                legendLabels = c('Not Significant','Significant logFC','Significant P-value ','Significant P-value & LogFC'),
                col = c('grey80', 'grey50', 'grey25', 'purple'),
                ylim = c(0, 5),
                caption = "Cutoffs: P < 0.05, |Log2FC| > 1.5"
) + theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.caption = element_text(size = 14))




#####################################################################################


##### Next step: t-test & linear test for candidate genes. Since we dont have anything significant if we look at adjusted p-value, looking at other statistics.

### 1) t-test for candidate genes.

# Initializing an empty data frame for storing t-test results.
all_candidates <- rownames_to_column(all_candidates, var = "Gene")

t_test_results <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(all_candidates)) {
  
  gene_id <- all_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting signficant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all c1 sample expression data placed in c1_group and all c2 sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata$TMMstatus == "ALT", drop = FALSE]
  c2_group <- gene_Expression[, metadata$TMMstatus == "TA", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results <- rbind(t_test_results, data.frame(
    Gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  # removing unrequired intermediates formed while forming the t-test table above.
  rm(gene_id)
  rm(gene_Expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
  rm(p_value_t_test)
}


# Adding ALT or TA gene to the table.
t_test_results <- merge(t_test_results, all_candidates[, c("Gene", "Gene Status")],  by.x = "Gene",
                        by.y = "Gene", all.x = TRUE)


# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig <- t_test_results[round(t_test_results$p_value_t_test, 2) <= 0.01, ]

# The adjusted p-values for t-test would not be significant either. Now, running regression test; any gene significant in either will be considered for future analysis.

### 2) Linear regression test for candidate genes.
regression_results <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene)) {
  
  data <- data.frame(Expression = dge_gene[gene, ], 
                     TMMstatus = metadata$TMMstatus)
  
  # Fit linear model.
  model <- lm(Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results <- rbind(regression_results, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # removing intermediates.
  rm(data)
  rm(summary_model)
  rm(gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# adjusted p-values.
regression_results$Adj.P.Value <- p.adjust(regression_results$p_value, method = "fdr")

# Filtering for r-squared > 0.3 & p-value < 0.01.
regression_results_sig <- regression_results %>%
  filter(round(p_value, 2) <= 0.01 & round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.3)
regression_results_sig <- merge(regression_results_sig, all_candidates[, c("Gene", "Gene Status")],  by.x = "Gene",
                                by.y = "Gene", all.x = TRUE)


# boxplots of only the candidate genes.

# layout for subplots.
num_genes <- nrow(regression_results_sig)

# Saving the plot as a PDF.
pdf("Output/rna-seq-regression_results.pdf", width = 10, height = 50)

# layout and margins
par(mfrow = c(ceiling(num_genes / 3), 3), mar = c(6, 6, 6, 1))

regression_test_candidates <- dge_gene[rownames(dge_gene) %in% regression_results_sig$Gene, ]

for (gene in rownames(regression_test_candidates)) {
  plot_data <- data.frame(Expression = regression_test_candidates[gene, ], 
                          Class = metadata$TMMstatus)
  
  # Fit linear model.
  model <- lm(Expression ~ Class, data = plot_data)
  summary_model <- summary(model)
  
  
  # boxplot.
  boxplot(Expression ~ Class, data = plot_data, 
          col = "lightblue", main = gene, 
          ylab = "Expression", xlab = "Class")
  
  # Adding R-squared value.
  legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                    "\np-value =", round(summary_model$coefficients[2, 4], 4)),
         bty = "n", cex = 0.8)
  
  rm(summary_model)
  rm(plot_data)
  
}

dev.off()

### Now merging candidates significant in both t-test & regression test.
geneSignature_candidates <- merge(t_test_results_sig, regression_results_sig, by = "Gene")
geneSignature_candidates <- geneSignature_candidates[, c("Gene", "Gene Status.x")]
colnames(geneSignature_candidates)[2] <- "gene_status"

saveRDS(geneSignature_candidates, "Output/geneSignatureCandidates.rds")



########## heatmap + gsva annotation for ALT and TA-known genes.


gene_alts <- list("ALT Genes" = c("FANCM", "FANCD2", "BLM", "SMARCA5", "CHD4", "BRCA2", "TOP3A"))
gene_tas <- list("TA Genes" = c("TERT", "DKC1", "GAR1", "MYC"))

gene_combined <- list(
  "ALT Genes" = c("FANCM", "FANCD2", "BLM", "SMARCA5", "CHD4", "BRCA2", "TOP3A"),
  "TA Genes" = c("TERT", "DKC1", "GAR1", "MYC"))


gene_list <- unlist(gene_combined, use.names = FALSE)
gene_categories <- unlist(gene_combined)
row_split_vector <- setNames(
  rep(names(gene_combined), times = lengths(gene_combined)),
  unlist(gene_combined, use.names = FALSE)
)

geneExpression_scaled_ordered <- t(scale(t(geneExpression[rownames(geneExpression) %in% gene_list, ,drop = FALSE])))
geneExpression_scaled_ordered <- geneExpression_scaled_ordered[gene_list, ,drop = FALSE]

# running gsva for alt genes first.
alt_gsva <- gsvaParam(geneExpression, gene_alts, kcdf = "Gaussian")
GSVA_result_alt <- gsva(alt_gsva)
GSVA_df_alt <- as.data.frame(GSVA_result_alt)
GSVA_long_alt <- pivot_longer(GSVA_df_alt, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_alt <- merge(GSVA_long_alt, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_alt <- GSVA_long_alt %>%
  arrange(TMMstatus)
GSVA_long_alt <- GSVA_long_alt[match(colnames(geneExpression_scaled_ordered), GSVA_long_alt$SampleID), , drop = FALSE]

# to match the order as per the data, need to convert SampleID into factor.
GSVA_long_alt$SampleID <- factor(GSVA_long_alt$SampleID, levels = colnames(geneExpression_scaled_ordered))
GSVA_long_alt$Color <- ifelse(GSVA_long_alt$GSVA_Score > 0, "Positive", "Negative")


# running gsva for ta genes.
ta_gsva <- gsvaParam(geneExpression, gene_tas, kcdf = "Gaussian")
GSVA_result_ta <- gsva(ta_gsva)
GSVA_df_ta <- as.data.frame(GSVA_result_ta)
GSVA_long_ta <- pivot_longer(GSVA_df_ta, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_ta <- merge(GSVA_long_ta, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_ta <- GSVA_long_ta %>%
  arrange(TMMstatus)
GSVA_long_ta <- GSVA_long_ta[match(colnames(geneExpression_scaled_ordered), GSVA_long_ta$SampleID), , drop = FALSE]

# to match the order as per the data, need to convert SampleID into factor.
GSVA_long_ta$SampleID <- factor(GSVA_long_ta$SampleID, levels = colnames(geneExpression_scaled_ordered))
GSVA_long_ta$Color <- ifelse(GSVA_long_ta$GSVA_Score > 0, "Positive", "Negative")

# Creating bar plot for GSVA.
combined_annotation <- HeatmapAnnotation(
  "ALT Genes GSVA Score" = anno_barplot(
    GSVA_long_alt$GSVA_Score, 
    gp = gpar(fill = ifelse(GSVA_long_alt$Color == "Positive", "blue", "red"), border =NA,lty="blank"),
    height = unit(1,"cm"),
    direction = "horizontal",
    axis_param = list(at = c(0), labels = c("0")),
  ),

  "TA Genes GSVA Score" = anno_barplot(
    GSVA_long_ta$GSVA_Score, 
    gp = gpar(fill = ifelse(GSVA_long_ta$Color == "Positive", "blue", "red"), border =NA,lty="blank"),
    height = unit(1,"cm"),
    axis_param = list(at = c(0), labels = c("0")),
    direction = "horizontal"
  ),
  gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

# Final heatmap + barPlots.
h1 <- Heatmap(
  geneExpression_scaled_ordered,
  border = TRUE,
  row_labels = rownames(geneExpression_scaled_ordered), 
  column_labels = GSVA_long_ta$cell_line_display_name,#colnames(geneExpression_scaled_ordered),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_split = metadata$TMMstatus,
  row_split = row_split_vector,
  row_title_gp = gpar(fontsize = 14, col = "black", fontface = "bold"),
  row_order = rownames(geneExpression_scaled_ordered),
  column_order = colnames(geneExpression_scaled_ordered),
  height = unit(6, "cm"),
  width = unit(8, "cm"),
  heatmap_legend_param = list(labels_gp=gpar(fontsize = 7),
                              grid_width = unit(0.5, "cm"),
                              #title_gp = gpar(fontsize = 10),
                              legend_height = unit(2, "cm"), legend_width = unit(0.5,"cm"),
                              direction= "vertical"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "black", fill = NA, lwd = 0.8))
  }
)


draw(h1, heatmap_legend_side = "right", annotation_legend_side = "right", 
     merge_legends = TRUE, 
     padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()

########## heatmap + gsva annotation for our signature -- from the best_signature function done at the end of this file.

gene_alts <- list("ALT Genes" = c("GASK1B", "SYCP2", "H2BC11", "MAB21L2", "CLSTN3", 
                                  "ELOVL4", "MAP9", "COL24A1", "PGM2L1", "STK32A", "CRYAB", 
                                  "ARMH4", "LFNG", "SYTL5", "CCL28"))
gene_tas <- list("TA Genes" = c("TERT", "CCDC8", "ALDH1A3", 
                                "ERICH1", "TDRP", "SHFL", 
                                "DNPH1", "IRX2"))

gene_combined <- list(
  "ALT Genes" = c("GASK1B", "SYCP2", "H2BC11", "MAB21L2", "CLSTN3", 
                  "ELOVL4", "MAP9", "COL24A1", "PGM2L1", "STK32A", "CRYAB", 
                  "ARMH4", "LFNG", "SYTL5", "CCL28"),
  "TA Genes" = c("TERT", "CCDC8", "ALDH1A3", 
                 "ERICH1", "TDRP", "SHFL", 
                 "DNPH1", "IRX2"))


gene_list <- unlist(gene_combined, use.names = FALSE)
gene_categories <- unlist(gene_combined)
row_split_vector <- setNames(
  rep(names(gene_combined), times = lengths(gene_combined)),
  unlist(gene_combined, use.names = FALSE)
)

geneExpression_scaled_ordered <- t(scale(t(geneExpression[rownames(geneExpression) %in% gene_list, ,drop = FALSE])))
geneExpression_scaled_ordered <- geneExpression_scaled_ordered[gene_list, ,drop = FALSE]

# running gsva for alt genes first.
alt_gsva <- gsvaParam(geneExpression, gene_alts, kcdf = "Gaussian")
GSVA_result_alt <- gsva(alt_gsva)
GSVA_df_alt <- as.data.frame(GSVA_result_alt)
GSVA_long_alt <- pivot_longer(GSVA_df_alt, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_alt <- merge(GSVA_long_alt, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_alt <- GSVA_long_alt %>%
  arrange(TMMstatus)
GSVA_long_alt <- GSVA_long_alt[match(colnames(geneExpression_scaled_ordered), GSVA_long_alt$SampleID), , drop = FALSE]

# to match the order as per the data, need to convert SampleID into factor.
GSVA_long_alt$SampleID <- factor(GSVA_long_alt$SampleID, levels = colnames(geneExpression_scaled_ordered))
GSVA_long_alt$Color <- ifelse(GSVA_long_alt$GSVA_Score > 0, "Positive", "Negative")


# running gsva for ta genes.
ta_gsva <- gsvaParam(geneExpression, gene_tas, kcdf = "Gaussian")
GSVA_result_ta <- gsva(ta_gsva)
GSVA_df_ta <- as.data.frame(GSVA_result_ta)
GSVA_long_ta <- pivot_longer(GSVA_df_ta, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_ta <- merge(GSVA_long_ta, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_ta <- GSVA_long_ta %>%
  arrange(TMMstatus)
GSVA_long_ta <- GSVA_long_ta[match(colnames(geneExpression_scaled_ordered), GSVA_long_ta$SampleID), , drop = FALSE]

# to match the order as per the data, need to convert SampleID into factor.
GSVA_long_ta$SampleID <- factor(GSVA_long_ta$SampleID, levels = colnames(geneExpression_scaled_ordered))
GSVA_long_ta$Color <- ifelse(GSVA_long_ta$GSVA_Score > 0, "Positive", "Negative")

# Creating bar plot for GSVA.
combined_annotation <- HeatmapAnnotation(
  "ALT Genes GSVA Score" = anno_barplot(
    GSVA_long_alt$GSVA_Score, 
    gp = gpar(fill = ifelse(GSVA_long_alt$Color == "Positive", "blue", "red"), border =NA,lty="blank"),
    height = unit(1,"cm"),
    direction = "horizontal",
    axis_param = list(at = c(0), labels = c("0")),
  ),
  
  "TA Genes GSVA Score" = anno_barplot(
    GSVA_long_ta$GSVA_Score, 
    gp = gpar(fill = ifelse(GSVA_long_ta$Color == "Positive", "blue", "red"), border =NA,lty="blank"),
    height = unit(1,"cm"),
    axis_param = list(at = c(0), labels = c("0")),
    direction = "horizontal"
  ),
  gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

# Final heatmap + barPlots.
h1 <- Heatmap(
  geneExpression_scaled_ordered,
  border = TRUE,
  row_labels = rownames(geneExpression_scaled_ordered), 
  column_labels = GSVA_long_ta$cell_line_display_name,#colnames(geneExpression_scaled_ordered),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_split = metadata$TMMstatus,
  row_split = row_split_vector,
  row_title_gp = gpar(fontsize = 14, col = "black", fontface = "bold"),
  row_order = rownames(geneExpression_scaled_ordered),
  column_order = colnames(geneExpression_scaled_ordered),
  height = unit(12, "cm"),
  width = unit(10, "cm"),
  heatmap_legend_param = list(labels_gp=gpar(fontsize = 7),
                              grid_width = unit(0.5, "cm"),
                              #title_gp = gpar(fontsize = 10),
                              legend_height = unit(2, "cm"), legend_width = unit(0.5,"cm"),
                              direction= "vertical"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "black", fill = NA, lwd = 0.8))
  }
)


draw(h1, heatmap_legend_side = "right", annotation_legend_side = "right", 
     merge_legends = TRUE, 
     padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()




##################################################################################



###################################################################################

# ##### ssGSEA test.
# 
# ### First, we will do ssGSEA test from all ALT genes and all Telomerase genes.
# 
# #  Telomerase set.
#gene_set_list_ta <- list(geneSignature_candidates$Gene[geneSignature_candidates$gene_status == "Telomerase"])
# ta_ssgsea <- ssgseaParam(geneExpression, gene_set_list_ta, normalize = TRUE)
# ssGSEA_result_ta <- gsva(ta_ssgsea)
# ssGSEA_df_ta <- as.data.frame(ssGSEA_result_ta)
# ssGSEA_long_ta <- pivot_longer(ssGSEA_df_ta, cols = everything(), names_to = "depmap_id", values_to = "ssGSEA_Score")
# ssGSEA_long_ta <- merge(ssGSEA_long_ta, metadata, by = "depmap_id", all.x = TRUE)
# 
# 
# # ALT set.
#gene_set_list_alt <- list(geneSignature_candidates$Gene[geneSignature_candidates$gene_status == "ALT"])
# alt_ssgsea <- ssgseaParam(geneExpression, gene_set_list_alt, normalize = TRUE)
# ssGSEA_result_alt <- gsva(alt_ssgsea)
# ssGSEA_df_alt <- as.data.frame(ssGSEA_result_alt)
# ssGSEA_long_alt <- pivot_longer(ssGSEA_df_alt, cols = everything(), names_to = "depmap_id", values_to = "ssGSEA_Score")
# ssGSEA_long_alt <- merge(ssGSEA_long_alt, metadata, by = "depmap_id", all.x = TRUE)
# 
# rm(ssGSEA_result_ta)
# rm(ssGSEA_result_alt)
# 
# ### Now, we will do ssGSEA test for all possible combination of 3 genes. Making combinations using combn function.
# 
# #combn does not take list.
# gene_set_list_ta <- unlist(gene_set_list_ta)
# 
# # 3 gene combination (first for TA genes).
# gene_combinations3_ta <- combn(gene_set_list_ta, 3, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations3_ta, 
#                       sapply(gene_combinations3_ta, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets, normalize = TRUE)
# 
# # ssGSEA scores.
# ssGSEA_results_ta_3Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets) # removing this intermediate.
# rm(ssgsea_params)
# 
# # Combination of 3 genes in ALT.
# gene_set_list_alt <- unlist(gene_set_list_alt)
# 
# gene_combinations3_alt <- combn(gene_set_list_alt, 3, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations3_alt, 
#                       sapply(gene_combinations3_alt, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets)
# 
# # ssGSEA scores.
# ssGSEA_results_alt_3Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets)
# rm(ssgsea_params)
# 
# ### All possible combination of 4 genes.
# 
# # 4 gene combination for TA.
# gene_combinations4_ta <- combn(gene_set_list_ta, 4, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations4_ta, 
#                       sapply(gene_combinations4_ta, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets, normalize = TRUE)
# 
# # ssGSEA scores.
# ssGSEA_results_ta_4Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets) # removing this intermediate.
# rm(ssgsea_params)
# 
# # 4 gene combination for ALT.
# gene_combinations4_alt <- combn(gene_set_list_alt, 4, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations4_alt, 
#                       sapply(gene_combinations4_alt, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets)
# 
# # ssGSEA scores.
# ssGSEA_results_alt_4Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets)
# rm(ssgsea_params)
# 
# ### All possible combination of 5 genes.
# 
# # 5 gene combination for TA.
# gene_combinations5_ta <- combn(gene_set_list_ta, 5, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations5_ta, 
#                       sapply(gene_combinations5_ta, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets, normalize = TRUE)
# 
# # ssGSEA scores.
# ssGSEA_results_ta_5Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets) # removing this intermediate.
# rm(ssgsea_params)
# 
# 
# # 5 gene combination for ALT.
# gene_combinations5_alt <- combn(gene_set_list_alt, 5, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations5_alt, 
#                       sapply(gene_combinations5_alt, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets)
# 
# # ssGSEA scores.
# ssGSEA_results_alt_5Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets)
# rm(ssgsea_params)
# 
# ### All possible combination of 6 genes.
# 
# # 6 gene combination for TA.
# gene_combinations6_ta <- combn(gene_set_list_ta, 6, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations6_ta, 
#                       sapply(gene_combinations6_ta, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets, normalize = TRUE)
# 
# # ssGSEA scores.
# ssGSEA_results_ta_6Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets) # removing this intermediate.
# rm(ssgsea_params)
# 
# # 6 gene combination for ALT.
# gene_combinations6_alt <- combn(gene_set_list_alt, 6, simplify = FALSE)
# 
# gene_sets <- setNames(gene_combinations6_alt, 
#                       sapply(gene_combinations6_alt, paste, collapse = "_"))
# 
# # ssGSEA parameter object.
# ssgsea_params <- ssgseaParam(geneExpression,
#                              gene_sets)
# 
# # ssGSEA scores.
# ssGSEA_results_alt_6Genes <- gsva(ssgsea_params)
# 
# rm(gene_sets)
# rm(ssgsea_params)
# 
# #################################################################################################
# 
# ###### We have a lot of combinations and we need to narrow it down.
# 
# ### finding which gene combination has the highest mean difference of ssgsea score for two conditions. Function to do so. (this is for ALT.)
# top_ssgsea_combination <- function(ssgsea_table, metadata) {
#   alt_cols <- list()
#   ta_cols <- list()
#   
#   # making two tables based on conditions.
#   for (i in colnames(ssgsea_table)) {
#     if (metadata$TMMstatus[metadata$depmap_id == i] == 'ALT') {
#       alt_cols[[i]] <- ssgsea_table[, i] 
#     }
#     
#     else if (metadata$TMMstatus[metadata$depmap_id == i] == 'TA') {
#       ta_cols[[i]] <- ssgsea_table[, i]
#     }
#   }
#   
#   # list into matrices.
#   alt_scores <- do.call(cbind, alt_cols)
#   ta_scores <- do.call(cbind, ta_cols)
#   
#   # subtracting mean from each gene set.
#   alt_means <- rowMeans(alt_scores, na.rm = TRUE)
#   ta_means <- rowMeans(ta_scores, na.rm = TRUE)
#   
#   mean_differences <- alt_means - ta_means
#   
#   # Create result to store data frame: mean difference of ssGSEA for each gene combination.
#   ssGSEA_difference_alt <- data.frame(
#     gene_combination = names(mean_differences),
#     alt_mean = alt_means,
#     ta_mean = ta_means,
#     mean_difference = mean_differences,
#     stringsAsFactors = FALSE
#   )
#   
#   # Sort by absolute difference (descending)
#   ssGSEA_difference_alt <- ssGSEA_difference_alt[order(-ssGSEA_difference_alt$mean_difference), ]
#   
#   return(ssGSEA_difference_alt)
#   
#   
# }
# 
# alt_5Genes <- top_ssgsea_combination(ssGSEA_results_alt_5Genes, metadata)
# alt_4Genes <- top_ssgsea_combination(ssGSEA_results_alt_4Genes, metadata)
# alt_3Genes <- top_ssgsea_combination(ssGSEA_results_alt_3Genes, metadata)
# alt_6Genes <- top_ssgsea_combination(ssGSEA_results_alt_6Genes, metadata)
# 
# 
# ### finding which gene combination has the highest mean difference of ssgsea score for two conditions. Function to do so. (this is for TA.)
# top_ssgsea_combination_ta <- function(ssgsea_table, metadata) {
#   alt_cols <- list()
#   ta_cols <- list()
#   
#   # making two tables based on conditions.
#   for (i in colnames(ssgsea_table)) {
#     if (metadata$TMMstatus[metadata$depmap_id == i] == 'ALT') {
#       alt_cols[[i]] <- ssgsea_table[, i]
#     }
#     
#     
#     else if (metadata$TMMstatus[metadata$depmap_id == i] == 'TA') {
#       ta_cols[[i]] <- ssgsea_table[, i]
#     }
#   }
#   
#   # list into matrices.
#   alt_scores <- do.call(cbind, alt_cols)
#   ta_scores <- do.call(cbind, ta_cols)
#   
#   # subtracting mean from each gene set.
#   alt_means <- rowMeans(alt_scores, na.rm = TRUE)
#   ta_means <- rowMeans(ta_scores, na.rm = TRUE)
#   
#   
#   mean_differences <- ta_means - alt_means
#   
#   
#   # Create result to store data frame: mean difference of ssGSEA for each gene combination.
#   ssGSEA_difference_ta <- data.frame(
#     gene_combination = names(mean_differences),
#     alt_mean = alt_means,
#     ta_mean = ta_means,
#     mean_difference = mean_differences,
#     stringsAsFactors = FALSE
#   )
#   
#   # Sort by absolute difference (descending)
#   ssGSEA_difference_ta <- ssGSEA_difference_ta[order(-ssGSEA_difference_ta$mean_difference), ]
#   
#   return(ssGSEA_difference_ta)
#   
#   
# }
# 
# ta_3Genes <- top_ssgsea_combination_ta(ssGSEA_results_ta_3Genes, metadata)
# ta_4Genes <- top_ssgsea_combination_ta(ssGSEA_results_ta_4Genes, metadata)
# ta_5Genes <- top_ssgsea_combination_ta(ssGSEA_results_ta_5Genes, metadata)
# ta_6Genes <- top_ssgsea_combination_ta(ssGSEA_results_ta_6Genes, metadata)
# 
# 
# ### Now, we will bind all combinations together for a phenotype. Then, calculate p-value from t-test, filter for significant difference across two phenotypes.
# 
# ta_allGenes <- rbind(ta_3Genes, ta_4Genes, ta_5Genes, ta_6Genes)
# ta_allGenes <- ta_allGenes[order(-ta_allGenes$mean_difference), ]
# 
# 
# ssGSEA_results_ta_allGenes <- rbind(ssGSEA_results_ta_3Genes, ssGSEA_results_ta_4Genes, ssGSEA_results_ta_5Genes,ssGSEA_results_ta_6Genes)
# 
# 
# alt_allGenes <- rbind(alt_3Genes, alt_4Genes, alt_5Genes, alt_6Genes)
# alt_allGenes <- alt_allGenes[order(-alt_allGenes$mean_difference), ]
# 
# ssGSEA_results_alt_allGenes <- rbind(ssGSEA_results_alt_3Genes, ssGSEA_results_alt_4Genes, ssGSEA_results_alt_5Genes,ssGSEA_results_alt_6Genes)
# 
# # We have a lot of ALT gene combinations. Picking up top 50000 gene combinations with highest mean difference.
# alt_allGenes <- alt_allGenes[1:50000, ]
# ssGSEA_results_alt_allGenes <- ssGSEA_results_alt_allGenes[rownames(ssGSEA_results_alt_allGenes) %in% rownames(alt_allGenes), ]
# 
# 
# ### t-test for each combination.
# 
# # first for ta group.
# ssgsea_t_test_ta <- data.frame(
#   gene_combination = character(),
#   alt_mean = numeric(),
#   ta_mean = numeric(),
#   mean_difference = numeric(),
#   p_value_t_test = numeric(),
#   fdr_t_test = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # all_p_values <- numeric() # storage for fdr calculation.
# 
# # Looping thru each gene combination:
# for (i in 1:nrow(ta_allGenes)) {
#   gene_combination_ta <- rownames(ta_allGenes)[i]
#   alt_mean <- ta_allGenes[i, "alt_mean"]
#   ta_mean <- ta_allGenes[i, "ta_mean"]
#   mean_difference_ta <- ta_allGenes[i, "mean_difference"]
#   ssgsea_values <- ssGSEA_results_ta_allGenes[rownames(ssGSEA_results_ta_allGenes) == gene_combination_ta, ,drop = FALSE]
#   
#   # alt and ta group.
#   alt_group <- ssgsea_values[, metadata$TMMstatus == "ALT", drop = FALSE]
#   ta_group <- ssgsea_values[, metadata$TMMstatus == "TA", drop = FALSE]
#   
#   
#   # t-test.
#   t_test <- t.test(ta_group, alt_group)
#   p_value_t_test <- t_test$p.value
#   # all_p_values <- c(all_p_values, p_value_t_test)
#   
#   # storing results.
#   ssgsea_t_test_ta <- rbind(ssgsea_t_test_ta, data.frame(
#     gene_combination = gene_combination_ta,
#     alt_mean = alt_mean,
#     ta_mean = ta_mean,
#     mean_difference = mean_difference_ta,
#     p_value_t_test = p_value_t_test
#   ))
#   
#   rm(gene_combination_ta)
#   rm(mean_difference_ta)
#   rm(ssgsea_values)
#   rm(alt_group)
#   rm(ta_group)
#   rm(t_test)
#   rm(p_value_t_test)
# }
# 
# # Filtering with p < 0.01 for t-test.
# # ssgsea_t_test_c2$fdr_t_test <- p.adjust(all_p_values, method = "BH")
# ssgsea_t_test_ta <- ssgsea_t_test_ta[round(ssgsea_t_test_ta$p_value_t_test, 2) <= 0.01,  ] # filter: p-value <  0.01.
# 
# 
# # Now t-test for ALT group.
# ssgsea_t_test_alt <- data.frame(
#   gene_combination = character(),
#   alt_mean = numeric(),
#   ta_mean = numeric(),
#   mean_difference = numeric(),
#   p_value_t_test = numeric(),
#   fdr_t_test = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # all_p_values <- numeric() 
# 
# # Looping thru each gene combination:
# for (i in 1:nrow(alt_allGenes)) {
#   gene_combination_alt <- rownames(alt_allGenes)[i]
#   alt_mean <- alt_allGenes[i, "alt_mean"]
#   ta_mean <- alt_allGenes[i, "ta_mean"]
#   mean_difference_alt <- alt_allGenes[i, "mean_difference"]
#   ssgsea_values <- ssGSEA_results_alt_allGenes[rownames(ssGSEA_results_alt_allGenes) == gene_combination_alt, ,drop = FALSE]
#   
#   
#   # alt and ta group.
#   alt_group <- ssgsea_values[, metadata$TMMstatus == "ALT", drop = FALSE]
#   ta_group <- ssgsea_values[, metadata$TMMstatus == "TA", drop = FALSE]
#   
#   # t-test.
#   t_test <- t.test(ta_group, alt_group)
#   p_value_t_test <- t_test$p.value
#   # all_p_values <- c(all_p_values, p_value_t_test)
#   
#   # storing results.
#   ssgsea_t_test_alt <- rbind(ssgsea_t_test_alt, data.frame(
#     gene_combination = gene_combination_alt,
#     alt_mean = alt_mean,
#     ta_mean = ta_mean,
#     mean_difference = mean_difference_alt,
#     p_value_t_test = p_value_t_test
#   ))
#   
#   rm(gene_combination_alt)
#   rm(mean_difference_alt)
#   rm(ssgsea_values)
#   rm(alt_group)
#   rm(ta_group)
#   rm(t_test)
#   rm(p_value_t_test)
# }
# 
# # ssgsea_t_test_c1$fdr_t_test <- p.adjust(all_p_values, method = "BH")
# ssgsea_t_test_alt <- ssgsea_t_test_alt[round(ssgsea_t_test_alt$p_value_t_test, 2) <= 0.01, ] # filter: p-value < 0.01.
# 
# # all gene combinations are significant. Moving on to regression test and keeping top 10 combinations for each phenotype.
# 
# 
# ##### Linear regression test for candidate genes.
# 
# ### First, ALT combinations.
# 
# ssgsea_regression_alt <- data.frame(
#   GeneCombination = character(),
#   estimate = numeric(),
#   p_value = numeric(),
#   R.squared = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # for each gene in the candidate list, updating regression_results.
# for (gene in rownames(ssGSEA_results_alt_allGenes)) {
#   
#   data <- data.frame(Expression = ssGSEA_results_alt_allGenes[gene, ], 
#                      Class = metadata$TMMstatus)
#   
#   # Fit linear model.
#   model <- lm(Expression ~ Class, data = data)
#   summary_model <- summary(model)
#   
#   
#   # Store regression results.
#   ssgsea_regression_alt <- rbind(ssgsea_regression_alt, data.frame(
#     geneCombination = gene,
#     estimate = summary_model$coefficients[2, 1],
#     p_value = summary_model$coefficients[2, 4],
#     R.Squared = summary_model$r.squared
#   ))
#   
#   # removing intermediates.
#   rm(data)
#   rm(summary_model)
#   rm(gene)
#   rm(estimate)
#   rm(p_value)
#   rm(R.squared)
# }
# 
# # order by r-squared value.
# ssgsea_regression_alt <- ssgsea_regression_alt[round(ssgsea_regression_alt$p_value, 2) <= 0.01, ]
# ssgsea_regression_alt <- merge(ssgsea_regression_alt, alt_allGenes, by.x = "geneCombination",
#                               by.y = "gene_combination", all.x = TRUE)
# ssgsea_regression_alt <- ssgsea_regression_alt[order(-ssgsea_regression_alt$R.Squared), ]
# 
# # making boxplots. ssGSEA vs. Condition for different gene combinations. (Top 10 for ALT).
# 
# ssgsea_alt_top10 <- ssgsea_regression_alt[1:10,]
# 
# 
# # layout for subplots.
# num_genes <- nrow(ssgsea_alt_top10)
# 
# # Saving the plot as a PDF.
# pdf("Output/ssgsea_results_alt.pdf", width = 10, height = 10)
# 
# # layout and margins
# par(mfrow = c(ceiling(num_genes / 3), 3), mar = c(4, 4, 2, 1))
# 
# ssgsea_alt_candidates <- ssGSEA_results_alt_allGenes[rownames(ssGSEA_results_alt_allGenes) %in% ssgsea_alt_top10$geneCombination, ]
# 
# for (geneCombination in rownames(ssgsea_alt_candidates)) {
#   plot_data <- data.frame(ssGSEA = ssGSEA_results_alt_allGenes[geneCombination, ], 
#                           Class = metadata$TMMstatus)
#   
#   # Fit linear model.
#   model <- lm(ssGSEA ~ Class, data = plot_data)
#   summary_model <- summary(model)
#   
#   # boxplot.
#   boxplot(ssGSEA ~ Class, data = plot_data, 
#           col = "lightblue", main = geneCombination, 
#           ylab = "ssGSEA Score", xlab = "Class", cex.main = 0.9)
#   
#   # Adding R-squared value.
#   legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
#                                     "\np-value =", round(summary_model$coefficients[2, 4], 4)),
#          bty = "n", cex = 0.8)
#   
#   rm(summary_model)
#   rm(plot_data)
#   
# }
# 
# dev.off()
# 
# 
# 
# ### Now, TA combinations.
# 
# ssgsea_regression_ta <- data.frame(
#   GeneCombination = character(),
#   estimate = numeric(),
#   p_value = numeric(),
#   R.squared = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # for each gene in the candidate list, updating regression_results.
# for (gene in rownames(ssGSEA_results_ta_allGenes)) {
#   
#   data <- data.frame(Expression = ssGSEA_results_ta_allGenes[gene, ], 
#                      Class = metadata$TMMstatus)
#   
#   # Fit linear model.
#   model <- lm(Expression ~ Class, data = data)
#   summary_model <- summary(model)
#   
#   
#   # Store regression results.
#   ssgsea_regression_ta <- rbind(ssgsea_regression_ta, data.frame(
#     geneCombination = gene,
#     estimate = summary_model$coefficients[2, 1],
#     p_value = summary_model$coefficients[2, 4],
#     R.Squared = summary_model$r.squared
#   ))
#   
#   # removing intermediates.
#   rm(data)
#   rm(summary_model)
#   rm(gene)
#   rm(estimate)
#   rm(p_value)
#   rm(R.squared)
# }
# 
# # order by r-squared value.
# ssgsea_regression_ta <- ssgsea_regression_ta[round(ssgsea_regression_ta$p_value, 2) <= 0.01, ]
# ssgsea_regression_ta <- merge(ssgsea_regression_ta, ta_allGenes, by.x = "geneCombination",
#                               by.y = "gene_combination", all.x = TRUE)
# ssgsea_regression_ta <- ssgsea_regression_ta[order(-ssgsea_regression_ta$R.Squared), ]
# 
# # making boxplots. ssGSEA vs. Condition for different gene combinations. (Top 10 for TA).
# ssgsea_ta_top10 <- ssgsea_regression_ta[1:10,]
# ssgsea_ta_top10
# 
# 
# # layout for subplots.
# num_genes <- nrow(ssgsea_ta_top10)
# 
# # Saving the plot as a PDF.
# pdf("Output/ssgsea_results_ta.pdf", width = 10, height = 10)
# 
# # layout and margins
# par(mfrow = c(ceiling(num_genes / 3), 3), mar = c(4, 4, 2, 1))
# 
# ssgsea_ta_candidates <- ssGSEA_results_ta_allGenes[rownames(ssGSEA_results_ta_allGenes) %in% ssgsea_ta_top10$geneCombination, ]
# 
# for (geneCombination in rownames(ssgsea_ta_candidates)) {
#   plot_data <- data.frame(ssGSEA = ssGSEA_results_ta_allGenes[geneCombination, ], 
#                           Class = metadata$TMMstatus)
#   
#   # Fit linear model.
#   model <- lm(ssGSEA ~ Class, data = plot_data)
#   summary_model <- summary(model)
#   
#   # boxplot.
#   boxplot(ssGSEA ~ Class, data = plot_data, 
#           col = "lightblue", main = geneCombination, 
#           ylab = "ssGSEA Score", xlab = "Class", cex.main = 0.9)
#   
#   # Adding R-squared value.
#   legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
#                                     "\np-value =", round(summary_model$coefficients[2, 4], 4)),
#          bty = "n", cex = 0.8)
#   
#   rm(summary_model)
#   rm(plot_data)
#   
# }
# 
# dev.off()
# 
# 
# 
# #################################################################################
# 
# 
# ##### running GSVA with the gene combinations that gave the highest correlation value.
# # Selecting highest R-squared.
# 
# ### for ALT: highest R-squared: CLSTN3_CRYAB_GASK1B_MAB21L2_MAGEA6_SYTL5.
# # Running gsva on ALT set first.
# gene_set_list_gsva <- list("CLSTN3_CRYAB_GASK1B_MAB21L2_MAGEA6_SYTL5" = c("CLSTN3", "CRYAB", "GASK1B", "MAB21L2", "MAGEA6", "SYTL5"))
# alt_gsva <- gsvaParam(geneExpression, gene_set_list_gsva, kcdf = "Gaussian")
# GSVA_result_alt <- gsva(alt_gsva)
# GSVA_df_alt <- as.data.frame(GSVA_result_alt)
# GSVA_long_alt <- pivot_longer(GSVA_df_alt, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
# GSVA_long_alt <- merge(GSVA_long_alt, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
# GSVA_long_alt <- GSVA_long_alt %>%
#   arrange(TMMstatus)
# 
# 
# ### for TA: highest R-squared: SHFL_TBX1_TDRP_TERT
# gene_set_list_gsva <- list("SHFL_TBX1_TDRP_TERT" = c("SHFL", "TBX1", "TDRP", "TERT"))
# ta_gsva <- gsvaParam(geneExpression, gene_set_list_gsva, kcdf = "Gaussian")
# GSVA_result_ta <- gsva(ta_gsva)
# GSVA_df_ta <- as.data.frame(GSVA_result_ta)
# GSVA_long_ta <- pivot_longer(GSVA_df_ta, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
# GSVA_long_ta <- merge(GSVA_long_ta, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
# GSVA_long_ta <- GSVA_long_ta %>%
#   arrange(TMMstatus)
# 
# 
# ### bar graphs for GSVA: ALT and TA.
# 
# # First ALT.
# GSVA_long_alt$SampleID <- factor(GSVA_long_alt$SampleID, levels = unique(GSVA_long_alt$SampleID))
# ggplot(GSVA_long_alt, aes(x = SampleID, y = GSVA_Score, fill = TMMstatus)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   labs(title = "GSVA Scores by Sample for ALT genes CLSTN3_CRYAB_GASK1B_MAB21L2_MAGEA6_SYTL5",
#        x = "Sample ID",
#        y = "GSVA Score",
#        color = "Class")
# 
# # Now, Telomerase samples.
# GSVA_long_ta$SampleID <- factor(GSVA_long_ta$SampleID, levels = unique(GSVA_long_ta$SampleID))
# ggplot(GSVA_long_ta, aes(x = SampleID, y = GSVA_Score, fill = TMMstatus)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   labs(title = "GSVA Scores by Sample for TA genes SHFL_TBX1_TDRP_TERT",
#        x = "Sample ID",
#        y = "GSVA Score",
#        color = "Class")
# 
# #########################################################################################
# 
# 

# ### Instead of boxplots, making a line graph. Literature OS genes.
# genes = unique(c("TP53","RB1",
#                  "ATRX",
#                  "MDM2",
#                  "TERT","MYCN", "DLG2", "PTPRQ", "ZFHX4",
#                  "TP53", "ALK", "FLG", "KMT2D", "EWSR1", "RB1", "WT1", "CDKN2A", "KMT2B", "NF1",
#                  "MAP3K4", "ATRX", "CIC", "ASPCR1", "STAG2", "IGF1R", "PTEN", "CDK4", "PDGFRA2",
#                  "MYC", "BRAF", "KRAS", "NRAS", "IGF1R", "KMT2C", "PALB2", "CCNE1", "COPS3", "DLEU1", "KDR"))
# 
# num_genes = length(genes)
# mutationVSexpression <- geneExpression[rownames(geneExpression) %in% genes, ]
# 
# expr_long <- melt(mutationVSexpression)
# colnames(expr_long) <- c("Gene", "Sample", "Expression")
# 
# 
# # Sort samples for consistent plotting
# expr_long$Sample <- factor(expr_long$Sample, levels = unique(expr_long$Sample))
# 
# # Create and save to PDF
# pdf("LineGraphs_OnePerGene.pdf", width = 10, height = 5)
# # Loop through each gene and plot separately
# for (gene in genes) {
#   p <- ggplot(subset(expr_long, Gene == gene), 
#               aes(x = Sample, y = Expression, group = 1)) +
#     geom_line(color = "steelblue", size = 1) +
#     geom_point(color = "black", size = 2) +
#     theme_minimal() + theme_classic2() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),
#           plot.title = element_text(hjust = 0.5)) +
#     labs(title = paste("Expression of", gene),
#          x = "Sample", y = "Expression Level")
#   
#   print(p)
# }
# 
# dev.off()
# 
# # layout for subplots.
# num_genes <- length(genes)
# 
# # Saving the plot as a PDF.
# pdf("mutationVScopyNumber4.pdf", width = 10, height = 20)
# 
# # layout and margins
# par(mfrow = c(ceiling(num_genes / 3), 3), mar = c(4, 4, 2, 1))
# 
# mutationVSexpression <- geneExpression[rownames(geneExpression) %in% genes, ]
# 
# for (gene in rownames(mutationVSexpression)) {
#   plot_data <- data.frame(geneInterest = mutationVSexpression[gene, ], 
#                           Class = metadata$Class)
#   
#   # Fit linear model.
#   model <- lm(geneInterest ~ Class, data = plot_data)
#   summary_model <- summary(model)
#   
#   # boxplot.
#   boxplot(geneInterest ~ Class, data = plot_data, 
#           col = "lightblue", main = gene, 
#           ylab = "Expression", xlab = "Class", cex.main = 0.9)
#   
#   # Adding R-squared value.
#   legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
#                                     "\np-value =", round(summary_model$coefficients[2, 4], 4)),
#          bty = "n", cex = 0.8)
#   
#   rm(summary_model)
#   rm(plot_data)
#   
# }
# 
# dev.off()
# 
# 
# genes = c("MUC19", "MUC4", "C17orf97", "KRTAP9-6", "LILRB3", "MT-CYB", "MT-ND5", "MUC16", "MUC21", "MUC3A", "MYOM3", "TTN", "ADAM29", "AHNAK2", "BTAF1", "CCDC88B", "EP400", "KMT2D", "MUC6", "XIRP2", "HLA-DRB1", "CACNA1B", "DSPP", 
#           "FCGBP", "MDC1", "RB1", "ABCA1", "ARHGAP17", "ARHGEF18", "FAM8A1", "OGFR", "USH2A")
# 
# num_genes = length(genes)
# mutationVSexpression <- geneExpression[rownames(geneExpression) %in% genes, ]
# library(reshape2)
# expr_long <- melt(mutationVSexpression)
# colnames(expr_long) <- c("Gene", "Sample", "Expression")
# 
# 
# # Sort samples for consistent plotting
# expr_long$Sample <- factor(expr_long$Sample, levels = unique(expr_long$Sample))
# 
# # Create and save to PDF
# pdf("LineGraphs_OnePerTop20.pdf", width = 10, height = 5)
# # Loop through each gene and plot separately
# for (gene in genes) {
#   p <- ggplot(subset(expr_long, Gene == gene), 
#               aes(x = Sample, y = Expression, group = 1)) +
#     geom_line(color = "steelblue", size = 1) +
#     geom_point(color = "black", size = 2) +
#     theme_minimal() + theme_classic2() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),
#           plot.title = element_text(hjust = 0.5)) +
#     labs(title = paste("Expression of", gene),
#          x = "Sample", y = "Expression Level")
#   
#   print(p)
# }
# 
# dev.off()



#################################################################################

# ssgsea function: finding the best signature by exhaustive method.

gene_set_list_ta <- list(geneSignature_candidates$Gene[geneSignature_candidates$gene_status == "Telomerase"])
gene_set_list_alt <- list(geneSignature_candidates$Gene[geneSignature_candidates$gene_status == "ALT"])


# making the function first.
best_signature <- function(pivot_gene, candidate_genes, expression_matrix, phenotype_vector, 
                                    phenotype_labels = c("ALT", "TA")) {
  
  # initializing starting conditions.
  current_signature <- pivot_gene  #starting with just 1 gene.
  best_pval <- Inf # highest value - gene combination with lower p-value replaces this every loop.
  all_results <- data.frame(Genes = I(list()), p_value = numeric(0))
  used_genes <- current_signature
  
  repeat {
    test_results <- list() # list to store list of signature gene tried and its p-value for each round.
    
    # Genes not used yet for current iteration.
    remaining_genes <- setdiff(candidate_genes, used_genes)
    if (length(remaining_genes) == 0) 
      break # all combinations have been tried..
    
    for (gene in remaining_genes) {
      gene <- as.character(gene)
      test_signature <- c(current_signature, gene)
      gene_set <- list(sig = test_signature) #gsva package takes a list; so making a list.
      
      # running ssGSEA.
      ssgsea_params <- ssgseaParam(expression_matrix,
                                   gene_set, normalize = TRUE)
      
      # ssGSEA scores.
      ssgsea_scores <- gsva(ssgsea_params)
      
      scores_alt <- ssgsea_scores[, phenotype_vector == phenotype_labels[1]]
      scores_ta  <- ssgsea_scores[, phenotype_vector == phenotype_labels[2]]
      
      pval <- t.test(scores_alt, scores_ta)$p.value
      test_results[[gene]] <- list(sig = test_signature, pval = pval)
    }
    
    # Tracking all combinations tested in the current iteration.
    for (combinations in test_results) {
      all_results <- rbind(all_results, data.frame(Genes = I(list(combinations$sig)), p_value = combinations$pval))
    }
    
    # Checking if any combination has an improved p-value from the t-test.
    min_pval <- min(sapply(test_results, `[[`, "pval"))
    best_gene <- names(which.min(sapply(test_results, `[[`, "pval")))
    best_this_round <- test_results[[best_gene]]
    
    if (min_pval < best_pval) {
      current_signature <- best_this_round$sig
      best_pval <- min_pval
      used_genes <- unique(c(used_genes, best_gene))
    } 
    
    else {
      break  # No improvement in the current round.
    }
  }
  
  # finally, sorting all combinations by increasing p-value.
  all_results <- all_results[order(all_results$p_value), ]
  best_results <- all_results %>%
    mutate(GeneList = sapply(Genes, paste, collapse = ", ")) %>% ## from the list of all results, picking the best 50.
    dplyr::select(GeneList, p_value) %>%
    arrange(p_value) %>%
    head(50) 
  
  return(best_results)
  
}  

gene_set_list_ta <- unlist(gene_set_list_ta) ## we make the list inside the function. Not unlisting gave an error.
gene_set_list_alt <- unlist(gene_set_list_alt)

# use of function on TA phenotype.
result_ta <- best_signature(
  pivot_gene = "TERT",
  candidate_genes = gene_set_list_ta,
  expression_matrix = geneExpression,
  phenotype_vector = metadata$TMMstatus,
  phenotype_labels = c("ALT", "TA")
)

## looking at gsva scores for this combination.
gene_set_list_gsva <- list("TERT_CCDC8_ALDH1A3_ERICH1_TDRP_SHFL_DNPH1_IRX2" = c("TERT", "CCDC8", "ALDH1A3", "ERICH1", "TDRP", "SHFL", "DNPH1", "IRX2"))

ta_gsva <- gsvaParam(geneExpression, gene_set_list_gsva, kcdf = "Gaussian")
GSVA_result_ta <- gsva(ta_gsva)
GSVA_df_ta <- as.data.frame(GSVA_result_ta)
GSVA_long_ta <- pivot_longer(GSVA_df_ta, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_ta <- merge(GSVA_long_ta, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_ta <- GSVA_long_ta %>%
  arrange(TMMstatus)

# barplot.
GSVA_long_ta$SampleID <- factor(GSVA_long_ta$SampleID, levels = unique(GSVA_long_ta$SampleID))
ggplot(GSVA_long_ta, aes(x = SampleID, y = GSVA_Score, fill = TMMstatus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.title = element_text(size = 14)) + 
  labs(x = "Sample ID",
       y = "GSVA Score",
       color = "Class")

# linear test.
model <- lm(GSVA_Score ~ TMMstatus, data = GSVA_long_ta)
summary_model <- summary(model)

# violin plot.
r_text <- paste0(
  "R² = ", round(summary_model$adj.r.squared, 3))

ggplot(GSVA_long_ta, aes(x = TMMstatus, y = GSVA_Score, color = TMMstatus, fill = TMMstatus)) +
  geom_violin(size = 0.2, draw_quantiles=c(0.5),alpha=0.5) +
  geom_point(position = position_jitter(width = .08), size = 3)  +
  scale_fill_manual(values=c("ALT"="lightblue","TA"="lightpink2")) + 
  scale_color_manual(values=c("ALT"="blue", "TA"="darkred"))+ 
  theme_classic() + 
  labs(x = "Class", y = "GSVA Score") +
  theme(
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 15, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"), 
    legend.position = "none") +
  stat_compare_means(comparisons = list(c("ALT","TA")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 8, tip.length = 0.01) +  
  annotate("text", x = 2, y = 0, label = r_text, size = 8, fontface = "bold") 

# use of function of ALT phenotype.
result_alt <- best_signature(
  pivot_gene = "GASK1B",
  candidate_genes = gene_set_list_alt,
  expression_matrix = geneExpression,
  phenotype_vector = metadata$TMMstatus,
  phenotype_labels = c("ALT", "TA")
)

## looking at gsva scores for this combination.
gene_set_list_gsva <- list("GASK1B_SYCP2_H2BC11_MAB21L2_CLSTN3_ELOVL4_MAP9_COL24A1_PGM2L1_STK32A_CRYAB_ARMH4_LFNG_SYTL5_CCL28" = c("GASK1B", "SYCP2", "H2BC11", "MAB21L2", "CLSTN3", "ELOVL4", "MAP9", "COL24A1", "PGM2L1", "STK32A", "CRYAB", "ARMH4", "LFNG", "SYTL5", "CCL28"))

alt_gsva <- gsvaParam(geneExpression, gene_set_list_gsva, kcdf = "Gaussian")
GSVA_result_alt <- gsva(alt_gsva)
GSVA_df_alt <- as.data.frame(GSVA_result_alt)
GSVA_long_alt <- pivot_longer(GSVA_df_alt, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
GSVA_long_alt <- merge(GSVA_long_alt, metadata, by.x = "SampleID", by.y = "depmap_id", all.x = TRUE)
GSVA_long_alt <- GSVA_long_alt %>%
  arrange(TMMstatus)

# linear test.
model <- lm(GSVA_Score ~ TMMstatus, data = GSVA_long_alt)
summary_model <- summary(model)

# violin plot.
r_text <- paste0(
  "R² = ", round(summary_model$adj.r.squared, 3))

ggplot(GSVA_long_alt, aes(x = TMMstatus, y = GSVA_Score, color = TMMstatus, fill = TMMstatus)) +
  geom_violin(size = 0.2, draw_quantiles=c(0.5),alpha=0.5) +
  geom_point(position = position_jitter(width = .08), size = 3)  +
  scale_fill_manual(values=c("ALT"="lightblue","TA"="lightpink2")) + 
  scale_color_manual(values=c("ALT"="blue", "TA"="darkred"))+ 
  theme_classic() + 
  labs(x = "Class", y = "GSVA Score") +
  theme(
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 15, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"), 
    legend.position = "none") +
  stat_compare_means(comparisons = list(c("ALT","TA")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 8, tip.length = 0.01) +  
  annotate("text", x = 2, y = 0, label = r_text, size = 8, fontface = "bold") 


# barplot.
GSVA_long_alt$SampleID <- factor(GSVA_long_alt$SampleID, levels = unique(GSVA_long_alt$SampleID))
ggplot(GSVA_long_alt, aes(x = SampleID, y = GSVA_Score, fill = TMMstatus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.title = element_text(size = 14)) + 
  labs(x = "Sample ID",
       y = "GSVA Score",
       color = "Class")






# ggplot(data = ta_expression, aes(x = logFC, y = -log10(P.Value), color = `Gene Status`)) +
#   geom_point() + 
#   scale_color_manual(values = c("Upregulated in TA samples." = "red", "Downregulated in TA samples." = "blue")) +
#   theme_classic() + geom_text_repel(data = top_labels, 
#                                     aes(label = Gene),
#                                     vjust = 0.2, hjust = 0.2, size = 4,
#                                     color = "black", box.padding = 0.5,
#                                     point.padding = 0.5, max.overlaps = Inf) +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 15, face = "bold"), 
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 13, face = "bold"))
# 
# 
# rm(top_labels)


##########################################
