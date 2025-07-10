library(limma)
library(tidyverse)
library(ggrepel)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(maftools)


### loading necessary rds files and metadata.
geneExpression <- readRDS("Output/geneExpression.rds") # from rnaSeq.R
copyNumber <- readRDS("OmicsCNGene.RDS")
metadata <- read_delim("Osteosarcoma_CellLines.txt")

metadataModel <- read_csv("Model.csv") # this particular metadata contains information like ageCategory, sex and sampleClassification (primary or metastatis).
metadataModel <- metadataModel[metadataModel$ModelID %in% metadata$depmap_id, ,drop = FALSE]
metadataModel <- metadataModel[match(metadata$depmap_id, metadataModel$ModelID), ]

copyNumber <- copyNumber[metadata$depmap_id, ]


copyNumber <- copyNumber[match(metadata$depmap_id, rownames(copyNumber)), ]


# arranging the data for copyNumber as per the metadata.
copyNumber <- t(copyNumber)

# 
# ### using edgeR for differential number of copy numbers across each groups.
# 
# # Building model matrix.
# 
# # different factors of the samples.
# class1 <- as.factor(metadata$TMMstatus)
# 
# # model matrix ~ without an intercept term.
# design <- model.matrix(~class1+0)
# 
# # differential expression object.
# copyNumber <- na.omit(copyNumber)
# copyNumber_DGE <- DGEList(counts=copyNumber, group=class1)
# 
# 
# # TMM normalization for calculating normalization factor.
# copyNumber_DGE <- calcNormFactors(copyNumber_DGE)
# 
# 
# # Calculating dispersion and fitting the model.
# d <- estimateDisp(copyNumber_DGE, design, verbose=TRUE)
# fit <- glmQLFit(d, design)
# 
# # contrast parameter (TA-ALT).
# contrast <- makeContrasts(class1TA-class1ALT, levels=design)
# 
# # differential expression test.
# fit2 <- glmQLFTest(fit, contrast = contrast)
# 
# # Adjusted p-value (False discovery rate correction.)
# fit2$table$fdr <- p.adjust(fit2$table$PValue, method ="BH") 
# 
# # Now, in fit2$table, gene name is a rowname. Making it a column of its own called gene.
# fit_expression <- fit2$table
# fit_expression$gene <- rownames(fit_expression)
# 
# # filtering for candidate genes. C2 refers to genes upregulated in TA and C1 refers to ALT phenotype.
# candidate_genes_C2 <- fit_expression[c(fit_expression$PValue <= 0.05 & fit_expression$logFC > 0.5), ]
# # candidate_genes_C1
# 
# candidate_genes_C1 <- fit_expression[c(fit_expression$PValue <= 0.05 & fit_expression$logFC < -0.5), ]
# # candidate_genes_C2
# 
# # binding candidate genes together.
# all_candidates <- rbind(candidate_genes_C1, candidate_genes_C2)
# all_candidates$gene <- rownames(all_candidates)
# 
# mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=113)
# 
# # Retrieving the Ensembl gene IDs and gene biotype (protein coding)
# protein_coding_genes <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
#                               filters="biotype",
#                               values = "protein_coding",
#                               mart = mart)
# 
# # Filtering genes to only include protein-coding genes.
# all_candidates2 <- all_candidates %>%
#   filter(gene %in% protein_coding_genes$hgnc_symbol)
# 
# # genes from C1 will have gene status ALT and genes from C2 will have gene status TA.
# all_candidates2 <- all_candidates2 %>%
#   mutate(`Gene Status` = case_when(
#     gene %in% candidate_genes_C1$gene ~ "Low Copy Number in TA Samples",
#     gene %in% candidate_genes_C2$gene ~ "High Copy Number in TA Samples",
#     TRUE ~ "Not Significant"))
# 
# 
# # doing this for the fit_expression table just for the volcano plot, where all genes are included.
# fit_expression <- fit_expression %>%
#   mutate(`Gene Status` = case_when(
#     gene %in% candidate_genes_C1$gene ~ "Low Copy Number in TA Samples",
#     gene %in% candidate_genes_C2$gene ~ "High Copy Number in TA Samples",
#     TRUE ~ "Not Significant"))
# 
# # Top 10 Genes in terms of fold change (to label in the volcano plot).
# ALT_top5 <- all_candidates2 %>%
#   arrange(desc(logFC)) %>%
#   slice_head(n = 5)
# 
# TA_top5 <- all_candidates2 %>%
#   arrange(logFC) %>%
#   slice_head(n= 5)
# 
# 
# #Volcano plot ~ all samples.
# ggplot(data = fit_expression, aes(x = logFC, y = -log10(PValue), color = `Gene Status`)) +
#   geom_point() + 
#   scale_color_manual(values = c("High Copy Number in TA Samples" = "red", "Low Copy Number in TA Samples" = "blue")) +
#   theme_classic() + geom_text_repel(data = bind_rows(ALT_top5, TA_top5),
#                                   aes(label = gene),
#                                    vjust = 0.5, hjust = 0.5, size = 3,
#                                   color = "black", box.padding = 0.5,
#                                   point.padding = 0.5, max.overlaps = Inf)
# 
# 
# copyNumberCandidates <- copyNumber[rownames(copyNumber) %in% all_candidates2$gene, ]

### using limma for differential number of copy numbers across each groups: data is continuous.

# creating design matrix for limma.
metadata$TMMstatus <- as.factor(metadata$TMMstatus)
design <- model.matrix(~ 0 + TMMstatus, data = metadata)
colnames(design) <- levels(metadata$TMMstatus)


# Estimating array weights to account for sample-specific variability in library sizes.
weights <- arrayWeights(copyNumber, design = design)

# Fitting linear model using limma.
fit <- lmFit(copyNumber, design, weights = weights)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

# Contrasts: ALT vs TA && TA Vs. ALT.
contrast.matrix <- makeContrasts(
  altVSta = ALT - TA,
  taVSalt = TA - ALT,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# top results for ALT.
alt_results <- topTable(fit2, coef = "altVSta", number = Inf, adjust = "fdr")

# top results for TA.
ta_results <- topTable(fit2, coef = "taVSalt", number = Inf, adjust = "fdr")


# Filtering for p-value < 0.05.
alt_candidates <- alt_results[alt_results$P.Value < 0.05 & alt_results$logFC > 0.75, ] # only taking positive values (genes upregulated in ALT.)
ta_candidates <- ta_results[ta_results$P.Value < 0.05 & ta_results$logFC > 0.75, ] # only taking positive values (genes upregulated in TA.)


# binding candidate genes together.
all_candidates <- rbind(alt_candidates, ta_candidates)
all_candidates$gene <- rownames(all_candidates)
# 
# mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=113)
# 
# # Retrieving the Ensembl gene IDs and gene biotype (protein coding)
# protein_coding_genes <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
#                               filters="biotype",
#                               values = "protein_coding",
#                               mart = mart)
# 
# # Filtering genes to only include protein-coding genes.
# all_candidates2 <- all_candidates %>%
#   filter(gene %in% protein_coding_genes$hgnc_symbol)

# genes from ALT will have gene status ALT and genes from TA will have gene status TA.
all_candidates <- all_candidates %>%
  mutate("Gene Status" = case_when(
    rownames(all_candidates) %in% rownames(alt_candidates) ~ "ALT",
    rownames(all_candidates) %in% rownames(ta_candidates) ~ "Telomerase"))

# copy number table of just the candidates. Log fold change of 0.75 and p-value less than 0.05 considered.
copyNumberCandidates <- copyNumber[rownames(copyNumber) %in% all_candidates$gene, ]

# ##### Making volcano plot for TA vs. ALT.

ta_expression <- ta_results[, c("logFC", "P.Value")]
alt_expression <- alt_results[, c("logFC", "P.Value")]
# 
# ta_expression <- ta_expression %>%
#   mutate("Gene Status" = case_when(
#     ta_results$logFC > 0.75 & ta_results$P.Value < 0.05 ~ "High Copy Number in TA samples.",
#     ta_results$logFC < -0.75 & ta_results$P.Value < 0.05 ~ "Low Copy Number in TA samples.",
#     TRUE ~"Not Significant"
#   ))

# Top 5 Genes in terms of fold change (to label in the volcano plot). Red denotes most copy number in TA and blue denotes least copy number in ALT phenotype.
ta_top5 <- ta_expression %>%
  filter(P.Value < 0.05) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5)

ta_bottom5 <- ta_expression %>%
  filter(P.Value < 0.05) %>%
  arrange(logFC) %>%
  slice_head(n = 5)

top_labels <- bind_rows(
  ta_top5 %>% rownames_to_column("Gene"),
  ta_bottom5 %>% rownames_to_column("Gene")
)

# .
# ggplot(data = ta_expression, aes(x = logFC, y = -log10(P.Value), color = `Gene Status`)) +
#   geom_point() + 
#   scale_color_manual(values = c("High Copy Number in TA samples." = "red", "Low Copy Number in TA samples." = "blue")) +
#   theme_classic() + geom_text_repel(data = top_labels,
#                                     aes(label = Gene),
#                                     vjust = 0.2, hjust = 0.2, size = 4.0,
#                                     color = "black", box.padding = 0.3,
#                                     point.padding = 0.3, max.overlaps = Inf) +
#   theme(
#     legend.position = c(0.85, 0.5),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     #legend.box.margin = margin(0, 2, 0, 2),
#     #legend.spacing = unit(1, "mm"),
#     axis.title.x = element_text(size = 18),
#     axis.title.y = element_text(size = 18),
#     axis.text.x  = element_text(size = 12, face = "bold"),
#     axis.text.y  = element_text(size = 12, face = "bold")
#   )

# Volcano plot ~ all samples
EnhancedVolcano(ta_expression,
                lab = rownames(ta_expression),
                x = 'logFC',
                y = 'P.Value',
                selectLab = top_labels$Gene,  # highlighting signature genes.
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P-value'),
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.05,
                FCcutoff = 0.75,
                pointSize = 2.0,
                arrowheads = FALSE,
                labFace = 'bold',
                boxedLabels = TRUE,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colAlpha = 0.8,
                ylim = c(0, 7),
                legendLabels = c('Not Significant','Significant logFC','Significant P-value ','Significant P-value & LogFC'),
                col = c('grey80', 'grey50', 'grey25', 'purple'),
                caption = "Cutoffs: P < 0.05, |Log2FC| > 0.75"
) + theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.caption = element_text(size = 14))

rm(top_labels)


### Next step: t-test & linear test for candidate genes. Since we dont have anything significant if we look at adjusted p-value, looking at other statistics.

### 1) t-test for candidate genes.

# Initializing an empty data frame for storing t-test results.
t_test_results <- data.frame(
  gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(all_candidates)) {
  
  gene_id <- all_candidates$gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting signficant genes (all samples for this gene).
  copy_Number <- copyNumberCandidates[rownames(copyNumberCandidates) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all c1 sample expression data placed in c1_group and all c2 sample expression data placed in c2_group.
  c1_group <- copy_Number[, metadata$TMMstatus == "ALT", drop = FALSE]
  c2_group <- copy_Number[, metadata$TMMstatus == "TA", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results <- rbind(t_test_results, data.frame(
    gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  # removing unrequired intermediateries formed while forming the t-test table above.
  rm(gene_id)
  rm(gene_Expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
  rm(p_value_t_test)
}


# Adding c1 or c2 gene to the table.
t_test_results <- merge(t_test_results, all_candidates[, c("gene", "Gene Status")], by = "gene", all.x = TRUE)

# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig <- t_test_results[t_test_results$p_value_t_test < 0.01, ]


### 2) Linear regression test for candidate genes.
regression_results <- data.frame(
  gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)


# layout for subplots.
num_genes <- nrow(copyNumberCandidates)

# Saving the plot as a PDF.
pdf("Output/copyNumberBoxplots.pdf", width = 20, height = 600)

# layout and margins
par(mfrow = c(ceiling(num_genes / 3), 3), mar = c(3, 3, 2, 1))

# for each gene in the candidate list, update regression_results and boxplot.
for (gene in rownames(copyNumberCandidates)) {
  
  plot_data <- data.frame(Expression = as.numeric(copyNumberCandidates[gene, ]), 
                          Class = metadata$TMMstatus)
  
  # Fit linear model.
  model <- lm(Expression ~ Class, data = plot_data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results <- rbind(regression_results, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # boxplot.
  boxplot(Expression ~ Class, data = plot_data, 
          col = "lightblue", main = gene, 
          ylab = "Expression", xlab = "Class")
  
  # Adding R-squared value.
  legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                    "\np-value =", round(summary_model$coefficients[2, 4], 3)),
         bty = "n", cex = 0.8)
  
  rm(plot_data)
  rm(summary_model)
  rm(Gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# Close the PDF device
dev.off()

# adjusted p-values.
regression_results$Adj.P.Value <- p.adjust(regression_results$p_value, method = "fdr")

regression_results_sig <- regression_results %>%
  filter(p_value < 0.01 & R.Squared > 0.3)

regression_results_sig <- merge(regression_results_sig, all_candidates[, c("gene", "Gene Status")], by.x = "Gene", by.y = "gene", all.x = TRUE)

# only using genes passing both regression and t-test. In our data, all regression test significant genes are significant as per t-test.
regression_results_sig <- regression_results_sig[regression_results_sig$Gene %in% t_test_results_sig$gene, ,drop = FALSE]

#####################################################################


##### gene expression of these differential copy number.
expression_matrix <- geneExpression[rownames(geneExpression) %in% regression_results_sig$Gene, ,drop = FALSE] # out of 22, only 9 are present in expression dataset consisting of protein coding genes.


### add t-test for expression difference of these expression matrix.

# Initializing an empty data frame for storing t-test results.
t_test_results_expression <- data.frame(
  gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(expression_matrix)) {
  
  gene_id <- rownames(expression_matrix)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting signficant genes (all samples for this gene).
  expression <- expression_matrix[rownames(expression_matrix) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all c1 sample expression data placed in c1_group and all c2 sample expression data placed in c2_group.
  c1_group <- expression[, metadata$TMMstatus == "ALT", drop = FALSE]
  c2_group <- expression[, metadata$TMMstatus == "TA", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results_expression <- rbind(t_test_results_expression, 
                                      data.frame(gene = gene_id, p_value_t_test = p_value_t_test))
  
  # removing unrequired intermediateries formed while forming the t-test table above.
  rm(gene_id)
  rm(expression)
  rm(c1_group)
  rm(c2_group)
  rm(t_test)
  rm(p_value_t_test)
}

# this is for heatmap.
expression_matrix <- scale(expression_matrix)

copy_number_matrix <- copyNumber[rownames(copyNumber) %in% rownames(expression_matrix), ,drop = FALSE]
shortlisted_genes <- rownames(copy_number_matrix)


# Output PDF
pdf("Output/gene_expression_cnv_heatmaps.pdf", width = 5, height = 4)



for (gene in shortlisted_genes) {
  # Subset expression and copy number for the gene
  expr_vec <- as.numeric(expression_matrix[gene, ])
  cn_vec <- as.numeric(copy_number_matrix[gene, ])
  
  # Turn into matrix for heatmap (expression)
  expr_mat <- matrix(expr_vec, nrow = 1)
  rownames(expr_mat) <- gene
  colnames(expr_mat) <- colnames(expression_matrix)
  
  # Get p-values from t-test tables
  expr_pval <- t_test_results_expression$p_value_t_test[t_test_results_expression$gene == gene]
  cnv_pval <- t_test_results$p_value_t_test[t_test_results$gene == gene]
  
  
  
  # Define top barplot annotation for CN
  bar_anno <- columnAnnotation(
    CopyNumber = anno_barplot(cn_vec, 
                              border = FALSE,
                              gp = gpar(fill = "darkblue")),
    annotation_name_side = "left"
  )
  
  
  
  
  # Create the heatmap
  ht <- Heatmap(expr_mat,
                name = "Expression",
                top_annotation = bar_anno,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                column_split = metadata$TMMstatus,
                column_title_gp = gpar(fontsize = 5),
                # column_title = paste0(
                #   "RNA-Seq Expression p-value: ", signif(expr_pval, 3), "\n",
                #   "CN difference p-value: ", signif(cnv_pval, 3), "\n"
                # ),
                column_names_gp = gpar(fontsize = 4),
                col = colorRamp2(c(min(expr_vec), mean(expr_vec), max(expr_vec)), c("blue", "white", "red")),
                height = unit(1, "cm"))
  
  # Draw heatmap
  draw(ht)
  # Adding custom title above the plot.
  grid.text(
    label = paste0(
      "RNA-Seq Expression p-value: ", signif(expr_pval, 3), "\n",
      "CN difference p-value: ", signif(cnv_pval, 3)
    ),
    x = unit(0.45, "npc"), y = unit(0.80, "npc"), just = "top",
    gp = gpar(fontsize = 5))
  
}

dev.off()

# Looking at t-test: ERICH1 has significant difference between ALT and TA in terms of both expression and copy number.

###########################################


##### looking at copy number + expression of genes from OS + cancer literature.
OS_gene <- c("TP53", "RB1", "MDM2", "ATRX", "TERT", "MYCN", "DLG2", "PTPRQ", "ZFHX4", "ALK",
             "FLG", "KMT2D", "EWSR1", "WT1", "CDKN2A", "KMT2B", "NF1",
             "STAG2", "IGF1R", "PTEN", "CDK4", "MYC", "KRAS", "PALB2", "CCNE1",
             "COPS3", "DLEU1", "KDR")


expression_matrix <- geneExpression[rownames(geneExpression) %in% OS_gene, ,drop = FALSE]
copy_number_matrix <- copyNumber[rownames(copyNumber) %in% rownames(expression_matrix), ,drop = FALSE]

# t-test for copy number of differential expressed genes.
t_test_results_cn <- data.frame(
  gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(copy_number_matrix)) {
  
  gene_id <- rownames(copy_number_matrix)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting signficant genes (all samples for this gene).
  copy_Number <- copy_number_matrix[rownames(copy_number_matrix) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all c1 sample expression data placed in c1_group and all c2 sample expression data placed in c2_group.
  c1_group <- copy_Number[, metadata$TMMstatus == "ALT", drop = FALSE]
  c2_group <- copy_Number[, metadata$TMMstatus == "TA", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results_cn <- rbind(t_test_results_cn, data.frame(
    gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  # removing unrequired intermediateries formed while forming the t-test table above.
  rm(gene_id)
  rm(copy_Number)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
}

# Adding t-test p-values for expression_matrix.
# Initializing an empty data frame for storing t-test results.
t_test_results_expression <- data.frame(
  gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(expression_matrix)) {
  
  gene_id <- rownames(expression_matrix)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting signficant genes (all samples for this gene).
  expression <- expression_matrix[rownames(expression_matrix) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all c1 sample expression data placed in c1_group and all c2 sample expression data placed in c2_group.
  c1_group <- expression[, metadata$TMMstatus == "ALT", drop = FALSE]
  c2_group <- expression[, metadata$TMMstatus == "TA", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results_expression <- rbind(t_test_results_expression, data.frame(
    gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  # removing unrequired intermediateries formed while forming the t-test table above.
  rm(gene_id)
  rm(expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
}

expression_matrix <- scale(expression_matrix)




# Output PDF
pdf("Output/gene_expression_cnv_heatmaps2.pdf", width = 5, height = 3)



for (gene in rownames(copy_number_matrix)) {
  # Subset expression and copy number for the gene
  expr_vec <- as.numeric(expression_matrix[gene, ])
  cn_vec <- as.numeric(copy_number_matrix[gene, ])
  
  # Turn into matrix for heatmap (expression)
  expr_mat <- matrix(expr_vec, nrow = 1)
  rownames(expr_mat) <- gene
  colnames(expr_mat) <- colnames(expression_matrix)
  
  # Get p-values from t-test tables
  expr_pval <- t_test_results_expression$p_value_t_test[t_test_results_expression$gene == gene]
  cnv_pval <- t_test_results_cn$p_value_t_test[t_test_results_cn$gene == gene]
  
  # Define top barplot annotation for CN
  bar_anno <- columnAnnotation(
    CopyNumber = anno_barplot(cn_vec, 
                              border = FALSE,
                              gp = gpar(fill = "darkblue")),
    annotation_name_side = "left"
  )
  
  # Create the heatmap
  ht <- Heatmap(expr_mat,
                name = "Expression",
                top_annotation = bar_anno,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                column_split = metadata$TMMstatus,
                column_title_gp = gpar(fontsize = 5),
                # #column_title = paste0(
                #   "RNA-Seq Expression p-value: ", signif(expr_pval, 3), "\n",
                #   "CN difference p-value: ", signif(cnv_pval, 3), "\n"
                # ),
                column_names_gp = gpar(fontsize = 4),
                col = colorRamp2(c(min(expr_vec), mean(expr_vec), max(expr_vec)), c("blue", "white", "red")),
                height = unit(1, "cm"))
  
  # Draw heatmap
  draw(ht)
  
  # Add custom title above the plot
  grid.text(
    label = paste0(
      "RNA-Seq Expression p-value: ", signif(expr_pval, 3), "\n",
      "CN difference p-value: ", signif(cnv_pval, 3)
    ),
    x = unit(0.45, "npc"), y = unit(0.80, "npc"), just = "top",
    gp = gpar(fontsize = 5))
}

dev.off()

# CDKN2A (overexpression), PTEN (underexpression), MDM2(overexpression), COPS3(overexpression), CDK4 (overexpression), MYC (overexpression) matches with changes in copy number. 

##################################################################

## similar to what we just did, but adding heatmaps for everything: expression, fusion, mutation, copy number.

##### looking at copy number + expression of genes + fusion + mutation from OS literature with heatmap.
OS_gene <- c("TP53", "RB1", "MDM2", "MYCN", "DLG2", "PTPRQ", "ZFHX4", "ALK",
             "FLG", "KMT2D", "EWSR1", "WT1", "CDKN2A", "KMT2B", "NF1",
             "STAG2", "IGF1R", "PTEN", "CDK4", "MYC", "KRAS", "PALB2", "CCNE1",
             "COPS3", "DLEU1", "KDR")

### expression and copy number data.
expression_matrix <- t(scale(t(geneExpression[rownames(geneExpression) %in% OS_gene, ,drop = FALSE])))

# replacing sampleID with cell line ID.
colnames(expression_matrix) <- metadata$cell_line_display_name[
  match(colnames(expression_matrix), metadata$depmap_id)]

copy_number_matrix <- t(scale(t(copyNumber[rownames(copyNumber) %in% rownames(expression_matrix), ,drop = FALSE])))
# replacing sampleID with cell line ID.
colnames(copy_number_matrix) <- metadata$cell_line_display_name[
  match(colnames(copy_number_matrix), metadata$depmap_id)]


expression_matrix <- expression_matrix[rownames(expression_matrix) %in% rownames(copy_number_matrix), ,drop = FALSE]


### gene fusion data.
geneFusion <- readRDS("Output/geneFusion.rds") # from geneFusion.R.

# Converting geneFusion to appropriate format.
geneFusion <- geneFusion %>%
  mutate(Present = 1) %>%
  distinct(Gene, SampleID, .keep_all = TRUE) %>%
  pivot_wider(names_from = SampleID, values_from = Present, values_fill = 0)

geneFusion <- geneFusion %>%
  dplyr::select(2:last_col())

geneFusion <- geneFusion[geneFusion$Gene %in% rownames(expression_matrix), ,drop = FALSE]
missing_genes <- setdiff(rownames(expression_matrix), geneFusion$Gene)

# creating a 0-filled dataframe for the missing genes: not doing this gave error.
zero_rows <- matrix(0, nrow = length(missing_genes), ncol = ncol(geneFusion) - 1) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(geneFusion)[-1]  # match sample column names
zero_rows <- cbind(Gene = missing_genes, zero_rows)

geneFusion <- bind_rows(geneFusion, zero_rows)
geneFusion <- geneFusion %>%
  group_by(Gene) %>%
  summarise(across(everything(), max), .groups = "drop")

geneFusion <- as.data.frame(geneFusion)
rownames(geneFusion) <- geneFusion$Gene
geneFusion$Gene = NULL

# replacing sampleID with cell line ID.
colnames(geneFusion) <- metadata$cell_line_display_name[
  match(colnames(geneFusion), metadata$depmap_id)]



### mutation data.
somaticMutation <- readRDS("Output/somaticMutation.rds")
somaticMutation <- somaticMutation@data
somaticMutation <- somaticMutation %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)

somaticMutation <- somaticMutation %>%
  mutate(Variant_Classification = as.character(Variant_Classification)) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode,
              values_from = Variant_Classification,
              values_fill = "None")

somaticMutation <- as.data.frame(somaticMutation)
rownames(somaticMutation) <- somaticMutation$Hugo_Symbol
somaticMutation$Hugo_Symbol = NULL

somaticMutation <- somaticMutation[rownames(somaticMutation) %in% rownames(copy_number_matrix), ,drop = FALSE]

# creating a "None"-filled dataframe for the missing genes

missing_genes <- setdiff(rownames(expression_matrix), rownames(somaticMutation))
zero_rows <- matrix("None", nrow = length(missing_genes), ncol = ncol(somaticMutation)) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(somaticMutation)
rownames(zero_rows) <- missing_genes

zero_rows <- cbind(zero_rows)

somaticMutation <- bind_rows(somaticMutation, zero_rows)


### matching the rows and columns before making the heatmap.
copy_number_matrix <- copy_number_matrix[rownames(expression_matrix), colnames(expression_matrix)]
somaticMutation <- somaticMutation[rownames(expression_matrix), colnames(expression_matrix)]
geneFusion <- geneFusion[rownames(expression_matrix), colnames(expression_matrix)]


## for adding Male or Female category.
sexMetadata <- as.data.frame(metadataModel[, c("StrippedCellLineName", "Sex"), drop = FALSE])  #subsetting required columns.
rownames(sexMetadata) <- sexMetadata$StrippedCellLineName
sexMetadata$StrippedCellLineName = NULL
sexMetadata = t(sexMetadata)

## adding categories: primary or metastasis.
sampleClassification <- as.data.frame(metadataModel[, c("StrippedCellLineName", "PrimaryOrMetastasis"), drop = FALSE])  #subsetting required columns.
rownames(sampleClassification) <- sampleClassification$StrippedCellLineName
sampleClassification$StrippedCellLineName = NULL
sampleClassification = t(sampleClassification)
rownames(sampleClassification)[rownames(sampleClassification) == "PrimaryOrMetastasis"] <- "Sample Classification" # renaming the rowname.

# ageCategory <- as.data.frame(metadataModel[, c("ModelID", "AgeCategory"), drop = FALSE])  #subsetting required columns.
# rownames(ageCategory) <- ageCategory$ModelID
# ageCategory$ModelID = NULL
# ageCategory = t(ageCategory)
# rownames(ageCategory)[rownames(ageCategory) == "AgeCategory"] <- "Age Category" # renaming the rowname.

fusion_colors <- c("0" = "gray90", "1" = "black")

mutation_colors <- c(
  "Missense_Mutation" = "green",
  "Nonsense_Mutation" = "red",
  "Frame_Shift_Del" = "blue",
  "Multi_Hit" = "black",
  "Splice_Site" = "orange",
  "None" = "gray90"  # no mutation
)

sexColors <- c(
  "Male" = "black",
  "Female" = "gray90"
)

sampleClassification_colors <- c(
  "Primary" = "darkgreen",
  "Metastatic" = "red"
)




ht1 <- Heatmap(expression_matrix,
                      name = "Expression",
                      cluster_columns = FALSE,
                      column_title = "Expression",
                      cluster_rows = FALSE,
                      show_column_names = TRUE,
                      row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                      column_split = metadata$TMMstatus,
                      column_title_gp = gpar(fontsize = 0),
                      column_names_gp = gpar(fontsize = 6),
                      height = unit(6, "cm"),
                      width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht2 <- Heatmap(copy_number_matrix,
               name = "Copy Number",
               cluster_columns = FALSE,
               column_title = "Copy Number",
               cluster_rows = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht3 <- Heatmap(geneFusion, name = "Fusion", 
               col = fusion_colors, 
               cluster_rows = FALSE, 
               column_title = "Fusion",
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })
              
ht4 <- Heatmap(somaticMutation, name = "Mutations", 
               col = mutation_colors, 
               cluster_rows = FALSE, 
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })
          

ht5 <- Heatmap(sexMetadata, name = "Sex",
               col = sexColors,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })

ht6 <- Heatmap(sampleClassification, name = "Sample Classification",
               col = sampleClassification_colors,
               na_col = "gray90",
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })
               
               
               
               
ht_combined <-  ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6

# exporting to pdf.
pdf("Output/ExpressionCnMutationFusions.pdf", width = 8, height = 12)
draw(ht_combined, heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(0, 0, 0, 0), "mm"))
dev.off()

# trying different grids here...
ht1g <- grid.grabExpr(draw(ht1 + ht2, newpage = FALSE))
ht2g <- grid.grabExpr(draw(ht3 + ht4, newpage = FALSE))


# Arrange them in a 2x1 vertical layout
pdf("Output/2x2_heatmap_grid.pdf", width = 14, height = 11)
grid.arrange(ht1g, ht2g, ncol = 1,
             top = textGrob("Heatmap For Expression, Copy Number, Mutation and Gene Fusions", gp = gpar(fontsize = 16, fontface = "bold"), just = "center"))
dev.off()







########################################################################################

### Now making heatmaps with ALT & TA genes.
gene_alts <- c("ATRX", "DAXX", "SLX4", "MUS81", "PML", "RAD51", "ASF1B", "CHEK1", "FANCM", "FANCD2", "BLM", "SMARCA5", "CHD4", "BRCA2", "TOP3A")
gene_tas <- c("TERT", "DKC1", "GAR1", "MYC")

gene_combined <-c("ATRX", "DAXX", "SLX4", "MUS81", "PML", "RAD51", "ASF1B", "CHEK1", "FANCM", "FANCD2", "BLM", "SMARCA5", "CHD4", "BRCA2", "TOP3A", "TERT", "DKC1", "GAR1", "MYC")


### scaling, matching order and finally heatmap.
expression_matrix <- t(scale(t(geneExpression[rownames(geneExpression) %in% gene_combined, ,drop = FALSE])))
# replacing sampleID with cell line ID.
colnames(expression_matrix) <- metadata$cell_line_display_name[
  match(colnames(expression_matrix), metadata$depmap_id)]

copy_number_matrix <- t(scale(t(copyNumber[rownames(copyNumber) %in% rownames(expression_matrix), ,drop = FALSE])))
# replacing sampleID with cell line ID.
colnames(copy_number_matrix) <- metadata$cell_line_display_name[
  match(colnames(copy_number_matrix), metadata$depmap_id)]

expression_matrix <- expression_matrix[rownames(expression_matrix) %in% rownames(copy_number_matrix), ,drop = FALSE]

gene_combined <- gene_combined[gene_combined %in% rownames(expression_matrix)]


### gene fusion data.
geneFusion <- readRDS("Output/geneFusion.rds") # from geneFusion.R.

# Converting geneFusion to appropriate format.
geneFusion <- geneFusion %>%
  mutate(Present = 1) %>%
  distinct(Gene, SampleID, .keep_all = TRUE) %>%
  pivot_wider(names_from = SampleID, values_from = Present, values_fill = 0)

geneFusion <- geneFusion %>%
  dplyr::select(2:last_col())

geneFusion <- geneFusion[geneFusion$Gene %in% rownames(expression_matrix), ,drop = FALSE]
missing_genes <- setdiff(rownames(expression_matrix), geneFusion$Gene)

# creating a 0-filled dataframe for the missing genes: not doing this gave error.
zero_rows <- matrix(0, nrow = length(missing_genes), ncol = ncol(geneFusion) - 1) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(geneFusion)[-1]  # match sample column names
zero_rows <- cbind(Gene = missing_genes, zero_rows)

geneFusion <- bind_rows(geneFusion, zero_rows)
geneFusion <- geneFusion %>%
  group_by(Gene) %>%
  summarise(across(everything(), max), .groups = "drop")

geneFusion <- as.data.frame(geneFusion)
rownames(geneFusion) <- geneFusion$Gene
geneFusion$Gene = NULL

# replacing sampleID with cell line ID.
colnames(geneFusion) <- metadata$cell_line_display_name[
  match(colnames(geneFusion), metadata$depmap_id)]

### mutation data.
somaticMutation <- readRDS("Output/somaticMutation.rds")
somaticMutation <- somaticMutation@data
somaticMutation <- somaticMutation %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)

somaticMutation <- somaticMutation %>%
  mutate(Variant_Classification = as.character(Variant_Classification)) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode,
              values_from = Variant_Classification,
              values_fill = "None")

somaticMutation <- as.data.frame(somaticMutation)
rownames(somaticMutation) <- somaticMutation$Hugo_Symbol
somaticMutation$Hugo_Symbol = NULL

somaticMutation <- somaticMutation[rownames(somaticMutation) %in% rownames(copy_number_matrix), ,drop = FALSE]

# creating a "None"-filled dataframe for the missing genes

missing_genes <- setdiff(rownames(expression_matrix), rownames(somaticMutation))
zero_rows <- matrix("None", nrow = length(missing_genes), ncol = ncol(somaticMutation)) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(somaticMutation)
rownames(zero_rows) <- missing_genes

zero_rows <- cbind(zero_rows)

somaticMutation <- bind_rows(somaticMutation, zero_rows)


### matching the rows and columns before making the heatmap.
copy_number_matrix <- copy_number_matrix[rownames(expression_matrix), colnames(expression_matrix)]
somaticMutation <- somaticMutation[rownames(expression_matrix), colnames(expression_matrix)]
geneFusion <- geneFusion[rownames(expression_matrix), colnames(expression_matrix)]

row_split <- ifelse(gene_combined %in% gene_alts, "ALT", "TA") # for heatmap.
row_split <- factor(row_split, levels = c("ALT", "TA"))
names(row_split) <- gene_combined


ht1 <- Heatmap(expression_matrix,
               name = "Expression",
               row_order = names(row_split),
               #row_split = row_split,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 20),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht2 <- Heatmap(copy_number_matrix,
               name = "Copy Number",
               row_order = names(row_split),
              #row_split = row_split,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 20),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht3 <- Heatmap(geneFusion, name = "Fusion", 
               col = fusion_colors, 
               row_order = names(row_split),
               #row_split = row_split,
               cluster_rows = FALSE, 
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })

ht4 <- Heatmap(somaticMutation, name = "Mutations", 
               col = mutation_colors, 
               row_order = names(row_split),
               #row_split = row_split,
               cluster_rows = FALSE, 
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht5 <- Heatmap(sexMetadata, name = "Sex",
               col = sexColors,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })

ht6 <- Heatmap(sampleClassification, name = "Sample Classification",
               col = sampleClassification_colors,
               na_col = "gray90",
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })




ht_combined <-  ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6

# exporting to pdf.
pdf("Output/ExpressionCnMutationFusions2.pdf", width = 6, height = 12)
draw(ht_combined)
dev.off()


####################################################################################

### Now making heatmaps with ALT & TA signature genes.
#gene_alts <- c("ATRX", "DAXX", "SLX4", "MUS81", "PML", "RAD51", "ASF1B", "CHEK1", "FANCM", "FANCD2", "BLM", "SMARCA5", "CHD4", "BRCA2", "TOP3A")
#gene_tas <- c("TERT", "DKC1", "GAR1", "MYC")

gene_combined <-c("GASK1B", "SYCP2", "H2BC11", "MAB21L2", "CLSTN3", "ELOVL4", "MAP9", "COL24A1", "PGM2L1", "STK32A", "CRYAB", "ARMH4", "LFNG", "SYTL5", "CCL28",
                "TERT", "CCDC8", "ALDH1A3", "ERICH1", "TDRP", "SHFL", "DNPH1", "IRX2")


### scaling, matching order and finally heatmap.
expression_matrix <- t(scale(t(geneExpression[rownames(geneExpression) %in% gene_combined, ,drop = FALSE])))
# replacing sampleID with cell line ID.
colnames(expression_matrix) <- metadata$cell_line_display_name[
  match(colnames(expression_matrix), metadata$depmap_id)]

copy_number_matrix <- t(scale(t(copyNumber[rownames(copyNumber) %in% rownames(expression_matrix), ,drop = FALSE])))
# replacing sampleID with cell line ID.
colnames(copy_number_matrix) <- metadata$cell_line_display_name[
  match(colnames(copy_number_matrix), metadata$depmap_id)]

expression_matrix <- expression_matrix[rownames(expression_matrix) %in% rownames(copy_number_matrix), ,drop = FALSE]

gene_combined <- gene_combined[gene_combined %in% rownames(expression_matrix)]


### gene fusion data.
geneFusion <- readRDS("Output/geneFusion.rds") # from geneFusion.R.

# Converting geneFusion to appropriate format.
geneFusion <- geneFusion %>%
  mutate(Present = 1) %>%
  distinct(Gene, SampleID, .keep_all = TRUE) %>%
  pivot_wider(names_from = SampleID, values_from = Present, values_fill = 0)

geneFusion <- geneFusion %>%
  dplyr::select(2:last_col())

geneFusion <- geneFusion[geneFusion$Gene %in% rownames(expression_matrix), ,drop = FALSE]
missing_genes <- setdiff(rownames(expression_matrix), geneFusion$Gene)

# creating a 0-filled dataframe for the missing genes: not doing this gave error.
zero_rows <- matrix(0, nrow = length(missing_genes), ncol = ncol(geneFusion) - 1) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(geneFusion)[-1]  # match sample column names
zero_rows <- cbind(Gene = missing_genes, zero_rows)

geneFusion <- bind_rows(geneFusion, zero_rows)
geneFusion <- geneFusion %>%
  group_by(Gene) %>%
  summarise(across(everything(), max), .groups = "drop")

geneFusion <- as.data.frame(geneFusion)
rownames(geneFusion) <- geneFusion$Gene
geneFusion$Gene = NULL

# replacing sampleID with cell line ID.
colnames(geneFusion) <- metadata$cell_line_display_name[
  match(colnames(geneFusion), metadata$depmap_id)]

### mutation data.
somaticMutation <- readRDS("Output/somaticMutation.rds")
somaticMutation <- somaticMutation@data
somaticMutation <- somaticMutation %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)

somaticMutation <- somaticMutation %>%
  mutate(Variant_Classification = as.character(Variant_Classification)) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode,
              values_from = Variant_Classification,
              values_fill = "None")

somaticMutation <- as.data.frame(somaticMutation)
rownames(somaticMutation) <- somaticMutation$Hugo_Symbol
somaticMutation$Hugo_Symbol = NULL

somaticMutation <- somaticMutation[rownames(somaticMutation) %in% rownames(copy_number_matrix), ,drop = FALSE]

# creating a "None"-filled dataframe for the missing genes

missing_genes <- setdiff(rownames(expression_matrix), rownames(somaticMutation))
zero_rows <- matrix("None", nrow = length(missing_genes), ncol = ncol(somaticMutation)) %>%
  as.data.frame()
colnames(zero_rows) <- colnames(somaticMutation)
rownames(zero_rows) <- missing_genes

zero_rows <- cbind(zero_rows)

somaticMutation <- bind_rows(somaticMutation, zero_rows)


### matching the rows and columns before making the heatmap.
copy_number_matrix <- copy_number_matrix[rownames(expression_matrix), colnames(expression_matrix)]
somaticMutation <- somaticMutation[rownames(expression_matrix), colnames(expression_matrix)]
geneFusion <- geneFusion[rownames(expression_matrix), colnames(expression_matrix)]

row_split <- ifelse(gene_combined %in% gene_alts, "ALT", "TA") # for heatmap.
row_split <- factor(row_split, levels = c("ALT", "TA"))
names(row_split) <- gene_combined


ht1 <- Heatmap(expression_matrix,
               name = "Expression",
               row_order = names(row_split),
               #row_split = row_split,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 20),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht2 <- Heatmap(copy_number_matrix,
               name = "Copy Number",
               row_order = names(row_split),
               #row_split = row_split,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 20),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht3 <- Heatmap(geneFusion, name = "Fusion", 
               col = fusion_colors, 
               row_order = names(row_split),
               #row_split = row_split,
               cluster_rows = FALSE, 
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })

ht4 <- Heatmap(somaticMutation, name = "Mutations", 
               col = mutation_colors, 
               row_order = names(row_split),
               #row_split = row_split,
               cluster_rows = FALSE, 
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(6, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })


ht5 <- Heatmap(sexMetadata, name = "Sex",
               col = sexColors,
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })

ht6 <- Heatmap(sampleClassification, name = "Sample Classification",
               col = sampleClassification_colors,
               na_col = "gray90",
               cluster_rows = FALSE,
               show_column_names = TRUE,
               cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8, fontface = "bold"),
               column_split = metadata$TMMstatus,
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontsize = 6),
               height = unit(0.2, "cm"),
               width = unit(6, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "black", fill = NA, lwd = 0.75))
               })




ht_combined <-  ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6

# exporting to pdf.
pdf("Output/ExpressionCnMutationFusions3.pdf", width = 6, height = 12)
draw(ht_combined)
dev.off()



#####################################################################################

### gistic analysis: making gistic files.

omicsProfile <- read_csv("OmicsProfiles.csv")
omicsProfile <- omicsProfile[omicsProfile$ModelID %in% metadata$depmap_id, ,drop = FALSE]


omicsProfile <- omicsProfile %>%
  arrange(ModelID)

# using segment file from datatype wgs.
gistic_list <- omicsProfile$ProfileID[omicsProfile$ModelID %in% metadata$depmap_id & omicsProfile$Datatype == "wgs"]
gistic_list_C1 <- omicsProfile$ProfileID[omicsProfile$ModelID %in% metadata$depmap_id[metadata$TMMstatus == "ALT"] & omicsProfile$Datatype == "wgs"]
gistic_list_C2 <- omicsProfile$ProfileID[omicsProfile$ModelID %in% metadata$depmap_id[metadata$TMMstatus == "TA"] & omicsProfile$Datatype == "wgs"]


CNsegment <- read_csv("OmicsCNSegmentsProfile.csv")

CNsegmentC1 <- CNsegment[CNsegment$ProfileID %in% gistic_list_C1, ,drop = FALSE]
CNsegmentC2 <- CNsegment[CNsegment$ProfileID %in% gistic_list_C2, ,drop = FALSE]
CNsegment <- CNsegment[CNsegment$ProfileID %in% gistic_list, ,drop = FALSE]

write_csv(CNsegmentC1, "Output/segmentC1.csv")
write_csv(CNsegmentC2, "Output/segmentC2.csv")
write_csv(CNsegment,  "Output/CNsegment.csv")

########################################################





