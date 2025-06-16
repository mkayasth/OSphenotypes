library(tidyverse)
library(biomaRt)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(patchwork)
library(EXTEND)


##### Loading bulk-seq RNA & metadata file.
geneExpression <- readRDS("Output/geneExpression.rds")  #from rnaSeq.R.
metadata <- read_delim("Osteosarcoma_CellLines.txt")

##### running extend scores for the two phenotypes.
extendScore <- RunEXTEND(geneExpression)
file.rename("TelomeraseScores.txt", "Output/geneExpression_telomeraseScores.txt")
geneExpression_telomerase_scores <- read.table("Output/geneExpression_telomeraseScores.txt", header = TRUE, row.names = 1)

# in sample names, the - has been replaced by .
rownames(geneExpression_telomerase_scores) <- gsub("\\.", "-", rownames(geneExpression_telomerase_scores))

geneExpression_telomerase_scores$Class <- metadata$TMMstatus[match(rownames(geneExpression_telomerase_scores), metadata$depmap_id)]
ggplot(geneExpression_telomerase_scores, aes(x = reorder(rownames(geneExpression_telomerase_scores), NormEXTENDScores), y = NormEXTENDScores, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "EXTEND Scores Across OS Samples",
       x = "Sample ID",
       y = "Normalized EXTEND Scores") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# boxplot.
model <- lm(geneExpression_telomerase_scores$NormEXTENDScores ~ geneExpression_telomerase_scores$Class, data = geneExpression_telomerase_scores)
summary_model <- summary(model)

boxplot(geneExpression_telomerase_scores$NormEXTENDScores ~ geneExpression_telomerase_scores$Class, ylab = "EXTEND Scores", xlab = "Phenotype")
legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                  "\np-value =", round(summary_model$coefficients[2, 4], 4)),
       bty = "n", cex = 0.8)

rm(summary_model)
rm(model)


#################################################################################

### correlating EXTEND score with TERT expression for the samples.
geneExpressionTERT <- geneExpression[rownames(geneExpression) == "TERT", , drop = FALSE]
geneExpression_telomerase_scores <- read.table("Output/geneExpression_telomeraseScores.txt", header = TRUE, row.names = 1)
rownames(geneExpression_telomerase_scores) <- gsub("\\.", "-", rownames(geneExpression_telomerase_scores))
geneExpression_telomerase_scores$Class <- metadata$TMMstatus[match(rownames(geneExpression_telomerase_scores), metadata$depmap_id)]


geneExpressionTERT <- as.data.frame(scale(t(geneExpressionTERT)))
common_rownames <- intersect(rownames(geneExpressionTERT), rownames(geneExpression_telomerase_scores))
geneExpressionTERT <- geneExpressionTERT[common_rownames, , drop = FALSE]
geneExpression_telomerase_scores <- geneExpression_telomerase_scores[common_rownames, , drop = FALSE]
geneExpression_telomerase_scores <- cbind(geneExpression_telomerase_scores, geneExpressionTERT)
geneExpression_telomerase_scores <- sort_by(geneExpression_telomerase_scores, geneExpression_telomerase_scores$Class)

# scatterplot TERT vs. Extend.

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(geneExpression_telomerase_scores$NormEXTENDScores, geneExpression_telomerase_scores$TERT, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(geneExpression_telomerase_scores, aes(y= TERT, x= NormEXTENDScores, colour = Class)) + 
  geom_point(size = 3, alpha = 0.9) + 
  theme_classic() + geom_smooth(method = "lm", se = TRUE, color = "white", fill = "lightgray") + 
  annotate("text", x = Inf, y = Inf, label = r_text,
          hjust = 1.1, vjust = 1.5, size = 4)

#################################################################################

# # gene expression pathway analysis playground.
# 
# expr_data <- read.csv("rna-seq.csv", row.names = 1, check.names = FALSE) # from rnaSeq.R.
# 
# # Creating GCT header
# gct_header <- data.frame(Name = rownames(expr_data),
#                          Description = "NA",
#                          expr_data,
#                          check.names = FALSE)
# 
# 
# gct_file <- "rna-seq.gct"
# con <- file(gct_file, "w")
# 
# # header lines
# writeLines("#1.2", con)
# writeLines(paste(nrow(expr_data), ncol(expr_data), sep = "\t"), con)
# 
# 
# write.table(gct_header, con, sep = "\t", row.names = FALSE, quote = FALSE)
# 
# 
# close(con)

# now for making cls file.
# samples <- metadata$depmap_id
# phenotypes <- metadata$TMMstatus
# 
# # CLS output path
# cls_file <- "phenotype.cls"
# 
# con <- file(cls_file, "w")
# writeLines(paste(length(samples), length(unique(phenotypes)), 1), con)
# writeLines(paste("#", paste(unique(phenotypes), collapse = " ")), con)
# writeLines(paste(phenotypes, collapse = " "), con)
# close(con)
# 



