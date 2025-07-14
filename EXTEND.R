library(tidyverse)
library(biomaRt)
library(ComplexHeatmap)
library(ggpubr)
library(gghalves)
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
  labs(x = "Sample ID",
       y = "Normalized EXTEND Scores") +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
  axis.title.x = element_text(size = 24),
  axis.title.y = element_text(size = 24),
  axis.text.y  = element_text(size = 18, face = "bold"),
  legend.title = element_text(size = 18, face = "bold"),
  legend.text  = element_text(size = 16)
)



# boxplot.
model <- lm(geneExpression_telomerase_scores$NormEXTENDScores ~ geneExpression_telomerase_scores$Class, data = geneExpression_telomerase_scores)
summary_model <- summary(model)

boxplot(geneExpression_telomerase_scores$NormEXTENDScores ~ geneExpression_telomerase_scores$Class, ylab = "EXTEND Scores", xlab = "Phenotype")
legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                  "\np-value =", round(summary_model$coefficients[2, 4], 4)),
       bty = "n", cex = 0.8)

# violin plot.

# First, extracting the text for violin plot.
r_text <- paste0("R² = ", round(summary_model$r.squared, 3))

ggplot(geneExpression_telomerase_scores, aes(x = Class, y = NormEXTENDScores, fill = Class, color = Class)) +
  geom_violin(size=0.2, draw_quantiles=c(0.5),alpha=0.5) +
  geom_point(position = position_jitter(width = .08), size = 3)  + scale_fill_manual(values=c("ALT"="lightblue","TA"="lightpink2")) + 
  scale_color_manual(values=c("ALT"="blue", "TA"="darkred"))+ theme_classic() + 
  labs(x = "Class", y = "Normalized EXTEND Scores") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16, face = "bold"), legend.position = "none") +
  stat_compare_means(comparisons = list(c("ALT","TA")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 8, tip.length = 0.01,
                     label.y = 1.05) +  
  annotate("text", x = 2, y = 0.3,
            label = r_text, size = 8, fontface = "bold")  +
coord_cartesian(ylim = c(0, 1.1))


# violin plot with boxplot.
ggplot(geneExpression_telomerase_scores, aes(x = Class, y = NormEXTENDScores)) +
  # Half violin.
  geom_half_violin(
    aes(fill = Class),
    side = "r",
    alpha = 0.5,
    trim = TRUE
  ) +
  # Boxplot on left.
  geom_boxplot(
    aes(fill = Class),
    width = 0.1,
    position = position_nudge(x = -0.07),
    outlier.shape = NA,
    alpha = 0.5
  ) +
  # Jittered points.
  geom_point(
    aes(color = Class),
    position = position_jitter(width = 0.1),
    size = 3,
    alpha = 1
  ) +
  scale_fill_manual(values = c("ALT" = "lightblue", "TA" = "lightpink2")) +
  scale_color_manual(values = c("ALT" = "darkblue", "TA" = "darkred")) +
  theme_classic() +
  labs(x = "Class", y = "Normalized EXTEND Scores") +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none") +
  stat_compare_means(comparisons = list(c("ALT","TA")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 8, tip.length = 0.01,
                     label.y = 1.05) +  
  annotate("text", x = 2, y = 0.3,
           label = r_text, size = 8, fontface = "bold")  +
  coord_cartesian(ylim = c(0, 1.1))

rm(summary_model)
rm(r_text)
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

ggplot(geneExpression_telomerase_scores, aes(y = TERT, x = NormEXTENDScores, colour = Class)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "TA" = "darkred")) + 
  geom_point(size = 6, alpha = 0.9) + 
  labs(x = "Normalized EXTEND Scores", y = "TERT expression") + 
  geom_smooth(method = "lm", level = 0.90, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 2, label = r_text, size = 10) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(0.9, 0.2)
  )

 


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



