library(tidyverse)
library(biomaRt)
library(stringr)


# Loading necessary rds files.
fusion <- readRDS("OmicsFusionFiltered.RDS")
metadata <- read_delim("Osteosarcoma_CellLines.txt")

right_counts <- table(fusion$RightGene) # fusion while gene is on the right side of the fusion.
left_counts <- table(fusion$LeftGene) # fusion while gene is on the left side of the fusion.

all_genes <- unique(union(names(left_counts), names(right_counts)))

# Combining counts: eventually add left and right to make total fusions.
gene_counts <- data.frame(
  Gene = all_genes,
  CountInRight = as.integer(right_counts[all_genes]),
  CountInLeft = as.integer(left_counts[all_genes])
)

# Replacing NAs with 0 for genes that were only in one side
gene_counts$CountInLeft[is.na(gene_counts$CountInLeft)] <- 0
gene_counts$CountInRight[is.na(gene_counts$CountInRight)] <- 0

# Calculate total
gene_counts$Total <- gene_counts$CountInLeft + gene_counts$CountInRight
gene_counts <- gene_counts[order(-gene_counts$Total), ]

# top 10 genes with most fusions.
gene_counts_10 <- gene_counts %>%
  slice_head(n = 10)

gene_counts_10 <- gene_counts_10 %>%
  mutate(
    Gene = str_extract(Gene, "^[^ ]+") # gene names also contained Ensembl ID. Removing those for simplicity.
  )
  

ggplot(gene_counts_10, aes(x = reorder(Gene, -Total), y = Total)) + geom_bar(stat = "identity") +
         theme_classic() +
         labs(x = "Gene", y = "Total Fusions", title = paste("Top 10", "Genes by Fusion Count")) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
# this graph is faulty because of presence of duplicates (same fusion in the same sample).

#####################################################################################

### removing the duplicates (same fusion event in the same sample).
fusion_main <- fusion %>%
  mutate(
    LeftGene = str_extract(LeftGene, "^[^ ]+"),
    RightGene = str_extract(RightGene, "^[^ ]+")
  ) %>%
  distinct(SampleID, FusionName, LeftGene, RightGene)

## now counting mutations in different class for a bar graph.
fusion_annotated <- fusion_main %>%
  left_join(metadata, by = c("SampleID" = "depmap_id"))

fusion_annotated <- fusion_annotated %>%
  pivot_longer(cols = c(LeftGene, RightGene), names_to = "Side", values_to = "Gene") %>%
 dplyr::select(SampleID, TMMstatus, Gene)

saveRDS(fusion_annotated, "Output/geneFusion.rds")

# pivoting fusion_annotated to see number the number of gene fusions for a gene across the two phenotypes.
fusion_class_counts <- fusion_annotated %>%
  dplyr::select(Gene, TMMstatus) %>%
  count(Gene, TMMstatus) %>%
  pivot_wider(names_from = TMMstatus, values_from = n, values_fill = 0)
  
right_counts <- table(fusion_main$RightGene)
left_counts <- table(fusion_main$LeftGene)

all_genes <- unique(union(names(left_counts), names(right_counts)))

# Combine counts
gene_counts <- data.frame(
  Gene = all_genes,
  CountInRight = as.integer(right_counts[all_genes]),
  CountInLeft = as.integer(left_counts[all_genes])
)

# Replacing NAs with 0 for genes that were only in one side
gene_counts$CountInLeft[is.na(gene_counts$CountInLeft)] <- 0
gene_counts$CountInRight[is.na(gene_counts$CountInRight)] <- 0

# Calculate total
gene_counts$Total <- gene_counts$CountInLeft + gene_counts$CountInRight
gene_counts <- gene_counts[order(-gene_counts$Total), ]
gene_counts_10 <- gene_counts %>%
  slice_head(n = 10)


gene_counts_10 <- gene_counts_10 %>%
  left_join (fusion_class_counts, by = "Gene")

### stacked barplot.
gene_counts_long <- gene_counts_10 %>%
  dplyr::select(Gene, ALT, TA) %>%
  pivot_longer(cols = c("ALT", "TA"), names_to = "Class", values_to = "FusionCount")

ggplot(gene_counts_long, aes(x = reorder(Gene, -FusionCount), y = FusionCount, fill = Class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("ALT" = "blue", "TA" = "peachpuff4")) +
  labs(x = "Gene", y = "Fusion Count", fill = "Class",
       title = "Stacked Barplot of Gene Fusions by Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5))


###################################################################################

### separating files based on their classes.
fusionC1 <- fusion[fusion$SampleID %in% metadata$depmap_id[metadata$TMMstatus == "ALT"], ]
fusionC2 <- fusion[fusion$SampleID %in% metadata$depmap_id[metadata$TMMstatus == "TA"], ]

View(head(sort(table(fusionC1$FusionName), decreasing = TRUE), n = 11))
View(head(sort(table(fusionC2$FusionName), decreasing = TRUE), n = 14))

View(sort(table(fusionC1$SampleID), decreasing = TRUE))
View(sort(table(fusionC2$SampleID), decreasing = TRUE))


##### Differential gene fusion. Using Fisher's exact test.

fusionC1 <- fusionC1 %>%
  mutate(
    LeftGene = str_extract(LeftGene, "^[^ ]+"),
    RightGene = str_extract(RightGene, "^[^ ]+")
  ) %>%
  distinct(SampleID, FusionName, LeftGene, RightGene) # for this test, we will be counting number of samples a particular gene has fused on. If TP53 fuses twice on the same sample, only counted as 1 for this test.

fusionC2 <- fusionC2 %>%
  mutate(
    LeftGene = str_extract(LeftGene, "^[^ ]+"),
    RightGene = str_extract(RightGene, "^[^ ]+")
  ) %>%
  distinct(SampleID, FusionName, LeftGene, RightGene)

View(head(sort(table(fusionC1$FusionName), decreasing = TRUE), n = 11))
View(head(sort(table(fusionC2$FusionName), decreasing = TRUE), n = 14))

# the highest occurence of a fusion is 2 for a phenotype. So, obviously none of the fusion events will be signficiant for one phenotype over the other.

# bar graph showing distribution of fusion across samples.
sampleCountC1 <- sort(table(fusionC1$SampleID), decreasing = TRUE)
sampleCountC1 <- as.data.frame(sampleCountC1)
colnames(sampleCountC1) <- c("SampleID", "FusionCount")
sampleCountC1$Class <- "ALT"

sampleCountC2 <- sort(table(fusionC2$SampleID), decreasing = TRUE)
sampleCountC2 <- as.data.frame(sampleCountC2)
colnames(sampleCountC2) <- c("SampleID", "FusionCount")
sampleCountC2$Class <- "TA"


sampleCountCombined <- rbind(sampleCountC1, sampleCountC2)


ggplot(sampleCountCombined, aes(x = SampleID, y = FusionCount, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Fusion Counts per Sample",
       x = "Sample ID",
       y = "Number of Fusions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


### Instead of looking at a complete fusion event, looking at differential presence of a particular gene in fusion events.
## for ALT phenotype.
right_counts <- table(fusionC1$RightGene)
left_counts <- table(fusionC1$LeftGene)

all_genes <- unique(union(names(left_counts), names(right_counts)))

# Combine counts
gene_counts_C1 <- data.frame(
  Gene = all_genes,
  CountInRight = as.integer(right_counts[all_genes]),
  CountInLeft = as.integer(left_counts[all_genes])
)

# Replacing NAs with 0 for genes that were only in one side
gene_counts_C1$CountInLeft[is.na(gene_counts_C1$CountInLeft)] <- 0
gene_counts_C1$CountInRight[is.na(gene_counts_C1$CountInRight)] <- 0

# Calculate total
gene_counts_C1$Total <- gene_counts_C1$CountInLeft + gene_counts_C1$CountInRight
gene_counts_C1 <- gene_counts_C1[order(-gene_counts_C1$Total), ]



## for TA phenotype.
right_counts <- table(fusionC2$RightGene)
left_counts <- table(fusionC2$LeftGene)

all_genes <- unique(union(names(left_counts), names(right_counts)))

# Combine counts
gene_counts_C2 <- data.frame(
  Gene = all_genes,
  CountInRight = as.integer(right_counts[all_genes]),
  CountInLeft = as.integer(left_counts[all_genes])
)

# Replacing NAs with 0 for genes that were only in one side
gene_counts_C2$CountInLeft[is.na(gene_counts_C2$CountInLeft)] <- 0
gene_counts_C2$CountInRight[is.na(gene_counts_C2$CountInRight)] <- 0

# Calculate total
gene_counts_C2$Total <- gene_counts_C2$CountInLeft + gene_counts_C2$CountInRight
gene_counts_C2 <- gene_counts_C2[order(-gene_counts_C2$Total), ]

##### differential fusion event. Running Fisher's exact test.

fusion_filtered <- fusion %>%
  mutate(
    LeftGene = str_extract(LeftGene, "^[^ ]+"),
    RightGene = str_extract(RightGene, "^[^ ]+")
  ) %>%
  distinct(SampleID, FusionName, LeftGene, RightGene)
  

# Combining left and right fusion partners into a long format.
fusion_long <- fusion_filtered %>%
  dplyr::select(SampleID, LeftGene, RightGene) %>%
  pivot_longer(cols = c(LeftGene, RightGene), names_to = "Side", values_to = "Gene") %>%
  distinct(SampleID, Gene)


# binary presence matrix for Fisher's exact test.
fusion_matrix <- fusion_long %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = Gene, values_from = Present, values_fill = 0)

fusion_matrix <- fusion_matrix %>%
  left_join(metadata, by = c("SampleID" = "depmap_id"))

fusion_matrix <- fusion_matrix %>%
  dplyr::select(SampleID, TMMstatus, everything())
fusion_matrix$TMMstatus <- factor(fusion_matrix$TMMstatus, levels = c("ALT", "TA"))


gene_list <- colnames(fusion_matrix)[3:ncol(fusion_matrix)]

results <- list()

# making 2x2 contigency table to run the test.
for (gene in gene_list) {
  table_gene <- table(fusion_matrix[[gene]], fusion_matrix$TMMstatus)

  if (all(dim(table_gene) == c(2, 2))) {
    test <- fisher.test(table_gene)
    results[[gene]] <- c(
      p_value = test$p.value
    )
  }
}

df_results <- do.call(rbind, results) %>% 
  as.data.frame()

df_results$Gene <- rownames(df_results)
df_results <- df_results %>% 
  arrange(p_value)

df_results$FDR <- p.adjust(df_results$p_value, method = "BH")

## From these results, none of the genes are significiantly fused in one phenotype over the other.


##### genes frequently mutated in either phenotypes.
fusion_summary <- data.frame(Gene = gene_list)

# Count how many samples in each class have a fusion
fusion_summary$NumberOfAltSamplesWithFusion <- sapply(gene_list, function(gene) {
  sum(fusion_matrix[[gene]] == 1 & fusion_matrix$TMMstatus == "ALT")
})

fusion_summary$NumberOfTAsamplesWithFusion <- sapply(gene_list, function(gene) {
  sum(fusion_matrix[[gene]] == 1 & fusion_matrix$TMMstatus == "TA")
})

#Filter for fusions present in both
fusion_summary_shared <- fusion_summary %>%
  filter(NumberOfAltSamplesWithFusion >= 1 & NumberOfTAsamplesWithFusion >= 1)

# rank by average or minimum count.
fusion_summary_shared <- fusion_summary_shared %>%
  mutate(AvgCount = (NumberOfTAsamplesWithFusion + NumberOfAltSamplesWithFusion) / 2) %>%
  arrange(desc(AvgCount)) 


# TP53 the most common fusion across both phenotypes.


