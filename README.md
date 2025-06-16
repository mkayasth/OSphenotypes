# rna_Seq.R:

-	Filtered for genes: logTPM > 1 in at least 5% of the sample.
-	Differential gene expression using Limma (p-value < 0.05 & logFC > 1.5): 171 candidates. 125 differentially upregulated in ALT and 46 differentially upregulated in TA.
o	Volcano plot (TA vs. ALT).
-	T-test (p-value < 0.01) followed by regression test (p-value < 0.01 & R2 > 0.3):  51 candidates: 39 ALT and 12 TA.
o	For ALT and TA genes: heatmap + GSVA score.
-	ssGSEA Test: making all possible combination of 3, 4, 5 and 6 genes and running ssGSEA to see what combination of genes gives the best distinction between ALT and TA groups.
-	ssGSEA test (single sample gene enrichment method) to identify signature genes for each phenotype using R package GSVA from Bioconductor.
-	Made all possible combinations of 3 genes, 4 genes, 5 genes and 6 genes from the 51 differentially expressed genes using Combinat package in R.
-	Since we have a lot of genes overexpressed in ALT, there are a lot of possible combinations of ALT genes. Calculating mean difference between ALT and TA for the gene combinations and taking the top 50,000 with highest difference.
-	For each combination, performing t-test (p-value < 0.01) followed by regression test (p-value < 0.01): arranging by decreasing R2 value. Running gsva with the signature.
o	Bar graph + boxplot of ssGSEA & GSVA.

# somaticMutation.R:
-	Making maf file: same gene may have different mutations. Making different headers for them by separating with & into different lines.
-	Changing variant classification and variant type into maftools friendly format.
-	Make ALT and TA maf files.
o	maf summary for ALT vs. TA.
o	Sample Summary and Gene Summary
o	Transition vs. Transversion.
o	drug interaction.
o	heatmap for mutational signature.
o	Mutation Load.
o	Co-oncoplot + co-Barplot.
o	Oncogenic Pathways
o	Mutational Signatures in ALT vs. TA: using NMF to extract de novo mutational signature and comparing to SBS signatures. For TA, only keeping mutations that appears in at least 10% of the samples to avoid rare mutations interference while running NMF.


# copyNumber.R:
-	Limma for differential copy Number (p-value < 0.05, logFC > 0.75); volcano plot.  36 ALT, 18 Telomerase.
-	T-test (p-value < 0.01) and R2 > 0.3: 4 ALT, 18 Telomerase.
-	Looking at differential gene expression of differential copy number genes.
-	Integrated analysis of OS and ALT+TA genes (heatmaps of expression, copy number, fusion, mutation..)
-	Preparing for gistic files – gistic continued in cnvBurden.R.

# geneFusion.R:
-	Removed the same fusion event counted multiple times (same sample, same fusion event, almost same starting and ending point counted as 1 event).
o	Stacked bar plot for fusion across both phenotypes/ individually phenotypes.
-	Fisher’s test to show differential presence of fusion in one phenotype over the other: no significantly differential fusion.
-	Shared fusion: table showing fusions common across the samples of both phenotypes.

# cnvBurden.R:
-	Making GISTIC format ready files: removing sex chromosomes because they had a lot of U (unknown) values.
-	This data is from GATK and centers around 1. Gistic takes in data centered around 0. So, taking a log2. Also, our data is already labeled as +, - or 0. Based on this, deciding on threshold while converting to gistic file format: log2(0.85) and log2(1.13).
-	CNV Burden Scores: aneuploidy scores, number of segments, Fraction of genome altered.
o	Bar graphs + boxplot to show difference between two phenotypes for these CNV Burden Scores.
-	Gistic files: finding similarities and differences between ALT and TA samples. Run on output from gistic.

# EXTEND.R:
-	running extend scores for the two phenotypes.
o	Boxplot + bar plot to show difference in EXTEND score between two phenotypes.
-	correlating EXTEND score with TERT expression for the samples.
o	Using scatterplot with correlation value.








