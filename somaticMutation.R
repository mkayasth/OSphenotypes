library(tidyverse)
library(biomaRt)
library(grid)
library(maftools)
library("mclust")

# Loading somatic mutation & copy number rds.
somaticMutation <- readRDS("OmicsSomaticMutations.RDS")
metadata <- read_delim("Osteosarcoma_CellLines.txt")


### making maf object.

maf_file <- data.frame(
  Hugo_Symbol = somaticMutation$HugoSymbol,
  Chromosome = somaticMutation$Chrom,
  Start_Position = somaticMutation$Pos,
  End_Position = somaticMutation$Pos,
  Reference_Allele = somaticMutation$Ref,
  Tumor_Seq_Allele2 = somaticMutation$Alt,
  Variant_Type = somaticMutation$VariantType,
  Variant_Classification = somaticMutation$VariantInfo,
  Tumor_Sample_Barcode = somaticMutation$ModelID,
  AF = somaticMutation$AF,
  RefCount = somaticMutation$RefCount,
  AltCount = somaticMutation$AltCount,
  t_vaf = somaticMutation$AltCount / (somaticMutation$AltCount + somaticMutation$RefCount),
  Exon = somaticMutation$Exon,
  Intron = somaticMutation$Intron,
  GnomadeAF = somaticMutation$GnomadeAF,
  GnomadgAF = somaticMutation$GnomadgAF,
  Polyphen = somaticMutation$Polyphen,
  Sift = somaticMutation$Sift, 
  RevelScore = somaticMutation$RevelScore,
  OncogeneHighImpact = somaticMutation$OncogeneHighImpact,
  TumorSuppressorHighImpact = somaticMutation$TumorSuppressorHighImpact, 
  Hotspot = somaticMutation$Hotspot,
  Protein_Change = somaticMutation$ProteinChange
)


# Multiple genes have multiple mutation variant classification. Changing each to a unique header.
maf_file <- maf_file %>%
  separate_rows(Variant_Classification, sep = "&")

# Mapping to MAF-compatible values.
maf_file <- maf_file %>%
  mutate(Variant_Type = case_when(
    Variant_Type == "SNV" ~ "SNP",
    Variant_Type == "insertion" ~ "INS",
    Variant_Type == "deletion" ~ "DEL",
    Variant_Type == "substitution" & nchar(Reference_Allele) == 2 ~ "DNP",
    Variant_Type == "substitution" & nchar(Reference_Allele) == 3 ~ "TNP",
    Variant_Type == "substitution" & nchar(Reference_Allele) > 3 ~ "ONP",
    TRUE ~ Variant_Type  # Keep existing value if none of the above match.
  ))

maf_file <- maf_file %>%
  mutate(Variant_Classification = case_when(
    Variant_Classification == "missense_variant" ~ "Missense_Mutation",
    Variant_Classification == "frameshift_variant" & Variant_Type == "DEL" ~ "Frame_Shift_Del",
    Variant_Classification == "frameshift_variant" & Variant_Type == "INS" ~ "Frame_Shift_Ins",
    Variant_Classification == "stop_gained" ~ "Nonsense_Mutation",
    Variant_Classification == "splice_region_variant" ~ "Splice_Region",
    Variant_Classification == "inframe_deletion" ~ "In_Frame_Del",
    Variant_Classification == "splice_acceptor_variant" ~ "Splice_Site", 
    Variant_Classification == "splice_polypyrimidine_tract_variant" ~ "Splice_Region",
    Variant_Classification == "intron_variant" ~ "Intron",
    Variant_Classification == "non_coding_transcript_variant" ~ "RNA",
    Variant_Classification == "inframe_insertion" ~ "In_Frame_Ins",
    Variant_Classification == "splice_donor_variant" ~ "Splice_Site",
    Variant_Classification == "stop_lost" ~ "Nonstop_Mutation",
    Variant_Classification == "start_lost" ~ "Translation_Start_Site",
    Variant_Classification == "splice_donor_5th_base_variant" ~ "Splice_Region",
    Variant_Classification == "coding_sequence_variant" ~ "Coding_Sequence",
    Variant_Classification == "5_prime_UTR_variant" ~ "5'UTR",
    Variant_Classification == "splice_donor_region_variant" ~ "Splice_Region", 
    Variant_Classification == "3_prime_UTR_variant" ~ "3'UTR",
    Variant_Classification == "protein_altering_variant" ~ "Missense_Mutation",
    Variant_Classification == "non_coding_transcript_exon_variant" ~ "RNA",
    TRUE ~ Variant_Classification
  ))


# exporting a maf file.
write.table(maf_file, file = "Output/output.maf", sep = "\t", quote = FALSE, row.names = FALSE)

# c1 for ALT and c2 for TA phenotypes.
c1_samples <- metadata$depmap_id[metadata$TMMstatus == "ALT"]
c2_samples <- metadata$depmap_id[metadata$TMMstatus == "TA"]

maf_c1 <- maf_file[maf_file$Tumor_Sample_Barcode %in% c1_samples, ]
maf_c2 <- maf_file[maf_file$Tumor_Sample_Barcode %in% c2_samples, ]

write.table(maf_c1, file = "Output/ALTsamples.maf", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(maf_c2, file = "Output/TAsamples.maf", sep = "\t", quote = FALSE, row.names = FALSE)

# maf includes all samples.
maf <- read.maf("Output/output.maf")
plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
getSampleSummary(maf)

saveRDS(maf, "Output/somaticMutation.rds")

geneSummary <- getGeneSummary(maf)
ecdfGeneSummary <- ecdf(geneSummary$MutatedSamples)

# Subset of genes to highlight: genes common in cancer and osteosarcoma from literature.
highlight_genes <- unique(c("TP53","RB1",
                            "ATRX",
                            "MDM2",
                            "TERT","MYCN", "DLG2", "PTPRQ", "ZFHX4",
                            "TP53", "ALK", "FLG", "KMT2D", "EWSR1", "RB1", "WT1", "CDKN2A", "KMT2B", "NF1",
                            "MAP3K4", "ATRX", "CIC", "ASPCR1", "STAG2", "IGF1R", "PTEN", "CDK4", "PDGFRA2",
                            "MYC", "BRAF", "KRAS", "NRAS", "IGF1R", "KMT2C", "PALB2", "CCNE1", "COPS3", "DLEU1", "KDR"))

# Adding a color flag.
geneSummary$Group <- ifelse(geneSummary$Hugo_Symbol %in% highlight_genes, "Highlight", "Other")

# Plotting line graph: this shows the number of samples with mutations in these genes compared to the most mutated genes in the sample.
# compared to the most mutated genes, these genes are not very significant in terms of frequency of mutation. But the top mutation genes could probably be artifact/noise given not much is reported about them in the literature.
ggplot(geneSummary, aes(x = Hugo_Symbol, y = MutatedSamples, color = Group)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("Highlight" = "black", "Other" = scales::alpha("gray", 0.8))) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Gene Expression with Highlighted Genes", y = "Number of samples with Mutations", x = "Genes")


# c1: alt samples.
mafC1 <- read.maf("Output/ALTsamples.maf")
plotmafSummary(mafC1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
getSampleSummary(mafC1)
geneSummaryC1 <- getGeneSummary(mafC1)

# Adding a color flag
geneSummaryC1$Group <- ifelse(geneSummaryC1$Hugo_Symbol %in% highlight_genes, "Highlight", "Other")

# Plotting line graph: most mutated samples compared to most common OS genes.
ggplot(geneSummaryC1, aes(x = Hugo_Symbol, y = MutatedSamples, color = Group)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("Highlight" = "black", "Other" = scales::alpha("gray", 0.8))) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Gene Expression with Highlighted Genes", y = "Number of samples with Mutations", x = "Genes")


# C2: TA samples.
mafC2 <- read.maf("Output/TAsamples.maf")
plotmafSummary(mafC2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
getSampleSummary(mafC2)
geneSummaryC2 <- getGeneSummary(mafC2)

# Adding a color flag
geneSummaryC2$Group <- ifelse(geneSummaryC2$Hugo_Symbol %in% highlight_genes, "Highlight", "Other")

# Plotting line graph.
ggplot(geneSummaryC2, aes(x = Hugo_Symbol, y = MutatedSamples, color = Group)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("Highlight" = "black", "Other" = scales::alpha("gray", 0.8))) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Gene Expression with Highlighted Genes", y = "Number of samples with Mutations", x = "Genes")


### drawing oncoplot for top 20 genes.
pdf("output/OncoplotAllSamples.pdf")
par(oma = c(1, 1, 1, 1))
oncoplot(maf = maf, top = 20, sampleOrder = metadata$depmap_id)
dev.off()

oncoplot(maf = mafC1, top = 20)
oncoplot(maf = mafC2, top = 20)


# transitions and transversions for SNP.
maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
mafC1.titv = titv(maf = mafC1, plot = FALSE, useSyn = TRUE)
mafC2.titv = titv(maf = mafC2, plot = FALSE, useSyn = TRUE)

# plotting titv summary
pdf("output/titvAllSamples.pdf")
par(oma = c(1, 1, 1, 1))
plotTiTv(res = maf.titv)
dev.off()

pdf("output/titvALT.pdf")
par(oma = c(1, 1, 1, 1))
plotTiTv(res = mafC1.titv)
dev.off()

pdf("output/titvTA.pdf")
par(oma = c(1, 1, 1, 1))
plotTiTv(res = mafC2.titv)
dev.off()


# comparing mutation load against TCGA cohorts: pediatric vs adult samples?
pdf("output/mutLoad.pdf", width = 8, height = 5)
maf.mutload = tcgaCompare(maf = c(mafC1, mafC2), cohortName = c('OS-ALT', 'OS-TA'), logscale = TRUE, capture_size = 50)
dev.off()



# plotting vaf.
plotVaf(maf = maf, vafCol = 'i_TumorVAF_WU')
dev.off()
plotVaf(maf = mafC1, vafCol = 'i_TumorVAF_WU')
plotVaf(maf = mafC2, vafCol = 'i_TumorVAF_WU')

### somatic interaction analysis.
#exclusive/co-occurance event analysis on top 25 mutated genes.
pdf("output/somaticInteractions.pdf", width = 5, height = 7)
par(oma = c(1, 1, 1, 1))
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()

pdf("output/somaticInteractionsALT.pdf", width = 5, height = 7)
par(oma = c(1, 1, 1, 1))
somaticInteractions(maf = mafC1, top = 25, pvalue = c(0.05, 0.1))
dev.off()

pdf("output/somaticInteractionsTA.pdf", width = 5, height = 7)
par(oma = c(1, 1, 1, 1))
somaticInteractions(maf = mafC2, top = 25, pvalue = c(0.05, 0.1))
dev.off()


# co-oncoplot: top 20 most mutated genes in ALT and TA (most overlap so total is not 40).
genes = c("MUC19", "MUC4", "C17orf97", "KRTAP9-6", "LILRB3", "MT-CYB", "MT-ND5", "MUC16", "MUC21", "MUC3A", "MYOM3", "TTN", "ADAM29", "AHNAK2", "BTAF1", "CCDC88B", "EP400", "KMT2D", "MUC6", "XIRP2", "HLA-DRB1", "CACNA1B", "DSPP", "FCGBP", "MDC1", "RB1", "ABCA1", "ARHGAP17", "ARHGEF18", "FAM8A1", "OGFR", "USH2A")
coOncoplot(m1 = mafC1, m2 = mafC2, m1Name = 'ALT', m2Name = 'TA', genes = genes, removeNonMutated = TRUE)
pdf("output/coBarPlot.pdf")
par(oma = c(1.5, 1.5, 1.5, 1.5))
coBarplot(m1 = mafC1, m2 = mafC2, m1Name = "ALT", m2Name = "TA", genes = genes)
dev.off()


# co-oncoplot with genes known to mutate in OS from literature.
genes = unique(c("TP53","RB1",
                 "ATRX",
                 "MDM2",
                 "TERT","MYCN", "DLG2", "PTPRQ", "ZFHX4",
                 "TP53", "ALK", "FLG", "KMT2D", "EWSR1", "RB1", "WT1", "CDKN2A", "KMT2B", "NF1",
                 "MAP3K4", "ATRX", "CIC", "ASPCR1", "STAG2", "IGF1R", "PTEN", "CDK4", "PDGFRA2",
                 "MYC", "BRAF", "KRAS", "NRAS", "IGF1R", "KMT2C", "PALB2", "CCNE1", "COPS3", "DLEU1", "KDR", "DAXX"))

coOncoplot(m1 = mafC1, m2 = mafC2, m1Name = 'ALT', m2Name = 'TA', genes = genes, removeNonMutated = FALSE)

# coOncoplot, in default settings, removes genes with 0 mutation across the two phenotypes. coBarplot does not. Manually doing so.
genes_filtered <- genes[genes %in% mafC1@data$Hugo_Symbol | genes %in% mafC2@data$Hugo_Symbol]
coBarplot(m1 = mafC1, m2 = mafC2, m1Name = "ALT", m2Name = "TA", genes = genes)

# co-oncoplot with genes known to mutate in ALT/TA from literature.
genes = c("ATRX", "DAXX", "FANCM", "BLM", "PML", "FANCD2", "UBR5", "RAD21", "SENP5", "NSMCE2", "RFC4", "MCPH1", "RECQL4", "ARID1A", "CCT5", "LRATD2", "RAD1", "SP100", "TERT", "TERC", "POT1", "TINF2", "DKC1")

coOncoplot(m1 = mafC1, m2 = mafC2, m1Name = 'ALT', m2Name = 'TA', genes = genes, removeNonMutated = FALSE)

# coOncoplot, in default settings, removes genes with 0 mutation across the two phenotypes. coBarplot does not. Manually doing so.
genes_filtered <- genes[genes %in% mafC1@data$Hugo_Symbol | genes %in% mafC2@data$Hugo_Symbol]
coBarplot(m1 = mafC1, m2 = mafC2,
          m1Name = "ALT", m2Name = "TA", genes = genes_filtered)


# drug-drug interactions.
dgiC1 <- drugInteractions(maf = mafC1, fontSize = 0.75)
title("(ALT Osteosarcoma Samples)", cex.main = 0.75)
dev.off()

dgiC2 = drugInteractions(maf = mafC2, fontSize = 0.75)
title("(TA Osteosarcoma Samples)", cex.main = 0.75)
dev.off()


# oncogenic signaling pathways.
pwsC1 = pathways(maf = mafC1, plotType = 'treemap')
title("Oncogenic Signaling Pathways in ALT Osteosarcoma Samples", cex.main = 1.2)
dev.off()

pwsC2 = pathways(maf = mafC2, plotType = 'treemap')
title("Oncogenic Signaling Pathways in TA Osteosarcoma Samples", cex.main = 1.2)
dev.off()

# visualizing the results.
pdf("output/PlotPathways.pdf", width = 6, height = 4)
par(oma = c(1.5, 1.5, 1.5, 1.5))
plotPathways(maf = mafC1, pathlist = pwsC1, showTumorSampleBarcodes = TRUE)
title("Oncogenic Pathways in ALT OS Samples", cex.main = 1.2)
plotPathways(maf = mafC2, pathlist = pwsC2, showTumorSampleBarcodes = TRUE )
title("Oncogenic Pathways in TA OS Samples", cex.main = 1.2)
dev.off()

# heterogeneity in the samples.
pdf("output/altSampleHeterogeneity.pdf", width = 5, height = 6)
for (sample in c1_samples) {
  mafC1.het = inferHeterogeneity(maf = mafC1, tsb = sample, vafCol = 'i_TumorVAF_WU')
  print(mafC1.het$clusterMeans)
  plotClusters(clusters = mafC1.het)
}
dev.off()

pdf("output/taSampleHeterogeneity.pdf", width = 5, height = 6)
for (sample in c2_samples) {
  mafC2.het = inferHeterogeneity(maf = mafC2, tsb = sample, vafCol = 'i_TumorVAF_WU')
  print(mafC2.het$clusterMeans)
  plotClusters(clusters = mafC2.het)
}
dev.off()


### Mutational signatures for ALT and TA phenotypes. hg38 was used.
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

# first ALT samples.
mafC1.tnm = trinucleotideMatrix(maf = mafC1, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

library('NMF')
mafC1.sign = estimateSignatures(mat = mafC1.tnm, nTry = 6)
plotCophenetic(res = mafC1.sign)
mafC1.sig = extractSignatures(mat = mafC1.tnm, n = 5)

#Compare against original 30 signatures 
mafC1.og30.cosm = compareSignatures(nmfRes = mafC1.sig, sig_db = "SBS")

library('pheatmap')
pheatmap::pheatmap(mat = mafC1.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

# finally signatures.
maftools::plotSignatures(nmfRes = mafC1.sig, title_size = 1.2, sig_db = "SBS")


# TA samples.
mafC2.tnm = trinucleotideMatrix(maf = mafC2, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# plotApobecDiff(tnm = mafC2.tnm, maf = mafC2, pVal = 0.2)
# no sample enriched for APOBEC.

mutation_type_cutoff <- 0.1 # keeping mutation that occurs in at least 10% of the samples.
# only way to run NMF is by setting a threshold. Or else, rare mutations interfere with NMF in TA phentoype case.
non_zero_fraction <- colSums(mafC2.tnm$nmf_matrix > 0) / nrow(mafC2.tnm$nmf_matrix)
mafC2.tnm$nmf_matrix <- mafC2.tnm$nmf_matrix[, non_zero_fraction >= mutation_type_cutoff]

mafC2.sign = estimateSignatures(mat = mafC2.tnm, nTry = 5)
plotCophenetic(res = mafC2.sign)
mafC2.sig = extractSignatures(mat = mafC2.tnm, n = 4)

#Compare against original 30 signatures.
mafC2.og30.cosm = compareSignatures(nmfRes = mafC2.sig, sig_db = "SBS")


pheatmap::pheatmap(mat = mafC2.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

# finally signatures.
maftools::plotSignatures(nmfRes = mafC2.sig, title_size = 1.2, sig_db = "SBS")



##############

### Looking at mutation rates of differentially expressed genes.
differentialGenes <- readRDS("Output/geneSignatureCandidates.rds") # from rnaSeq.R.
genes = c(differentialGenes$Gene)
coOncoplot(m1 = mafC1, m2 = mafC2, m1Name = 'C1', m2Name = 'C2', genes = genes, removeNonMutated = FALSE)
coBarplot(m1 = mafC1, m2 = mafC2, m1Name = "C1", m2Name = "C2", genes = genes)

# coOncoplot, in default settings, removes genes with 0 mutation across the two phenotypes. coBarplot does not. Manually doing so.
genes_filtered <- genes[genes %in% mafC1@data$Hugo_Symbol | genes %in% mafC2@data$Hugo_Symbol]
coBarplot(m1 = mafC1, m2 = mafC2,
          m1Name = "ALT", m2Name = "TA", genes = genes_filtered)






