library(tidyverse)
library(dplyr)

# For certain likely high-ploidy samples, PureCN does not produce predictions without manual curation, so these samples are omitted in absolute copy number dataset (only 10 out of 16 samples present) even though we have WGS/WES data for them. 
# we will thus work with relative copy number segment data for this (from GATK).

# Loading segment file
seg <- read_csv("OmicsCNSegmentsProfile.csv")

# formatting the metadata to only include wgs.
metadata <- read_delim("Osteosarcoma_CellLines.txt")
profile <- read_csv("OmicsProfiles.csv") # this file is for mapping Profile ID to ModelID (SampleID/depmapID).
profile <- profile[profile$ModelID %in% metadata$depmap_id, ,drop = FALSE]
profile <-  profile[profile$Datatype == "wgs", ,drop = FALSE]

colnames(seg) <- c("sample", "chromosome", "start", "end", "segmentMean", "numProbes", "Status")


# replacing profile ID with sample ID from profile.
seg <- seg[seg$sample %in% profile$ProfileID, ,drop = FALSE]
seg <- left_join(seg, profile, by = c("sample" = "ProfileID"))
seg$sample <- seg$ModelID 
seg <- seg[, c("sample", "chromosome", "start", "end", "segmentMean", "numProbes", "Status")]

# X chromosomes have unknown status. Removing all sex chromosomes.
seg <- seg[!seg$chromosome %in% c("X", "Y"), , drop = FALSE]


# changing segment mean into gistic format. This data is from GATK and centers around 1. Gistic takes in data centered around 0.
seg <- seg %>%
  mutate(segmentMean = log2(segmentMean))

seg <- seg %>%
  group_by(sample, chromosome)


write.csv(seg, "Output/CNOmicsAll.csv", row.names = FALSE)

# defining c1(ALT) and c2(TA) also -- for gistic files.
segC1 <- seg[seg$sample%in% metadata$depmap_id[metadata$TMMstatus == "ALT"], ,drop = FALSE]
segC2 <- seg[seg$sample%in% metadata$depmap_id[metadata$TMMstatus == "TA"], ,drop = FALSE]

write.csv(segC1, "Output/CNOmicsALT.csv", row.names = FALSE)
write.csv(segC2, "Output/CNOmicsTA.csv", row.names = FALSE)

### making calculations from the Cell Reports Knijnenburg paper.

# 1) number of segments: high segment represents higher genome fragmentation. Normal genome has large contiguous regions with no CNA.
segment_counts <- seg %>%
  group_by(sample) %>%
  summarise(num_segments = n())

segment_counts <- segment_counts %>%
  arrange(desc(num_segments))

segment_counts <- segment_counts %>%
  left_join(metadata, by = c("sample" = "depmap_id"))

# visualization.
ggplot(segment_counts, aes(x = reorder(sample, -num_segments), y = num_segments, fill = TMMstatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Number of segments per sample",
       x = "Sample ID",
       y = "Number of segments") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

model <- lm(num_segments ~ TMMstatus, data = segment_counts)
summary_model <- summary(model)
boxplot(num_segments~TMMstatus, data = segment_counts, ylab = "Number of Segments", xlab = "Class")

# ggplot(segment_counts, aes(x = TMMstatus, y = num_segments)) +
#   geom_violin(trim = FALSE, fill = "lightblue", color = "black") +
#   geom_boxplot(width = 0.1, outlier.shape = NA) +
#   ylab("Number of Segments") + stat_summary(fun.y= mean, geom="point", shape= 21, size= 1) +
#   xlab("Class") +
#   theme_minimal() + annotate("text", x = Inf, y = Inf, 
#                              label = paste("R² =", round(summary_model$r.squared, 3),
#                                            "\np =", round(summary_model$coefficients[2, 4], 4)),
#                              hjust = 1.1, vjust = 1.5, size = 3)

# Adding R-squared value.
legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                  "\np-value =", round(summary_model$coefficients[2, 4], 4)),
       bty = "n", cex = 0.8)

### Difference between number of segments not significant. (higher in ALT phenotype though)


# 2) FGA: Fraction of genome altered. (Altered bases/Total bases)
fga <- seg %>%
  mutate(size = end - start + 1)

fga <- fga %>%
  group_by(sample) %>%
  summarise(
    total_bases = sum(size),
    altered_bases = sum(size[Status %in% c("+", "-")]),
    FGA = altered_bases / total_bases
  )

fga <- fga %>%
  arrange(desc(FGA))

fga <- fga %>%
  left_join(metadata,by = c("sample" = "depmap_id"))


# visualization.
ggplot(fga, aes(x = reorder(sample, -FGA), y = FGA, fill = TMMstatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Fraction of genome altered",
       x = "Sample ID",
       y = "Fraction of genome altered") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

model <- lm(FGA ~ TMMstatus, data = fga)
summary_model <- summary(model)
boxplot(FGA~TMMstatus, data = fga, ylab = "Fraction of genome altered", xlab = "Class")

# Adding R-squared value.
legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                  "\np-value =", round(summary_model$coefficients[2, 4], 4)),
       bty = "n", cex = 0.8)

### Fraction of genome altered: higher in ALT.

#### aneuploidy score: assign each chromosome arm (from hg38) gain or loss.
library(GenomicRanges)

# hg38 cytoband arm data
cytoband <- read.table("cytoBand.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(cytoband) <- c("chromosome", "arm_start", "arm_end", "arm", "gieStain")

# Removing "chr" prefix in chromosome.
cytoband$chromosome <- gsub("^chr", "", cytoband$chromosome)

# creating genomic ranges.
seg_gr <- GRanges(
  seqnames = seg$chromosome,
  ranges = IRanges(start = seg$start, end = seg$end),
  sample = seg$sample,
  status = seg$Status
)

arm_gr <- GRanges(
  seqnames = cytoband$chromosome,
  ranges = IRanges(start = cytoband$arm_start, end = cytoband$arm_end),
  arm = cytoband$arm
)

# overlaps in cytoband and seg.
hits <- findOverlaps(seg_gr, arm_gr)

overlaps <- data.frame(
  sample = mcols(seg_gr)$sample[queryHits(hits)],
  status = mcols(seg_gr)$status[queryHits(hits)],
  chrom = as.character(seqnames(arm_gr)[subjectHits(hits)]),
  arm = mcols(arm_gr)$arm[subjectHits(hits)],
  seg_index = queryHits(hits),
  arm_index = subjectHits(hits),
  stringsAsFactors = FALSE
)

# Adding genomic coordinates for overlap calculations
overlaps <- overlaps %>%
  mutate(
    seg_start = start(seg_gr)[seg_index],
    seg_end = end(seg_gr)[seg_index],
    arm_start = start(arm_gr)[arm_index],
    arm_end = end(arm_gr)[arm_index]
  ) %>%
  mutate(
    overlap_start = pmax(seg_start, arm_start),
    overlap_end = pmin(seg_end, arm_end),
    overlap_bases = pmax(overlap_end - overlap_start + 1, 0)
  ) %>%
  filter(overlap_bases > 0)

# assigning arm level alteration: if 50% or more of the arm is altered by gain or loss, it is counted.
arm_status <- overlaps %>%
  group_by(sample, chrom, arm) %>%
  summarise(
    arm_length = sum(overlap_bases),
    gain_bases = sum(overlap_bases[status == "+"]),
    loss_bases = sum(overlap_bases[status == "-"]),
    gain_frac = gain_bases / arm_length,
    loss_frac = loss_bases / arm_length,
    .groups = "drop"
  ) %>%
  mutate(
    arm_status = case_when(
      gain_frac >= 0.5 & loss_frac < 0.5 ~ "gain",
      loss_frac >= 0.5 & gain_frac < 0.5 ~ "loss",
      gain_frac >= 0.5 & loss_frac >= 0.5 ~ "both",
      TRUE ~ "neutral"
    )
  )

# counting altered arms.
aneuploidy_scores <- arm_status %>%
  filter(arm_status %in% c("gain", "loss", "both")) %>%
  group_by(sample) %>%
  summarise(
    altered_arms = n(),
    gained_arms = sum(arm_status == "gain"),
    lost_arms = sum(arm_status == "loss"),
    both = sum(arm_status == "both"),
    .groups = "drop"
  )

aneuploidy_scores <- aneuploidy_scores %>%
  arrange(desc(altered_arms))

aneuploidy_scores <- aneuploidy_scores %>%
  left_join(metadata, by = c("sample" = "depmap_id"))

# visualization.
ggplot(aneuploidy_scores, aes(x = reorder(sample, -altered_arms), y = altered_arms, fill = TMMstatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Aneuploidy Score",
       x = "Sample ID",
       y = "Aneuploidy Score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

model <- lm(altered_arms ~ TMMstatus, data = aneuploidy_scores)
summary_model <- summary(model)

boxplot(altered_arms~TMMstatus, data = aneuploidy_scores, ylab = "Aneuploidy score", xlab = "Class")

# Adding R-squared value.
legend("topright", legend = paste("R² =", round(summary_model$r.squared, 3), 
                                  "\np-value =", round(summary_model$coefficients[2, 4], 4)),
       bty = "n", cex = 0.8)



###################################################################################

# working with gistic files: finding similarities and differences between c1 and c2 samples in CNA.
gistic_c1 <- read_delim("gistic_data/CNc1/all_lesions.conf_75.txt")
gistic_c1 <- gistic_c1 %>%
  mutate(`Wide Peak Limits` = str_trim(str_replace(`Wide Peak Limits`, "\\s*\\([^\\)]+\\)", ""))) # remove the values inside brackets for comparison of c1 and c2.

gistic_c2 <- read_delim("gistic_data/CNc2/all_lesions.conf_75.txt")
gistic_c2 <- gistic_c2 %>%
  mutate(`Wide Peak Limits` = str_trim(str_replace(`Wide Peak Limits`, "\\s*\\([^\\)]+\\)", ""))) # remove the values inside brackets for comparison of c1 and c2.

gistic_c1 <- gistic_c1[1:72, ,drop = FALSE] # values got repeated so just taking unrepeated ones (hardcoding based on our data).
gistic_c2 <- gistic_c2[1:52, ,drop = FALSE]

# Creating composite keys for both datasets to find common ranges.
get_intersect_key <- function(df) {
  first_word <- word(df$`Unique Name`, 1)
  df$key <- paste(df$Descriptor, first_word, sep = "_")
  return(df)
}

# Applying to both dataframes
gistic_c1 <- get_intersect_key(gistic_c1)
gistic_c2 <- get_intersect_key(gistic_c2)

# intersecting keys
common_keys <- intersect(gistic_c1$key, gistic_c2$key)

# unqiue CNAs for each phenotype.
cytoband_c1 <- gistic_c1[!(gistic_c1$key %in% common_keys), ]
cytoband_c1 <- get_intersect_key(cytoband_c1)

cytoband_c2 <- gistic_c2[!(gistic_c2$key %in% common_keys), ]
cytoband_c2 <- get_intersect_key(cytoband_c2)

