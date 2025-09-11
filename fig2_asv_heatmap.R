library(dada2)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(phyloseq)
library(vegan)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(ggrepel)

source("utils.R")

working_dir <- "/data/local/jy1008/Allegretti"
meta_data <- read.table(file.path(working_dir, "metadata",
                "meta_data_subject.csv"), header = TRUE,
                sep = ",", stringsAsFactors = FALSE)

time_data <- read.table(file.path(working_dir, "metadata",
                "meta_time_series_mayo_scores.csv"), header = TRUE,
                sep = ",", stringsAsFactors = FALSE)
time_data <- time_data %>%
    dplyr::rename(sample_id = subject_id)
# careful: rename() gets masked by S4Vectors

meta_sub <- time_data[, c("sample_id", "total_partial_mayo")]
meta_sub <- meta_sub %>%
    extract(sample_id, into = c("subject_id", "timepoint"),
            regex = "(CP\\d+)([A-Za-z]+\\d+)", remove = FALSE)
meta_sub$subject_id <- str_replace(meta_sub$subject_id,
                                "^(CP\\d{3})(\\d{2})$", "\\1-\\2")
meta_sub <- merge(meta_sub, meta_data[, c("subject_id", "treatment")],
                  by = "subject_id", all.x = TRUE)

# convert all time points to weeks
meta_sub$time_numeric <- as.numeric(dplyr::recode(
  meta_sub$timepoint,
  "Day1" = "0",
  "Wk1" = "1",
  "Wk4" = "4",
  "Wk8" = "8",
  "Wk12" = "12",
  "Wk16" = "16",
  "Wk32" = "32",
  .default = NA_character_
))

# data and figure output dirs from dada2
dada2_data_out <- "test_data_output_v2"
dada2_figure_out <- "test_figure_output_v2"

seqtab.nochim <- readRDS(file.path(working_dir,
                    dada2_data_out, "seqtab_final.rds"))

# filter for time series samples
seqtab.filt <- seqtab.nochim[rownames(seqtab.nochim) %in%
                    time_data$sample_id, ]
# check what was removed
setdiff(rownames(seqtab.nochim), rownames(seqtab.filt))

# assign ASVs to taxa (family or genus)
taxa <- read.table(
  file.path(working_dir, dada2_data_out, "taxa.csv"),
  sep = " ",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "\"",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE
)
# colnames(taxa)[1] <- "ASV"

#
# Filtering ASVs
#


# filter out ASVs which are not assigned to any taxa
taxa_filtered <- taxa[!taxa$label %in% c("NA_NA_NA", NA), ]
seqtab.filt2 <- seqtab.filt[, colnames(seqtab.filt) %in%
                    rownames(taxa_filtered)]
# get rid of all-zero ASVs
seqtab.filt2 <- seqtab.filt2[, colSums(seqtab.filt2) > 0]
taxa_filtered <- taxa_filtered[colnames(seqtab.filt2), , drop = FALSE]

# plot elbow plots (from utils.R)
plot_qc_elbow(seqtab_filtered)


# # total counts per ASV (across all samples)
# asv_totals <- colSums(seqtab.filt2)
# # Filter: Keep ASVs with at least 10 reads
# seqtab.filt3 <- seqtab.filt2[, asv_totals >= 10]

# asv_prevalence <- colSums(seqtab.filt3 > 0)
# # Keep ASVs present in at least a % of samples
# seqtab.filt4 <- seqtab.filt3[,
#                     asv_prevalence >= 0.2 * nrow(seqtab.filt3)]

# Step 1: relative abundances per sample
relabund <- sweep(seqtab.filt2, 1, rowSums(seqtab.filt2), FUN = "/")

# Step 2: check threshold (0.0015 = 0.15%)
pass_threshold <- colSums(relabund >= 0.001) >= 2

# Step 3: filter ASVs
seqtab.filt3 <- seqtab.filt2[, pass_threshold]

# Original ASV sequences as names
seqtab_to_use <- seqtab.filt3
asv_seqs <- colnames(seqtab_to_use)

taxa_filtered <- taxa_filtered[colnames(seqtab.filt3), , drop = FALSE]

# Assign short IDs: ASV1, ASV2, ...
asv_ids <- paste0("ASV", seq_along(asv_seqs))
names(asv_ids) <- asv_seqs  # Map from sequence to ID
colnames(seqtab_to_use) <- asv_ids[colnames(seqtab_to_use)]
# > colnames(seqtab_to_use)[1:3]
# TACGTAGGGGGCAAGCGTTATCC...  "ASV1"
# AACGTAGGTCACAAGCGTTGTCC...  "ASV2"
# AACGTAGGTCACAAGCGTTGTCC...  "ASV3"

# structure(
#   c("ASV1", "ASV2", "ASV3"),
#   names = c("TACGTAGGG...", "AACGTAGGT...", "AACGTAGGT...")
# )


asv_mapping <- data.frame(
  ASV_ID = asv_ids,
  Sequence = names(asv_ids),
  Taxonomy = taxa[names(asv_ids), "label"],
  stringsAsFactors = FALSE
)
get_asv_id <- function(sequence) {
  matched_row <- asv_mapping[asv_mapping$Sequence == sequence, ]
  if (nrow(matched_row) == 0) {
    return(NA)
  }
  return(paste(matched_row$Taxonomy, matched_row$ASV_ID, sep=" | "))
}
write.csv(asv_mapping,
          file.path(working_dir, dada2_data_out, "fig2_heatmap_filt_asv_mapping.csv"),
          row.names = FALSE)






# Keep top N most variable ASVs
# asv_variance <- apply(seqtab.filt4, 2, var)
# top_var_asvs <- names(sort(asv_variance, decreasing = TRUE))[1:50]
# seqtab.filt5 <- seqtab.filt4[, top_var_asvs]

# collapse into taxa at level of family
# Transpose seqtab.filt so ASVs are rows (matching taxa_df)
seqtab_df <- as.data.frame(t(seqtab_to_use))
seqtab_df$ASV <- rownames(seqtab_df)

# Join with family info
# convert rownames of taxa_filtered to ASV[Num]
temp_asvs <- get_asv_id(rownames(taxa_filtered))
temp_asvs <- sapply(strsplit(temp_asvs, " | "),
                    function(z) trimws(tail(z, 1)))
# rownames(taxa_filtered) <- temp_asvs
taxa_filtered$ASV <- temp_asvs

seqtab_taxa <- seqtab_df %>%
  left_join(taxa_filtered[, c("ASV", "Family")], by = "ASV")

# Sum counts per family
family_counts <- seqtab_taxa %>%
  group_by(Family) %>%
  summarise(across(-ASV, sum))  # sum across all ASVs in that family

# Compute total abundance across all samples
family_counts <- family_counts %>%
  mutate(Total = rowSums(across(-Family))) %>%
  arrange(desc(Total))  # sort descending by total abundance

# Select top N families
top_n <- 15
top_families <- family_counts$Family[1:top_n]

# --- Fold remaining families into "Other" ---
family_counts_collapsed <- family_counts %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other"))


# --- Sum counts per family, including "Other" ---
family_counts_collapsed_summed <- family_counts_collapsed %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), sum)) %>%  # sum all samples for duplicates (i.e., "Other")
  ungroup()

# --- Convert to matrix: rows = families, columns = samples ---
family_mat <- as.matrix(family_counts_collapsed_summed %>%
                           select(-Family, -Total))   # remove both Family and Total
rownames(family_mat) <- family_counts_collapsed_summed$Family

# --- Plot heatmap ---
ht <- plot_asv_heatmap2(
  t(family_mat),
  meta_sub,
  min_abundance = 0,
  min_prevalence = 0,
  rel_abund = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE
)
pdf(file.path(working_dir, dada2_figure_out, "fig2_family_heatmap.pdf"),
    width = 12, height = 8)
draw(ht, heatmap_legend_side = "bottom")
dev.off()



# stacked bar plots of family composition
# --- Reshape for plotting ---
df_long <- family_counts_collapsed %>%
  select(-Total) %>%
  pivot_longer(-Family, names_to = "Sample", values_to = "Abundance") %>%
  # Join metadata for subject info
  left_join(meta_sub %>% select(sample_id, subject_id), by = c("Sample" = "sample_id")) %>%
  # Compute relative abundance per sample
  group_by(Sample) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  # --- Extract last numeric value from sample names ---
  mutate(
    last_num = sapply(str_extract_all(Sample, "\\d+"), function(x) as.numeric(tail(x, 1))),
    time_num = ifelse(str_detect(Sample, "(?i)^Day\\d+"),
                      last_num / 7,
                      last_num),
    Time = str_replace(Sample, "^CP\\d+\\d+", "")
  ) %>%
  # --- Reorder samples within each subject by time_num ---
  group_by(subject_id) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(time_num)]))) %>%
  mutate(Time = factor(Time, levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16"))) %>%
  ungroup()
# --- Extract time point ---
# df_long <- df_long %>%
#   mutate(Time = str_replace(Sample, "^CP\\d+\\d+", ""))

# Compute total abundance per Family
family_totals <- df_long %>%
  group_by(Family) %>%
  summarise(total_abund = sum(Abundance), .groups = "drop")
# Order families by descending total abundance, exclude "Other"
ordered_families <- family_totals %>%
  filter(Family != "Other") %>%
  arrange(desc(total_abund)) %>%
  pull(Family)
# Add "Other" to the end
ordered_families <- c(ordered_families, "Other")
# Convert Family into factor with this order
df_long$Family <- factor(df_long$Family, levels = ordered_families)

# --- Plot faceted stacked barplots ---
sb <- ggplot(df_long, aes(x = Time, y = RelAbund, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Time", y = "Relative abundance", fill = "Family") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  facet_wrap(~ subject_id, ncol = 4, scales = "free_x") +
  scale_fill_viridis_d(option = "turbo")

pdf(file.path(working_dir, dada2_figure_out, "fig2_family_top15_stacked_bars.pdf"),
    width = 12, height = 8)
print(sb)
dev.off()







# ---------- Alpha Diversity Analysis ----------
# Alpha diversity using filtered ASVs (seqtab.filt3)
# seqtab.filt2 has only removal of unassigned ASVs (no read count or prevalence filtering)
# seqtab_to_use <- seqtab.filt2
seqtab_to_use <- seqtab.filt3
sample_order <- rownames(seqtab_to_use)
meta_all <- meta_sub[match(sample_order, meta_sub$sample_id), ]
rownames(meta_all) <- meta_all$sample_id
missing_ids <- setdiff(sample_order, meta_all$sample_id)
if (length(missing_ids) > 0) {
  warning("Some sample IDs not found in df: ", paste(missing_ids, collapse = ", "))
}
print("checking meta row order alignment")
all(rownames(seqtab_to_use) == meta_all$sample_id)

taxa <- read.csv(file.path(working_dir, dada2_data_out, "taxa.csv"), sep=' ')
taxa <- as.matrix(taxa)
# Create the taxa table for only ASVs in seqtab.filt5
# asv_ids5 <- colnames(seqtab.filt4) # AV1, AV2
# Step 2: Map ASV IDs back to ASV sequences
# asv_seqs5 <- names(asv_ids)[asv_ids %in% asv_ids5]
# Step 3: Subset the taxonomy table using the ASV sequences
# taxa_filtered <- taxa[asv_seqs5, , drop = FALSE]
# rownames(taxa_filtered) <- asv_ids[asv_seqs5]
taxa_filtered <- taxa[colnames(seqtab_to_use), , drop = FALSE]
print("checking taxa row order alignment")
all(colnames(seqtab_to_use) == rownames(taxa_filtered))
ps <- phyloseq(otu_table(seqtab_to_use, taxa_are_rows=FALSE),
               sample_data(meta_all),
               tax_table(taxa_filtered))

# NOTE: Chao1 returns stderr 0 for many samples because there are very few singletons
# ✅ Recommendations
# If you’ve filtered rare ASVs, use Shannon or Observed Richness, not Chao1.
# If your focus is on community structure or evenness, go with Shannon.
# If you care most about number of taxa, use Observed Richness.
# alpha_div <- estimate_richness(ps, measures = "Chao1")
alpha_div <- estimate_richness(ps, measures = c("Shannon"))
alpha_div$sample_id <- rownames(alpha_div)
alpha_div_merged <- merge(alpha_div, meta_all, by = "sample_id")
alpha_div_merged$timepoint <- factor(alpha_div_merged$timepoint, levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16", "Wk32"))
alpha_div_merged$treatment <- factor(alpha_div_merged$treatment, levels = c("placebo", "active"))

ad <- ggplot(alpha_div_merged, aes(x = time_numeric, y = Shannon, color = treatment, group = subject_id)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ subject_id, ncol = 4, scales = "free_y") +  # free_y gives separate y-axis per subject
  scale_x_continuous(name = "Time (weeks)") +
  scale_y_continuous(name = "Shannon Diversity") +
  scale_color_manual(values = c("placebo" = "#1f77b4", "active" = "#ff7f0e")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
pdf(file.path(working_dir, dada2_figure_out, "fig3_alpha_diversity_trend_per_subject.pdf"),
    width = 12, height = 8)
print(ad)
dev.off()


# Wilcoxon test for treatment differences at each timepoint
# Compute delta from baseline
library(rstatix)
alpha_diff <- alpha_div_merged %>%
  group_by(subject_id) %>%
  mutate(
    baseline = Shannon[time_numeric == 0][1],   # take the baseline Shannon for this subject
    delta_Shannon = Shannon - baseline
  ) %>%
  ungroup()

wilcox_results <- alpha_diff %>%
  filter(time_numeric != 0) %>%  # exclude baseline itself
  group_by(time_numeric) %>%
  wilcox_test(delta_Shannon ~ treatment) %>%  # rank-sum test
  adjust_pvalue(method = "BH") %>%           # optional multiple testing correction
  add_significance()

ggplot(alpha_diff %>% filter(time_numeric != 0),
       aes(x = time_numeric, y = delta_Shannon, color = treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = "Time (weeks)", y = "Change in Shannon Diversity", color = "Treatment") +
  theme_minimal()

# ----- pairwise alpha-diversity plot -----
library(cowplot)

# Ensure time points are numeric and ordered
time_points <- sort(unique(alpha_div_merged$time_numeric))

# Store plots in a list
plots <- list()

# Loop over successive pairs
for (i in 1:(length(time_points)-1)) {
  
  t1 <- time_points[i]
  t2 <- time_points[i+1]
  
  # Subset to just the two time points
  df_pair <- alpha_div_merged %>%
    filter(time_numeric %in% c(t1, t2))
  
  # Compute treatment averages
  avg_df <- df_pair %>%
    group_by(treatment, time_numeric) %>%
    summarize(mean_Shannon = mean(Shannon, na.rm = TRUE), .groups = "drop")
  
  time_labels <- unique(df_pair$timepoint[df_pair$time_numeric %in% c(t1, t2)])

  # Create plot
  p <- ggplot(df_pair, aes(x = time_numeric, y = Shannon, group = subject_id, color = treatment)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    # Add labels for each point
    # geom_text_repel(aes(label = subject_id),
    #                 size = 3,
    #                 max.overlaps = Inf,  # ensures all labels try to be drawn
    #                 show.legend = FALSE) +
    # Mean lines
    geom_line(data = avg_df, aes(x = time_numeric, y = mean_Shannon, 
                                 group = treatment, color = treatment),
              size = 1.2) +
    scale_x_continuous(breaks = c(t1, t2), labels = time_labels) +
    labs(title = paste0("Time pair: ", df_pair$timepoint[df_pair$time_numeric==t1][1],
                        " -> ", df_pair$timepoint[df_pair$time_numeric==t2][1]),
         x = "Time", y = "Shannon diversity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store plot
  plots[[i]] <- p
  
  # Optional: print immediately
  print(p)
}



# Combine plots in a grid
combined_plot <- plot_grid(plotlist = plots, nrow=1)
ggsave(file.path(working_dir, dada2_figure_out, "alpha_diversity_Shannon_pairwise_wilcox_vis.pdf"),
        plot = combined_plot, width = 16, height = 4, dpi = 300)

# ----- -----

# V1 plot
p <- ggplot(alpha_div_merged, aes(x = timepoint, y = Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5) +
  labs(
    title = "Alpha Diversity Across Time Points by Treatment",
    x = "Timepoint",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
ggsave(file.path(working_dir, dada2_figure_out, "alpha_diversity_Shannon_boxplot_vs_treatment_v1.pdf"),
        plot = p, width = 8, height = 6, dpi = 300)


# 9/4/2025
# very simple test: is the alpha diversity different between treatment groups at each time point?
result_alpha_per_time <- alpha_div_merged %>%
  group_by(timepoint) %>%
  summarise(
    test = list(
      broom::tidy(
        wilcox.test(Shannon ~ treatment, data = cur_data())
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest(test)

result_alpha_per_time

# # A tibble: 6 × 5
#   timepoint statistic p.value method                       alternative
#   <fct>         <dbl>   <dbl> <chr>                        <chr>      
# 1 Day1             26   0.867 Wilcoxon rank sum exact test two.sided  
# 2 Wk1              42   0.328 Wilcoxon rank sum exact test two.sided  
# 3 Wk4              12   1     Wilcoxon rank sum exact test two.sided  
# 4 Wk8              14   1     Wilcoxon rank sum exact test two.sided  
# 5 Wk12             16   0.788 Wilcoxon rank sum exact test two.sided  
# 6 Wk16             17   0.183 Wilcoxon rank sum exact test two.sided
# the statistic is the sum of the ranks in the first group (active)

# > table(alpha_div_merged$timepoint)
# Day1  Wk1  Wk4  Wk8 Wk12 Wk16 Wk32 
#   15   16   10   11   11   10    0


# V2 plot (with lines connecting points per subject)
# Ensure timepoint is a factor in the correct order
alpha_div_merged$timepoint <- factor(alpha_div_merged$timepoint, levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16", "Wk32"))
# Assign numeric values for x-axis positions
alpha_div_merged$timepoint_num <- as.numeric(alpha_div_merged$timepoint)
# Apply dodge offset (e.g., -0.2 for placebo, +0.2 for active)
alpha_div_merged$x_dodge <- alpha_div_merged$timepoint_num +
  ifelse(alpha_div_merged$treatment == "active", 0.2, -0.2)
p <- ggplot(alpha_div_merged) +
  # Boxplots
  geom_boxplot(
    aes(x = timepoint, y = Shannon, fill = treatment),
    outlier.shape = NA,
    alpha = 0.7,
    width=0.4,
    position = position_dodge(width = 0.4)
  ) +
  # Lines connecting points per subject (use x_dodge!)
  geom_line(
    aes(x = x_dodge, y = Shannon, group = subject_id, color = treatment),
    alpha = 0.4
  ) +
  # Points at dodged positions
  geom_point(
    aes(x = x_dodge, y = Shannon, fill = treatment),
    shape = 21, color = "black", size = 2
  ) +
  scale_x_discrete(
    name = "Timepoint",
    breaks = 1:length(levels(alpha_div_merged$timepoint)),
    labels = levels(alpha_div_merged$timepoint)
  ) +
  labs(
    title = "Alpha Diversity Across Time Points by Treatment",
    y = "Shannon Diversity Index"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
ggsave(file.path(working_dir, dada2_figure_out, "alpha_diversity_Shannon_boxplot_vs_treatment_v2.pdf"),
        plot = p, width = 8, height = 6, dpi = 300)








  # Wilcoxon test
  # This creates a list of comparisons like: list(c("active", "placebo"))
treatment_comparisons <- alpha_div_merged %>%
  filter(!is.na(treatment)) %>%
  distinct(timepoint) %>%
  pull(timepoint) %>%
  lapply(function(tp) {
    subset_df <- filter(alpha_div_merged, timepoint == tp)
    if (length(unique(subset_df$treatment)) == 2) {
      list(tp = tp, p = wilcox.test(Shannon ~ treatment, data = subset_df)$p.value)
    } else {
      NULL
    }
  }) %>%
  compact()
do.call(rbind, lapply(treatment_comparisons, as.data.frame))
# Shannon
#     tp         p
# 1 Day1 0.9550894
# 2  Wk1 0.3282051
# 3 Wk12 0.7878788
# 4 Wk16 0.1166667
# 5  Wk8 0.9272727
# 6  Wk4 0.9142857


# Find pairs of successive timepoints
timepoints <- sort(unique(alpha_div_merged$timepoint))
timepoint_pairs <- Map(c, timepoints[-length(timepoints)], timepoints[-1])

timepoint_comparisons <- lapply(unique(alpha_div_merged$treatment), function(trt) {
  lapply(timepoint_pairs, function(tp_pair) {
    df <- alpha_div_merged %>% filter(treatment == trt, timepoint %in% tp_pair)
    if (length(unique(df$timepoint)) == 2) {
      p <- wilcox.test(Shannon ~ timepoint, data = df)$p.value
      list(trt = trt, tp1 = tp_pair[1], tp2 = tp_pair[2], p = p)
    } else {
      NULL
    }
  }) %>% compact()
}) %>% unlist(recursive = FALSE)
do.call(rbind, lapply(timepoint_comparisons, as.data.frame))
# Shannon
#        trt  tp1  tp2         p
# 1   active Day1  Wk1 1.0000000
# 2   active  Wk1  Wk4 0.3449883
# 3   active  Wk4  Wk8 0.7307692
# 4   active  Wk8 Wk12 0.9015152
# 5   active Wk12 Wk16 0.9015152
# 6  placebo Day1  Wk1 0.4418026
# 7  placebo  Wk1  Wk4 0.4606061
# 8  placebo  Wk4  Wk8 0.6857143
# 9  placebo  Wk8 Wk12 0.6857143
# 10 placebo Wk12 Wk16 0.1142857



# beta diversity
# Bray-Curtis distance
bc_dist <- phyloseq::distance(ps, method = "bray")
set.seed(101)
adonis2(bc_dist ~ subject_id, data = data.frame(sample_data(ps)))
# adonis2(formula = bc_dist ~ subject_id, data = data.frame(sample_data(ps)))
#          Df SumOfSqs      R2      F Pr(>F)    
# Model    15  13.8480 0.75064 11.439  0.001 ***
# Residual 57   4.6002 0.24936                  
# Total    72  18.4482 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
set.seed(101)
adonis2(bc_dist ~ treatment, data = data.frame(sample_data(ps)))
# adonis2(formula = bc_dist ~ treatment, data = data.frame(sample_data(ps)))
#          Df SumOfSqs      R2      F Pr(>F)    
# Model     1   0.8096 0.04389 3.2589  0.001 ***
# Residual 71  17.6386 0.95611                  
# Total    72  18.4482 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
meta_beta_div <- data.frame(sample_data(ps))
set.seed(101)
adonis2(bc_dist ~ treatment * time_numeric,
          data = meta_beta_div, strata = meta_beta_div$subject_id,
          by = "term")
# By using strata = subject_id, you ensure that:
# Permutations are restricted within each subject,
# so that sample identities are preserved during testing.
# This avoids artificially inflating significance
# by comparing samples across individuals (e.g. in
# longitudinal studies where samples are not independent).

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata
# Permutation: free
# Number of permutations: 999
#
# adonis2(formula = bc_dist ~ treatment * time_numeric, data = meta_beta_div, by = "term", strata = meta_beta_div$subject_id)
#                        Df SumOfSqs      R2      F Pr(>F)   
# treatment               1   0.8096 0.04389 3.2563  0.016 * 
# time_numeric            1   0.2897 0.01570 1.1652  0.005 **
# treatment:time_numeric  1   0.1935 0.01049 0.7783  0.524   
# Residual               69  17.1553 0.92992                 
# Total                  72  18.4482 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# PCoA is like UMAP but linear (and thus interpretable - standard for beta diversity)
# Step 1: Ordinate
ordination <- ordinate(ps, method = "PCoA", distance = "bray")
# Step 2: Extract ordination as dataframe, including metadata
ord_df <- plot_ordination(ps, ordination, justDF = TRUE)
ord_df$timepoint_plot <- factor(ord_df$timepoint,
                          levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16", "Wk32"))
ord_df$timepoint_plot <- as.numeric(ord_df$timepoint_plot)
ord_df$subject_id <- as.character(ord_df$subject_id)
label_df <- ord_df %>%
  group_by(subject_id) %>%
  arrange(timepoint_plot) %>%
  filter(timepoint_plot == min(timepoint_plot))  # take first time point per subject
# Step 3: Plot
p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = treatment, group = subject_id)) +
  geom_point(aes(size = timepoint_plot), alpha = 0.8) +
  geom_path(alpha = 0.4) +  # connects samples from the same subject
  geom_text(
    data = label_df,
    aes(label = subject_id),
    size = 3, vjust = -1
  ) +
  scale_size_continuous(range = c(2, 8)) +  # adjust size scaling
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances",
    subtitle = "Larger points indicate later timepoints",
    x = "PCoA Axis 1", y = "PCoA Axis 2",
    size = "Time", color = "Treatment"
  )
ggsave(file.path(working_dir, dada2_figure_out, "beta_diversity_timecourse.pdf"), plot=p)
p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = treatment, group = subject_id)) +
  geom_point(alpha = 0.8) +
  geom_path(alpha = 0.4) +  # connects samples from the same subject
  geom_text(
    data = label_df,
    aes(label = subject_id),
    size = 3, vjust = -1
  ) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances",
    subtitle = "Larger points indicate later timepoints",
    x = "PCoA Axis 1", y = "PCoA Axis 2",
    size = "Time", color = "Treatment"
  )
ggsave(file.path(working_dir, dada2_figure_out, "beta_diversity_timecourse_nosize.pdf"), plot=p)

# which taxa contribute most to Axis.2?
# Extract taxa abundance matrix (samples x taxa)
ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))
# otu_mat <- as(otu_table(ps), "matrix")
otu_mat <- as(otu_table(ps_relabund), "matrix")
# Ensure samples match ord_df rows (same order)
all(rownames(otu_mat) == ord_df$SampleID)  # should be TRUE
# Extract Axis.2 sample scores
axis2 <- ord_df$Axis.2
taxa_cors <- apply(otu_mat, 2, function(x) cor(x, axis2, method = "spearman"))
# Rank taxa by absolute correlation strength
taxa_ranked <- sort(abs(taxa_cors), decreasing = TRUE)
taxa_pvals <- apply(otu_mat, 2, function(x) cor.test(x, axis2, method = "spearman")$p.value)
# Combine correlation and p-values
taxa_stats <- data.frame(
  taxon = names(taxa_cors),
  correlation = taxa_cors,
  p_value = taxa_pvals
)
# Adjust p-values for multiple testing
taxa_stats$padj <- p.adjust(taxa_stats$p_value, method = "BH")
taxa_stats <- taxa_stats %>%
  rowwise() %>%
  mutate(asv_id = get_asv_id(taxon)) %>%
  ungroup()
write.csv(taxa_stats,
          file.path(working_dir, dada2_data_out, "fig3_beta_div_taxa_stats.csv"),
          row.names = FALSE)

# Filter significant taxa
signif_taxa <- subset(taxa_stats, padj < 0.05)
# Sort by strength of correlation
signif_taxa <- signif_taxa[order(abs(signif_taxa$correlation), decreasing = TRUE), ]
pdf(file.path(working_dir, dada2_figure_out, "beta_div_asvs_relabund.pdf"), width = 5.5, height = 5)  # open multi-page PDF
for(sig_asv in as.character(signif_taxa$taxon[1:20])) {
  # asv_counts <- otu_table(ps)[, sig_asv]  # replace with your ASV ID
  # asv_counts_vec <- as.numeric(asv_counts)  # vector of counts
  # names(asv_counts_vec) <- sample_names(ps)
  asv_relabund <- otu_table(ps_relabund)[, sig_asv]
  asv_vec <- as.numeric(asv_relabund)
  names(asv_vec) <- sample_names(ps_relabund)

  # Add to ordination data frame
  # ord_df$ASV_count <- asv_counts_vec[ord_df$sample_id]
  ord_df$ASV_relabund <- asv_vec[ord_df$sample_id]
  sig_asv_id <- get_asv_id(sig_asv)
  print(sig_asv_id)
  if (is.na(sig_asv_id)) {
    # NOTE: sig_asv was built from seqtab.filt4, so it has a nz in >= 20% samples filter
    warning("ASV ID not found for sequence: ", sig_asv)
    next
  }
  # Plot with coloring by ASV abundance
  # p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = ASV_count, group = subject_id)) +
  p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = ASV_relabund, group = subject_id)) +
    geom_point(size = 3) +
    geom_path(alpha = 0.4) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal() +
  #   labs(color = paste("read count"),
  labs(color = paste("Relative Abundance"),
        title = paste("PCoA: ", sig_asv_id))
  print(p)
}
dev.off()

#
# ---------- Differential abundance analysis ----------
#


# Alternate: run DESeq2 pairwise
library(DESeq2)
library(apeglm)

# NOTE: for ASV mapping to work, make sure phyloseq object
# is generated from the same as asv_mapping (seqtab.filt3)

timepoints <- c("Wk1", "Wk4", "Wk8", "Wk12", "Wk16")
dds_list <- list()   # store full DESeqDataSet for each time point

for(tp in timepoints) {
  # Subset samples for current time point
  # samples_to_keep <- sample_data(ps)$timepoint %in% c("Day1", tp)
  samples_to_keep <- sample_data(ps)$timepoint %in% c(tp)
  ps_sub <- prune_samples(samples_to_keep, ps)

  # Convert variables to factors
  sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, levels = c("placebo", "active"))
  # sample_data(ps_sub)$time <- factor(sample_data(ps_sub)$timepoint, levels = c("Day1", tp))

  # Convert to DESeq2 object
  dds_sub <- phyloseq_to_deseq2(ps_sub, ~ treatment)

  # Run DESeq
  dds_sub <- DESeq(dds_sub)

  # Store full results object
  dds_list[[tp]] <- dds_sub
}


library(ggrepel)
output_dir <- "deseq2_results_v3"
dir.create(output_dir, showWarnings = FALSE)

for(tp in names(dds_list)) {
  dds_sub <- dds_list[[tp]]

  # Extract results
  res_df <- as.data.frame(results(dds_sub))
  res_df$ASV <- rownames(res_df)

  # Write CSV
  write.csv(res_df, file.path(output_dir, paste0("DESeq2_", tp, "_full_results.csv")), row.names = FALSE)

  # Prepare volcano plot data
  res_df <- res_df %>%
    mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      sig = padj < 0.05
    )

  # Count significant ASVs and total ASVs
  n_sig <- sum(res_df$sig, na.rm = TRUE)
  n_total <- nrow(res_df)
  label_text <- paste0("Significant: ", n_sig, "\nTotal: ", n_total)


  res_df <- res_df %>%
    mutate(ASV_label = sapply(ASV, get_asv_id))

  # Volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj_plot))) +
    geom_point(aes(color = sig), alpha = 0.7) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      data = subset(res_df, sig),
      aes(label = ASV_label),
      size = 3,
      max.overlaps = 20
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = label_text,
      hjust = 1.1, vjust = 1.1,
      size = 4
    ) +
    labs(
      title = paste("Volcano plot:", tp),
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    ) +
    theme_minimal()

  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0("Volcano_", tp, ".png")),
    plot = p,
    width = 7,
    height = 6
  )
}



sig_asvs_list <- lapply(dds_list, function(res) {
  df <- as.data.frame(results(res))
  df$ASV <- rownames(df)  # keep the ASV IDs
  df %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    pull(ASV)
})
all_sig_asvs <- unique(unlist(sig_asvs_list))

# Flatten list to a vector of all significant ASVs across comparisons
all_sig_asvs_flat <- unlist(sig_asvs_list)
# Count occurrences per ASV
asv_counts <- table(all_sig_asvs_flat)
# Keep ASVs that appear in at least 2 comparisons
all_sig_asvs_2plus <- names(asv_counts[asv_counts >= 2])

# sapply(all_sig_asvs, get_asv_id, USE.NAMES = FALSE)
seqtab.filt3.sig <- seqtab.filt3[, all_sig_asvs_2plus]
colnames(seqtab.filt3.sig) <-
  sapply(colnames(seqtab.filt3.sig), get_asv_id, USE.NAMES = FALSE)
ht <- plot_asv_heatmap2(seqtab.filt3.sig,
                    meta_sub,
                    min_abundance = 0,
                    min_prevalence = 0,
                    rel_abund = TRUE,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE)
pdf(file.path(working_dir, dada2_figure_out, "fig5_deseq2_pairedtime_sig_asvs_2plus.pdf"),
    width = 12, height = 8)
draw(ht, heatmap_legend_side = "bottom")
dev.off()




# visualize across time points
res_all <- c()
for(tp in names(dds_list)) {
  dds_sub <- dds_list[[tp]]

  # Extract results
  res_df <- as.data.frame(results(dds_sub))
  res_df$comparison <- tp
  res_df$ASV <- rownames(res_df)  # keep the ASV IDs
  res_df$ASV_label <- sapply(res_df$ASV, get_asv_id, USE.NAMES = FALSE)
  res_all <- rbind(res_all, res_df)
}
rownames(res_all) <- NULL
res_all$ASV <- NULL

write.csv(res_all, 
          "deseq2_res_all.csv", row.names = FALSE)


res_all_sig <- res_all[!is.na(res_all$padj) & res_all$padj < 0.05, ]
asv_label_set <- unique(res_all_sig$ASV_label)
res_all_disp <- res_all[res_all$ASV_label %in% asv_label_set, ]
res_all_disp$padj[is.na(res_all_disp$padj)] <- 1
res_all_disp$nlog10padj <- -log10(res_all_disp$padj)

# First, add columns for color (sign) and border (significance)
res_all_disp <- res_all_disp %>%
  mutate(
    sign = ifelse(log2FoldChange > 0, "up", "down"),
    sig_border = ifelse(padj < 0.05, "significant", "ns")
  )
res_all_disp$comparison <- factor(res_all_disp$comparison, levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16"))

asv_mat <- res_all_disp %>%
  select(ASV_label, comparison, log2FoldChange) %>%
  pivot_wider(
    names_from = comparison,
    values_from = log2FoldChange,
    values_fill = 0
  )
# Compute distance and clustering using the matrix without row names
asv_values <- asv_mat %>% select(-ASV_label)
dist_mat <- dist(as.matrix(asv_values))
hc <- hclust(dist_mat, method = "complete")
asv_order <- asv_mat$ASV_label[hc$order]

res_all_disp$ASV_label <- factor(res_all_disp$ASV_label, levels = asv_order)

# Dotplot
p <- ggplot(res_all_disp, aes(
  x = comparison,
  y = ASV_label,
  size = nlog10padj,   # size by magnitude of fold change
  fill = log2FoldChange,
  stroke = ifelse(sig_border == "significant", 1.5, 0)  # border thickness
)) +
  geom_point(shape = 21) +   # shape 21 allows border color & stroke
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,   # 0 fold change = neutral
    name = "log2FC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    x = "Time point",
    y = "ASV",
    color = "Direction",
    size = "-log10 padj"
  )
ggsave("dada2_deseq2_dotplot.pdf", plot = p, width = 6, height = 10)





# investigating extremely high fold changes
# Prevotellaceae_Prevotellaceae UCG-001_NA | ASV264
asv_str <- "TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGCAGGCCGCCGTGCAAGCGTGCCGTGAAAAGCAGCGGCCCAACCGCTGCCCTGCGGCGCGAACTGCTTGGCTTGAGTGCGCCGGAAGCGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCCCGCTGTGGCGCCACTGACGCTGAGGCTCGAAGGTGCGGGTATCGAACAGG"
temp <- seqtab.filt3[, asv_str]
print(temp[temp > 0])


# the significant asvs are significant due to extreme sparsity
# either:
# --- run DESeq2 but more aggressively filtering out sparse ASVs
# --- run a Fisher exact test (treatment vs no treatment) x (time 0 vs time t) (but this might be too simple, i.e. there's no reads in 1,1)

# Pick abundance threshold (e.g. 0.1% relative abundance)
threshold <- 0.001  

# Convert to relative abundances
rel <- sweep(seqtab.filt3, 1, rowSums(seqtab.filt3), "/")

# Count number of samples where ASV is above threshold
present_counts <- colSums(rel >= threshold)

# Keep ASVs present in at least 10 samples
keep_asvs <- names(present_counts[present_counts >= 10])

# Filter table
seqtab.filt4 <- seqtab.filt3[, keep_asvs]











# In each time point, run deseq2 vs treatment
library(DESeq2)
library(apeglm)

# NOTE: for ASV mapping to work, make sure phyloseq object
# is generated from the same as asv_mapping (seqtab.filt3)

timepoints <- c("Wk1", "Wk4", "Wk8", "Wk12", "Wk16")
dds_list <- list()   # store full DESeqDataSet for each time point

for(tp in timepoints) {
  # Subset samples for current time point
  # samples_to_keep <- sample_data(ps)$timepoint %in% c("Day1", tp)
  samples_to_keep <- sample_data(ps)$timepoint %in% c(tp)
  ps_sub <- prune_samples(samples_to_keep, ps)

  # Convert variables to factors
  sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, levels = c("placebo", "active"))
  # sample_data(ps_sub)$time <- factor(sample_data(ps_sub)$timepoint, levels = c("Day1", tp))

  # Convert to DESeq2 object
  dds_sub <- phyloseq_to_deseq2(ps_sub, ~ treatment)

  # Run DESeq
  dds_sub <- DESeq(dds_sub)

  # Store full results object
  dds_list[[tp]] <- dds_sub
}

for(tp in names(dds_list)) {
  dds_sub <- dds_list[[tp]]

  # Extract results
  res_df <- as.data.frame(results(dds_sub))
  res_df$ASV <- rownames(res_df)
}


