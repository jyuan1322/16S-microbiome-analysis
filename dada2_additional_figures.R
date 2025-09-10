library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

df <- read.table("/data/local/jy1008/Allegretti/test_data_output_v2/dada2_stats.tsv", header = TRUE, sep = "\t")
df$sample_name = rownames(df)
df$subject_id <- gsub("^CP101(\\d+)(Day|Wk).*", "\\1", df$sample_name)
df$time_point <- gsub("^CP101\\d+", "", df$sample_name)

p <- ggplot(df, aes(x = input, y = pct_retained, color = subject_id)) +
  geom_point(size = 3) +
  labs(x = "Input reads", y = "Percent retained", title = "Input reads vs Percent Retained") +
  theme_bw() +
  # Ensure the lowest points have ticks
  scale_x_continuous(expand = expansion(add = 10000)) +  # add tiny buffer
  scale_y_continuous(expand = expansion(add = 0.3))

ggsave("dada2_input_vs_pct_retained.pdf", plot = p, width = 6, height = 4)


# Visualize the read counts from dada2_stats.tsv in a bar plot. This is for checking the read quality.

# Order time_point
df$time_point <- factor(df$time_point, levels = c("Day1", "Wk1", "Wk4", "Wk8", "Wk12", "Wk16"))

df <- df[order(df$subject_id, df$time_point), ]
# df$subject_time_label <- paste(df$subject_id, df$time_point, sep = " ")

# Create subject_time_label as an ordered factor
df$subject_time_label <- factor(
  paste(df$subject_id, df$time_point, sep = " "),
  levels = unique(paste(df$subject_id, df$time_point, sep = " "))
)

# Grouped bar plot: x = time_point, color/fill = subject
p <- ggplot(df, aes(x = subject_time_label, y = input, fill = subject_id)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.9) +
  geom_hline(yintercept = 10000, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  labs(x = "Time Point", y = "Input Reads", fill = "Subject") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(option = "turbo")
ggsave("dada2_input_read_count.pdf", plot = p, width = 12, height = 4)


# bar plot with more metrics
# Pivot longer: convert counts into "metric"/"value"
df_long <- df %>%
  pivot_longer(
    cols = c(input, filtered, merged, nonchim), 
    names_to = "metric", 
    values_to = "count"
  )

# Make sure df_long$metric is ordered
df_long$metric <- factor(
  df_long$metric,
  levels = c("input", "filtered", "merged", "nonchim")
)

# Stacked bar plot
p <- ggplot(df_long, aes(x = subject_time_label, y = count, fill = metric)) +
  geom_col(position = "identity", alpha = 1.0) +
  theme_bw() +
  labs(x = "Sample", y = "Read Count", fill = "Processing Step") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("dada2_filtering_read_count_bars.pdf", plot = p, width = 12, height = 4)



# overall read quality per cycle
qc_df <- read.csv("/data/local/jy1008/Allegretti/test_data_output_v2/read_quality_summary.csv",
                  header = TRUE, stringsAsFactors = FALSE)
p <- ggplot(qc_df, aes(x = Cycle)) +
  geom_line(aes(y = Mean_All), color = "blue", size = 1) +                # mean line
  geom_ribbon(aes(ymin = Mean_All - SD_All, ymax = Mean_All + SD_All),   # shaded area
              fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  expand_limits(y = 20) +  # force y-axis to start at 0
  labs(x = "Cycle", y = "Quality Score", title = "Mean Quality ± SD Across Cycles")
ggsave("dada2_read_quality_across_cycles.pdf", plot = p, width = 6, height = 4)