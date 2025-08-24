library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)

working_dir <- "/data/local/jy1008/Allegretti"
meta_data <- read.table(file.path(working_dir, "metadata", 
                "meta_data_subject.csv"), header = TRUE,
                sep = ",", stringsAsFactors = FALSE)

time_data <- read.table(file.path(working_dir, "metadata", 
                "meta_time_series_mayo_scores.csv"), header = TRUE,
                sep = ",", stringsAsFactors = FALSE)
time_data <- time_data %>%
    rename(sample_id = subject_id)

meta_sub <- time_data[, c("sample_id", "total_partial_mayo")]
meta_sub <- meta_sub %>%
    extract(sample_id, into = c("subject_id", "timepoint"),
            regex = "(CP\\d+)([A-Za-z]+\\d+)", remove = FALSE)
meta_sub$subject_id <- str_replace(meta_sub$subject_id,
                                "^(CP\\d{3})(\\d{2})$", "\\1-\\2")
meta_sub <- merge(meta_sub, 
                  meta_data[, c("subject_id", "treatment", "age", "sex")], 
                  by = "subject_id", all.x = TRUE)

# convert all time points for display
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

# mean and sd per treatment and time
summary_meta <- meta_sub %>%
    group_by(treatment, time_numeric) %>%
    summarise(
        mean_score = mean(total_partial_mayo, na.rm = TRUE),
        sd_score = sd(total_partial_mayo, na.rm = TRUE),
        .groups = "drop"
    )


# Plot individual subjects + mean ± sd ribbon per treatment
p <- ggplot() +
  # Individual subject lines (optional)
  geom_line(data = meta_sub, aes(x = time_numeric, y = total_partial_mayo, group = subject_id, color = treatment), alpha = 0.3) +
  geom_point(data = meta_sub, aes(x = time_numeric, y = total_partial_mayo, group = subject_id, color = treatment), alpha = 0.3) +

  # Mean line per treatment
  geom_line(data = summary_meta, aes(x = time_numeric, y = mean_score, color = treatment), size = 1.2) +
  
  # Ribbon for mean ± sd
  geom_ribbon(data = summary_meta, aes(x = time_numeric,
                                     ymin = mean_score - sd_score,
                                     ymax = mean_score + sd_score,
                                     fill = treatment),
              alpha = 0.2, color = NA) +
  
  theme_minimal() +
  labs(
    title = "Total partial Mayo score over time by treatment",
    x = "Week",
    y = "Total partial Mayo score",
    color = "Treatment",
    fill = "Treatment"
  )
ggsave("fig1_time_series_overview.pdf",
        plot = p, width = 8, height = 6)

# -----
# Each subject's time series
meta_sub$subject_id_disp <- gsub("^CP101-", "", meta_sub$subject_id)

# Label positions
label_data <- meta_sub %>%
  group_by(subject_id_disp, treatment) %>%
  summarise(
    x = min(time_numeric, na.rm = TRUE),
    y = max(total_partial_mayo, na.rm = TRUE),
    .groups = "drop"
  )

# Create an ordering factor so placebo/active alternate columns
meta_sub <- meta_sub %>%
  arrange(treatment, subject_id_disp) %>%
  mutate(panel_id = forcats::fct_inorder(paste(treatment, subject_id_disp, sep = "_")))

label_data <- label_data %>%
  arrange(treatment, subject_id_disp) %>%
  mutate(panel_id = forcats::fct_inorder(paste(treatment, subject_id_disp, sep = "_")))

# Plot
p <- ggplot(meta_sub, aes(x = time_numeric, y = total_partial_mayo, color = treatment)) +
  geom_line(alpha = 1.0) +
  geom_point(alpha = 1.0) +
  geom_text(
    data = label_data,
    aes(x = Inf, y = Inf, label = subject_id_disp),
    inherit.aes = FALSE,
    hjust = 1.1,  # push slightly inside from right edge
    vjust = 1.1,  # push slightly inside from top edge
    size = 3,
    color = "black"
  ) +
  theme_bw() +
  facet_wrap(~ panel_id, ncol = 2) +
  labs(
    title = "Total partial Mayo score over time by treatment",
    x = "Week",
    y = "Total partial Mayo score",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text = element_blank()
  )
ggsave("fig1_time_series_individual_subjects.pdf",
        plot = p, width = 4, height = 4)
# -----


# subject-level regression
library(lme4)
library(lmerTest)  # for p-values

meta_sub$subject_id <- as.factor(meta_sub$subject_id)
meta_sub$treatment <- factor(meta_sub$treatment, 
                              levels = c("placebo", "active"))
# with random slope
# model <- lmer(total_partial_mayo ~ 
#                 time_numeric * treatment + 
#                 (1 + time_numeric | subject_id),
#               data = meta_sub)
# with random intercept only
model <- lmer(total_partial_mayo ~ 
                time_numeric * treatment + 
                (1 | subject_id),
              data = meta_sub)
summary(model)

# No significant interaction between time and treatment

# Linear mixed model fit by REML. t-tests use Satterthwaite's method [
# lmerModLmerTest]
# Formula: total_partial_mayo ~ time_numeric * treatment + (1 | subject_id)
#    Data: meta_sub
# 
# REML criterion at convergence: 290.1
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.99784 -0.69329  0.04912  0.65356  2.76381 
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  subject_id (Intercept) 1.561    1.249   
#  Residual               1.481    1.217   
# Number of obs: 77, groups:  subject_id, 16
# 
# Fixed effects:
#                               Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)                   4.868243   0.521614 19.124355   9.333 1.49e-08
# time_numeric                 -0.016757   0.004605 65.758917  -3.639  0.00054
# treatmentactive              -0.482508   0.728749 18.305022  -0.662  0.51615
# time_numeric:treatmentactive  0.005256   0.005901 64.770902   0.891  0.37637
#                                 
# (Intercept)                  ***
# time_numeric                 ***
# treatmentactive                 
# time_numeric:treatmentactive    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) tm_nmr trtmnt
# time_numerc -0.290              
# treatmntctv -0.716  0.208       
# tm_nmrc:trt  0.227 -0.780 -0.303
