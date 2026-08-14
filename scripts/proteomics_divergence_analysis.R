# ============================================================
# COMPLETE STANDALONE SCRIPT
# Proteomic Divergence Analysis – CCR7⁺ Dendritic Cells
# Nature Communications revision (Reviewer #2 response)
# ============================================================
# Generates:
#   1. All 4 main graphs (Y-axis fixed 0–7)
#   2. Supplementary Table with individual |\u0394log₂| values per protein
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(tibble)
library(ggsignif)
library(openxlsx)

# ====================== 1. LOAD DATA ======================
data_path <- "Proteomic data Luis.xlsx"

proteomics_data <- read_excel(data_path) %>%
  filter(is.na(`Potential contaminant`) | `Potential contaminant` != "YES") %>%
  select(`Gene names`, id, 
         starts_with("10_"), starts_with("21_"), starts_with("28_"),
         `ANOVA q-value`) %>%
  filter(complete.cases(across(where(is.numeric)))) %>%
  mutate(across(where(is.numeric), 
                ~ . - median(c_across(where(is.numeric)), na.rm = TRUE)))

protein_info <- proteomics_data %>%
  select(`Gene names`, id)

full_data <- proteomics_data %>%
  select(-`ANOVA q-value`) %>%
  unite("row_name", `Gene names`, id, sep = "_") %>%
  column_to_rownames("row_name") %>%
  as.matrix()

# ====================== 2. GROUPS ======================
pos_samples <- list(
  "10" = c("10_M_pos_001", "10_F_pos_003"),
  "21" = c("21_M_pos_005", "21_M_pos_009", "21_F_pos_007", "21_F_pos2_011"),
  "28" = c("28_M_pos_013", "28_F_pos_015", "28_F_pos2_017")
)

neg_groups <- list(
  "10" = c("10_M_neg_002", "10_F_neg_004"),
  "21" = c("21_M_neg_006", "21_M_neg2_010", "21_F_neg2_012"),
  "28" = c("28_M_neg_014", "28_F_neg_016", "28_F_neg2_018")
)

# ====================== 3. HELPER FUNCTIONS ======================
get_mean_profile <- function(sample_names, data_mat) {
  if (length(sample_names) > 1) {
    rowMeans(data_mat[, sample_names, drop = FALSE], na.rm = TRUE)
  } else {
    data_mat[, sample_names]
  }
}

get_star <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ====================== 4. CALCULATE ALL |\u0394log₂| VALUES ======================
GFP_pos_28vs10 <- abs(get_mean_profile(pos_samples[["28"]], full_data) - get_mean_profile(pos_samples[["10"]], full_data))
GFP_neg_28vs10 <- abs(get_mean_profile(neg_groups[["28"]], full_data) - get_mean_profile(neg_groups[["10"]], full_data))
GFP_pos_21vs10 <- abs(get_mean_profile(pos_samples[["21"]], full_data) - get_mean_profile(pos_samples[["10"]], full_data))
GFP_neg_21vs10 <- abs(get_mean_profile(neg_groups[["21"]], full_data) - get_mean_profile(neg_groups[["10"]], full_data))
GFP_pos_28vs21 <- abs(get_mean_profile(pos_samples[["28"]], full_data) - get_mean_profile(pos_samples[["21"]], full_data))
GFP_neg_28vs21 <- abs(get_mean_profile(neg_groups[["28"]], full_data) - get_mean_profile(neg_groups[["21"]], full_data))

Day10_pos_vs_neg <- abs(get_mean_profile(pos_samples[["10"]], full_data) - get_mean_profile(neg_groups[["10"]], full_data))
Day21_pos_vs_neg <- abs(get_mean_profile(pos_samples[["21"]], full_data) - get_mean_profile(neg_groups[["21"]], full_data))
Day28_pos_vs_neg <- abs(get_mean_profile(pos_samples[["28"]], full_data) - get_mean_profile(neg_groups[["28"]], full_data))

# ====================== 5. SUPPLEMENTARY TABLE ======================
cat("Creating Supplementary Table...\n")

suppl_table <- data.frame(
  Gene_names               = protein_info$`Gene names`,
  Protein_ID               = protein_info$id,
  GFP_pos_Day28_vs_10      = GFP_pos_28vs10,
  GFP_neg_Day28_vs_10      = GFP_neg_28vs10,
  GFP_pos_Day21_vs_10      = GFP_pos_21vs10,
  GFP_neg_Day21_vs_10      = GFP_neg_21vs10,
  GFP_pos_Day28_vs_21      = GFP_pos_28vs21,
  GFP_neg_Day28_vs_21      = GFP_neg_28vs21,
  GFP_pos_vs_GFP_neg_Day10 = Day10_pos_vs_neg,
  GFP_pos_vs_GFP_neg_Day21 = Day21_pos_vs_neg,
  GFP_pos_vs_GFP_neg_Day28 = Day28_pos_vs_neg
)

wb <- createWorkbook()
addWorksheet(wb, "Individual_Protein_Values")
writeData(wb, "Individual_Protein_Values", suppl_table)

headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
addStyle(wb, "Individual_Protein_Values", headerStyle, rows = 1, cols = 1:ncol(suppl_table), gridExpand = TRUE)
setColWidths(wb, "Individual_Protein_Values", cols = 1:ncol(suppl_table), widths = "auto")

saveWorkbook(wb, "Supplementary_Table_Individual_Protein_Values.xlsx", overwrite = TRUE)
cat("\u2705 Supplementary table saved: Supplementary_Table_Individual_Protein_Values.xlsx\n\n")

# ====================== 6. GRAPH 1: Combined ======================
cat("Generating Graph 1: Combined GFP+ vs GFP- Temporal...\n")

temporal_comparisons <- list(
  "GFP+ Day 28 vs 10" = GFP_pos_28vs10,
  "GFP- Day 28 vs 10" = GFP_neg_28vs10,
  "GFP+ Day 21 vs 10" = GFP_pos_21vs10,
  "GFP- Day 21 vs 10" = GFP_neg_21vs10,
  "GFP+ Day 28 vs 21" = GFP_pos_28vs21,
  "GFP- Day 28 vs 21" = GFP_neg_28vs21
)

plot_data1 <- bind_rows(lapply(names(temporal_comparisons), function(name) {
  tibble(Comparison = name, Value = temporal_comparisons[[name]])
}))

plot_data1$Comparison <- factor(plot_data1$Comparison, 
                                levels = c("GFP+ Day 28 vs 10", "GFP- Day 28 vs 10",
                                           "GFP+ Day 21 vs 10", "GFP- Day 21 vs 10",
                                           "GFP+ Day 28 vs 21", "GFP- Day 28 vs 21"))
plot_data1$Fill <- ifelse(grepl("^GFP\\+", plot_data1$Comparison), "#2E8B57", "white")

ks_gfp_vs_gfpneg <- c(
  ks.test(temporal_comparisons[["GFP+ Day 28 vs 10"]], temporal_comparisons[["GFP- Day 28 vs 10"]])$p.value,
  ks.test(temporal_comparisons[["GFP+ Day 21 vs 10"]], temporal_comparisons[["GFP- Day 21 vs 10"]])$p.value,
  ks.test(temporal_comparisons[["GFP+ Day 28 vs 21"]], temporal_comparisons[["GFP- Day 28 vs 21"]])$p.value
)

p1 <- ggplot(plot_data1, aes(x = Comparison, y = Value, fill = Fill)) +
  geom_jitter(width = 0.12, size = 0.22, alpha = 0.65, color = "#333333", stroke = 0.3) +
  geom_violin(trim = FALSE, alpha = 0.55, color = "black", linewidth = 0.4, width = 0.7) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, color = "black", linewidth = 0.5) +
  geom_signif(
    comparisons = list(
      c("GFP+ Day 28 vs 10", "GFP- Day 28 vs 10"),
      c("GFP+ Day 21 vs 10", "GFP- Day 21 vs 10"),
      c("GFP+ Day 28 vs 21", "GFP- Day 28 vs 21")
    ),
    annotations = sapply(ks_gfp_vs_gfpneg, get_star),
    y_position = c(5.5, 5.9, 6.3),
    tip_length = 0.015,
    textsize = 4.0
  ) +
  scale_fill_identity() +
  scale_y_continuous(breaks = 0:7, limits = c(0, 7), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_minimal(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(angle = 25, hjust = 1), 
        legend.position = "none") +
  labs(title = "GFP+ vs GFP\u2013 Temporal Comparison", 
       subtitle = "Full QC proteome (n = 3404 proteins)",
       y = "|\u0394log\u2082 fold change|")

ggsave("Figure_3E_GFP+_vs_GFP-_Temporal.png", p1, width = 13, height = 6.5, dpi = 300, bg = "white")
cat("Saved: Figure_3E_GFP+_vs_GFP-_Temporal.png\n")

# ====================== 7. GRAPH 2: GFP+ Temporal ======================
cat("Generating Graph 2: GFP+ Temporal Changes...\n")

gfp_temp <- list(
  "Day 28 vs 10" = GFP_pos_28vs10,
  "Day 21 vs 10" = GFP_pos_21vs10,
  "Day 28 vs 21" = GFP_pos_28vs21
)

plot_data2 <- bind_rows(lapply(names(gfp_temp), function(name) {
  tibble(Comparison = name, Value = gfp_temp[[name]])
}))
plot_data2$Comparison <- factor(plot_data2$Comparison, levels = c("Day 28 vs 10", "Day 21 vs 10", "Day 28 vs 21"))

ks_gfp <- c(
  ks.test(gfp_temp[["Day 28 vs 10"]], gfp_temp[["Day 21 vs 10"]])$p.value,
  ks.test(gfp_temp[["Day 28 vs 10"]], gfp_temp[["Day 28 vs 21"]])$p.value,
  ks.test(gfp_temp[["Day 21 vs 10"]], gfp_temp[["Day 28 vs 21"]])$p.value
)

p2 <- ggplot(plot_data2, aes(x = Comparison, y = Value, fill = "#2E8B57")) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.7, color = "#333333") +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.4, width = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.95, color = "black", linewidth = 0.5) +
  geom_signif(comparisons = list(c("Day 28 vs 10", "Day 21 vs 10"), 
                                 c("Day 28 vs 10", "Day 28 vs 21"),
                                 c("Day 21 vs 10", "Day 28 vs 21")),
              annotations = sapply(ks_gfp, get_star),
              y_position = c(5.5, 5.9, 6.3),
              tip_length = 0.015,
              textsize = 4.0) +
  scale_fill_identity() +
  scale_y_continuous(breaks = 0:7, limits = c(0, 7), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_minimal(base_size = 12) +
  labs(title = "GFP+ Temporal Changes", y = "|\u0394log\u2082 fold change|")

ggsave("Figure_GFP+_Temporal_Changes.png", p2, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: Figure_GFP+_Temporal_Changes.png\n")

# ====================== 8. GRAPH 3: GFP- Temporal (WHITE) ======================
cat("Generating Graph 3: GFP- Temporal Changes (White)...\n")

gfpneg_temp <- list(
  "Day 28 vs 10" = GFP_neg_28vs10,
  "Day 21 vs 10" = GFP_neg_21vs10,
  "Day 28 vs 21" = GFP_neg_28vs21
)

plot_data3 <- bind_rows(lapply(names(gfpneg_temp), function(name) {
  tibble(Comparison = name, Value = gfpneg_temp[[name]])
}))
plot_data3$Comparison <- factor(plot_data3$Comparison, levels = c("Day 28 vs 10", "Day 21 vs 10", "Day 28 vs 21"))

ks_gfpneg <- c(
  ks.test(gfpneg_temp[["Day 28 vs 10"]], gfpneg_temp[["Day 21 vs 10"]])$p.value,
  ks.test(gfpneg_temp[["Day 28 vs 10"]], gfpneg_temp[["Day 28 vs 21"]])$p.value,
  ks.test(gfpneg_temp[["Day 21 vs 10"]], gfpneg_temp[["Day 28 vs 21"]])$p.value
)

p3 <- ggplot(plot_data3, aes(x = Comparison, y = Value)) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.7, color = "#333333") +
  geom_violin(trim = FALSE, fill = "white", color = "black", alpha = 0.7, linewidth = 0.4, width = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", color = "black", alpha = 0.95, linewidth = 0.5) +
  geom_signif(comparisons = list(c("Day 28 vs 10", "Day 21 vs 10"), 
                                 c("Day 28 vs 10", "Day 28 vs 21"),
                                 c("Day 21 vs 10", "Day 28 vs 21")),
              annotations = sapply(ks_gfpneg, get_star),
              y_position = c(5.5, 5.9, 6.3),
              tip_length = 0.015,
              textsize = 4.0) +
  scale_y_continuous(breaks = 0:7, limits = c(0, 7), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_minimal(base_size = 12) +
  labs(title = "GFP- Temporal Changes", y = "|\u0394log\u2082 fold change|")

ggsave("Figure_GFP-_Temporal_Changes.png", p3, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: Figure_GFP-_Temporal_Changes.png\n")

# ====================== 9. GRAPH 4: GFP+ vs GFP- per Day ======================
cat("Generating Graph 4: GFP+ vs GFP- per Day (Lighter Green)...\n")

day_diffs <- list(
  "Day 10" = Day10_pos_vs_neg,
  "Day 21" = Day21_pos_vs_neg,
  "Day 28" = Day28_pos_vs_neg
)

plot_data4 <- bind_rows(lapply(names(day_diffs), function(name) {
  tibble(Day = name, Value = day_diffs[[name]])
}))
plot_data4$Day <- factor(plot_data4$Day, levels = c("Day 10", "Day 21", "Day 28"))

ks_day <- c(
  ks.test(day_diffs[["Day 10"]], day_diffs[["Day 21"]])$p.value,
  ks.test(day_diffs[["Day 21"]], day_diffs[["Day 28"]])$p.value,
  ks.test(day_diffs[["Day 10"]], day_diffs[["Day 28"]])$p.value
)

p4 <- ggplot(plot_data4, aes(x = Day, y = Value, fill = "#81C784")) +
  geom_jitter(width = 0.18, size = 0.25, alpha = 0.7, color = "#333333") +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", linewidth = 0.4, width = 0.75) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.95, color = "black", linewidth = 0.5) +
  geom_signif(comparisons = list(c("Day 10", "Day 21"), c("Day 21", "Day 28"), c("Day 10", "Day 28")),
              annotations = sapply(ks_day, get_star),
              y_position = c(5.5, 5.9, 6.3),
              tip_length = 0.015,
              textsize = 4.0) +
  scale_fill_identity() +
  scale_y_continuous(breaks = 0:7, limits = c(0, 7), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_minimal(base_size = 12) +
  labs(title = "GFP+ vs GFP\u2013 at Each Timepoint", y = "|\u0394log\u2082 fold change| (GFP+ vs GFP\u2013)")

ggsave("Figure_GFP+_vs_GFP-_per_Day.png", p4, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: Figure_GFP+_vs_GFP-_per_Day.png\n")

cat("\n\u2705 All outputs generated successfully!\n")
