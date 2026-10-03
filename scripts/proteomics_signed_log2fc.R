# ============================================================
# Signed log2 fold change — CCR7+ DC proteomics
# Same load, filters, and sample-wise median centering as
# proteomics_divergence_analysis.R
# ============================================================
# Sign convention (Group1 minus Group2):
#   Temporal: later day minus earlier day
#     positive = higher at the later day
#   Within day: GFP+ minus GFP-
#     positive = higher in GFP+
# Place next to: Proteomic data Luis.xlsx
# Output: Supplementary_Table_Signed_Log2FC.xlsx
# ============================================================

library(dplyr)
library(tidyr)
library(readxl)
library(tibble)
library(openxlsx)

data_path <- "Proteomic data Luis.xlsx"

raw <- read_excel(data_path) %>%
  filter(is.na(`Potential contaminant`) | `Potential contaminant` != "YES") %>%
  select(
    `Gene names`, id,
    starts_with("10_"), starts_with("21_"), starts_with("28_"),
    `ANOVA q-value`
  )

intensity_cols <- names(raw)[grepl("^(10_|21_|28_)", names(raw))]

proteomics_data <- raw %>%
  filter(if_all(all_of(intensity_cols), ~ !is.na(.))) %>%
  mutate(across(all_of(intensity_cols), ~ . - median(., na.rm = TRUE)))

protein_info <- proteomics_data %>% select(`Gene names`, id)

full_data <- proteomics_data %>%
  select(-`ANOVA q-value`) %>%
  unite("row_name", `Gene names`, id, sep = "_") %>%
  column_to_rownames("row_name") %>%
  as.matrix()

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

get_mean_profile <- function(sample_names, data_mat) {
  if (length(sample_names) > 1) {
    rowMeans(data_mat[, sample_names, drop = FALSE], na.rm = TRUE)
  } else {
    data_mat[, sample_names]
  }
}

signed <- function(later, earlier) later - earlier

GFP_pos_28vs10 <- signed(get_mean_profile(pos_samples[["28"]], full_data), get_mean_profile(pos_samples[["10"]], full_data))
GFP_neg_28vs10 <- signed(get_mean_profile(neg_groups[["28"]], full_data), get_mean_profile(neg_groups[["10"]], full_data))
GFP_pos_21vs10 <- signed(get_mean_profile(pos_samples[["21"]], full_data), get_mean_profile(pos_samples[["10"]], full_data))
GFP_neg_21vs10 <- signed(get_mean_profile(neg_groups[["21"]], full_data), get_mean_profile(neg_groups[["10"]], full_data))
GFP_pos_28vs21 <- signed(get_mean_profile(pos_samples[["28"]], full_data), get_mean_profile(pos_samples[["21"]], full_data))
GFP_neg_28vs21 <- signed(get_mean_profile(neg_groups[["28"]], full_data), get_mean_profile(neg_groups[["21"]], full_data))

Day10_pos_vs_neg <- signed(get_mean_profile(pos_samples[["10"]], full_data), get_mean_profile(neg_groups[["10"]], full_data))
Day21_pos_vs_neg <- signed(get_mean_profile(pos_samples[["21"]], full_data), get_mean_profile(neg_groups[["21"]], full_data))
Day28_pos_vs_neg <- signed(get_mean_profile(pos_samples[["28"]], full_data), get_mean_profile(neg_groups[["28"]], full_data))

signed_table <- data.frame(
  Gene_names = protein_info$`Gene names`,
  Protein_ID = protein_info$id,
  GFP_pos_Day28_minus_Day10 = GFP_pos_28vs10,
  GFP_neg_Day28_minus_Day10 = GFP_neg_28vs10,
  GFP_pos_Day21_minus_Day10 = GFP_pos_21vs10,
  GFP_neg_Day21_minus_Day10 = GFP_neg_21vs10,
  GFP_pos_Day28_minus_Day21 = GFP_pos_28vs21,
  GFP_neg_Day28_minus_Day21 = GFP_neg_28vs21,
  GFP_pos_minus_GFP_neg_Day10 = Day10_pos_vs_neg,
  GFP_pos_minus_GFP_neg_Day21 = Day21_pos_vs_neg,
  GFP_pos_minus_GFP_neg_Day28 = Day28_pos_vs_neg
)

wb <- createWorkbook()
addWorksheet(wb, "Signed_Log2FC")
writeData(wb, "Signed_Log2FC", signed_table)
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
addStyle(wb, "Signed_Log2FC", headerStyle, rows = 1, cols = 1:ncol(signed_table), gridExpand = TRUE)
setColWidths(wb, "Signed_Log2FC", cols = 1:ncol(signed_table), widths = "auto")
saveWorkbook(wb, "Supplementary_Table_Signed_Log2FC.xlsx", overwrite = TRUE)

message("Wrote Supplementary_Table_Signed_Log2FC.xlsx (", nrow(signed_table), " proteins)")
message("Sign: temporal = later minus earlier; within day = GFP+ minus GFP-")
