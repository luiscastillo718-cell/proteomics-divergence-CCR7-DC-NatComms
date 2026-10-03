# ============================================================
# Phase B — effect sizes + BH FDR for planned contrasts
# Proteomic divergence, CCR7+ DCs (Nat Comms R2)
# ============================================================
# Does NOT re-read the raw LFQ workbook.
# Input (same folder): Supplementary_Table_Individual_Protein_Values.xlsx
#   produced by the sample-wise median-centered script.
# Output: Proteomic_Divergence_EffectSize_BH.csv
#
# Effect size = difference of medians of |Δlog2| distributions
#   (first named group minus second).
# CI = percentile bootstrap over proteins (B = 2000).
# Proteins are correlated, so the CI is descriptive, not a
# design-based confidence interval. Primary evidence is the
# median difference; KS p and BH q are secondary.
#
# Biological pools (mice per pool; GFP+ and GFP- are sorts
# from these pools, not extra biological n):
#   Male   D10 n=19; D21-1 n=7; D21-2 n=6; D28 n=8
#   Female D10 n=18; D21-1 n=6; D21-2 n=3; D28-1 n=6; D28-2 n=6
# ============================================================

library(readxl)
library(openxlsx)

set.seed(20261003)

in_path  <- "Supplementary_Table_Individual_Protein_Values.xlsx"
out_path <- "Proteomic_Divergence_EffectSize_BH.csv"
B <- 2000L

tab <- as.data.frame(read_excel(in_path, sheet = 1))
stopifnot(nrow(tab) > 100)

median_diff <- function(a, b) median(a, na.rm = TRUE) - median(b, na.rm = TRUE)

bootstrap_ci <- function(a, b, B) {
  n <- length(a)
  stopifnot(length(b) == n)
  diffs <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    diffs[i] <- median_diff(a[idx], b[idx])
  }
  qs <- quantile(diffs, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(lo = qs[1], hi = qs[2])
}

# Planned contrasts, same pairs as the raincloud brackets.
# contrast_set groups tests that share a figure.
pairs <- list(
  list(set = "GFP+ vs GFP- temporal",
       label = "Day 28 vs 10: GFP+ minus GFP-",
       a = "GFP_pos_Day28_vs_10", b = "GFP_neg_Day28_vs_10"),
  list(set = "GFP+ vs GFP- temporal",
       label = "Day 21 vs 10: GFP+ minus GFP-",
       a = "GFP_pos_Day21_vs_10", b = "GFP_neg_Day21_vs_10"),
  list(set = "GFP+ vs GFP- temporal",
       label = "Day 28 vs 21: GFP+ minus GFP-",
       a = "GFP_pos_Day28_vs_21", b = "GFP_neg_Day28_vs_21"),
  list(set = "GFP+ temporal",
       label = "GFP+: (Day 28 vs 10) minus (Day 21 vs 10)",
       a = "GFP_pos_Day28_vs_10", b = "GFP_pos_Day21_vs_10"),
  list(set = "GFP+ temporal",
       label = "GFP+: (Day 28 vs 10) minus (Day 28 vs 21)",
       a = "GFP_pos_Day28_vs_10", b = "GFP_pos_Day28_vs_21"),
  list(set = "GFP+ temporal",
       label = "GFP+: (Day 21 vs 10) minus (Day 28 vs 21)",
       a = "GFP_pos_Day21_vs_10", b = "GFP_pos_Day28_vs_21"),
  list(set = "GFP- temporal",
       label = "GFP-: (Day 28 vs 10) minus (Day 21 vs 10)",
       a = "GFP_neg_Day28_vs_10", b = "GFP_neg_Day21_vs_10"),
  list(set = "GFP- temporal",
       label = "GFP-: (Day 28 vs 10) minus (Day 28 vs 21)",
       a = "GFP_neg_Day28_vs_10", b = "GFP_neg_Day28_vs_21"),
  list(set = "GFP- temporal",
       label = "GFP-: (Day 21 vs 10) minus (Day 28 vs 21)",
       a = "GFP_neg_Day21_vs_10", b = "GFP_neg_Day28_vs_21"),
  list(set = "GFP+ vs GFP- by day",
       label = "Within-day |GFP+ vs GFP-|: Day 10 minus Day 21",
       a = "GFP_pos_vs_GFP_neg_Day10", b = "GFP_pos_vs_GFP_neg_Day21"),
  list(set = "GFP+ vs GFP- by day",
       label = "Within-day |GFP+ vs GFP-|: Day 21 minus Day 28",
       a = "GFP_pos_vs_GFP_neg_Day21", b = "GFP_pos_vs_GFP_neg_Day28"),
  list(set = "GFP+ vs GFP- by day",
       label = "Within-day |GFP+ vs GFP-|: Day 10 minus Day 28",
       a = "GFP_pos_vs_GFP_neg_Day10", b = "GFP_pos_vs_GFP_neg_Day28")
)

need <- unique(unlist(lapply(pairs, function(p) c(p$a, p$b))))
missing <- setdiff(need, names(tab))
if (length(missing)) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

rows <- lapply(pairs, function(p) {
  a <- tab[[p$a]]
  b <- tab[[p$b]]
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]
  b <- b[ok]
  est <- median_diff(a, b)
  ci <- bootstrap_ci(a, b, B)
  ks <- ks.test(a, b)
  data.frame(
    contrast_set = p$set,
    contrast = p$label,
    n_proteins = length(a),
    median_A = median(a),
    median_B = median(b),
    delta_median = est,
    boot_ci_low = ci["lo"],
    boot_ci_high = ci["hi"],
    ks_D = unname(ks$statistic),
    ks_p_raw = ks$p.value,
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
out$ks_q_BH <- p.adjust(out$ks_p_raw, method = "BH")
out$ci_excludes_0 <- out$boot_ci_low > 0 | out$boot_ci_high < 0

# Stable column order
out <- out[, c(
  "contrast_set", "contrast", "n_proteins",
  "median_A", "median_B", "delta_median",
  "boot_ci_low", "boot_ci_high", "ci_excludes_0",
  "ks_D", "ks_p_raw", "ks_q_BH"
)]

write.csv(out, out_path, row.names = FALSE)
message("Wrote ", out_path, " (", nrow(out), " contrasts, B = ", B, ")")
print(out[, c("contrast", "delta_median", "boot_ci_low", "boot_ci_high", "ks_p_raw", "ks_q_BH")],
      row.names = FALSE)
