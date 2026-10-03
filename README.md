# Proteomic Divergence Analysis – CCR7⁺ Dendritic Cells

Revised proteomic divergence analysis for the Nature Communications manuscript (CCR7⁺ dendritic cells). Full QC proteome, no significance pre-filter.

## What this analysis does

- Sample-wise median centering of LFQ intensities (ANOVA q-value excluded).
- Per-protein `|Δlog2|` from day-mean profiles, plotted as rainclouds. Points are proteins.
- Signed log2FC table (later day minus earlier day; within day, GFP+ minus GFP−).
- Primary statistic: difference of median `|Δlog2|` with a protein bootstrap 95% CI (B = 2000).
- Secondary: two-sided KS tests, BH FDR across the 12 planned contrasts.

Each LFQ column is one donor pool. GFP+ and GFP− are sorts from that pool.

- Male: day 10 n=19; day 21-1 n=7; day 21-2 n=6; day 28 n=8
- Female: day 10 n=18; day 21-1 n=6; day 21-2 n=3; day 28-1 n=6; day 28-2 n=6

## Repository structure

```
scripts/proteomics_divergence_analysis.R      # centering, |Δlog2| table, four figures
scripts/proteomics_divergence_effect_sizes.R  # median differences, bootstrap CIs, BH
scripts/proteomics_signed_log2fc.R            # signed log2FC table
docs/Methods.md
```

## How to run

Place `Proteomic data Luis.xlsx` in the working directory.

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "readxl", "tibble", "ggsignif", "openxlsx"))
source("scripts/proteomics_divergence_analysis.R")
source("scripts/proteomics_signed_log2fc.R")
# effect sizes read the |Δlog2| workbook just written:
source("scripts/proteomics_divergence_effect_sizes.R")
```

## Outputs

- `Supplementary_Table_Individual_Protein_Values.xlsx` — per-protein `|Δlog2|`
- `Supplementary_Table_Signed_Log2FC.xlsx` — per-protein signed log2FC
- `Proteomic_Divergence_EffectSize_BH.csv` — 12 contrasts
- `Figure_3E_GFP+_vs_GFP-_Temporal.png`
- `Figure_GFP+_Temporal_Changes.png`
- `Figure_GFP-_Temporal_Changes.png`
- `Figure_GFP+_vs_GFP-_per_Day.png`

Figure brackets are still raw KS stars. Report the median differences and CIs from the CSV as the result.

## Author

Luis Castillo Montanez  
GitHub [@luiscastillo718-cell](https://github.com/luiscastillo718-cell)
