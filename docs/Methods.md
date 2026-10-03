# Methods – Proteomic Divergence Analysis

## Proteomic data processing
Label-free quantification intensities were filtered to remove potential contaminants. Proteins with a missing value in any sample were excluded (n = 3404). Log2 intensities were median-centered within each sample (column median of the LFQ intensity columns only). The ANOVA q-value was not used for filtering or for centering.

Each LFQ sample is one donor pool. GFP+ and GFP− are sorts from that pool, not additional mice. Pool sizes were: male day 10, 19; male day 21-1, 7; male day 21-2, 6; male day 28, 8; female day 10, 18; female day 21-1, 6; female day 21-2, 3; female day 28-1, 6; female day 28-2, 6.

## Calculation of proteomic divergence
Proteomic divergence was assessed on the full QC proteome, with no significance pre-filter.

### Absolute difference (|Δlog2|)
For each protein, the magnitude of change was the absolute difference of group-mean log2 intensities:

```
|Δlog2| = | mean(log2 Intensity_Group1) − mean(log2 Intensity_Group2) |
```

Day-mean profiles were used. This is the y-axis of the raincloud plots. Points are proteins, not biological replicates. This metric has no direction.

### Signed log2 fold change
Signed log2 fold change used the same group means without the absolute value:

```
log2FC = mean(log2 Intensity_Group1) − mean(log2 Intensity_Group2)
```

Temporal contrasts are later day minus earlier day (positive = higher at the later day). Within-day contrasts are GFP+ minus GFP− (positive = higher in GFP+). Per-protein signed values are in Supplementary_Table_Signed_Log2FC.xlsx. Signed distributions were not plotted.

## Statistical analysis
The primary summary is the difference in median |Δlog2| between two distributions, with a percentile bootstrap 95% interval over proteins (B = 2000, seed 20261003). Proteins are correlated, so the interval is descriptive, not a mouse-level confidence interval.

Two-sided Kolmogorov–Smirnov tests were also run. Benjamini–Hochberg FDR was applied across the 12 planned contrasts. KS q-values are secondary: with n = 3404 they are sensitive to small shifts. Raw KS stars on the figures are not the primary evidence. Results are in Proteomic_Divergence_EffectSize_BH.csv.

GFP+ median |Δlog2| exceeded GFP− at every temporal window: day 21 vs 10, +0.302 (0.276–0.328); day 28 vs 10, +0.225 (0.202–0.242); day 28 vs 21, +0.183 (0.167–0.201). Within GFP+, the largest step was day 21 vs 10 (median 0.604), then day 28 vs 10 (0.475), then day 28 vs 21 (0.351). Within-day |GFP+ vs GFP−| was higher at day 21 (median 0.512) than at day 10 (0.351) or day 28 (0.346). Day 10 vs day 28 was null (+0.005; −0.017 to +0.022).

## Visualization
Raincloud plots (violin, boxplot, protein points) were made with ggplot2 and ggsignif. The y-axis is |Δlog2 fold change|, fixed at 0–7. Brackets on the current figures still show raw KS p-value stars.
