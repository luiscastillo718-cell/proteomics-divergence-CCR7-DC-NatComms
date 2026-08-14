# Methods – Proteomic Divergence Analysis

## Proteomic data processing
Label-free quantification proteomics data were filtered to remove potential contaminants. Proteins with missing values across samples were excluded. Log\u2082-transformed intensities were median-centered across all samples to correct for technical variation.

## Calculation of proteomic divergence
Proteomic divergence was assessed at the distribution level without pre-filtering on statistical significance. Two complementary metrics were used:

### Absolute difference (|\Delta log\u2082|)
For each protein, the magnitude of change between two groups was calculated as the absolute difference between their mean log\u2082 intensities:

```
|\Delta log\u2082| = | mean(log\u2082 Intensity_Group1) − mean(log\u2082 Intensity_Group2) |
```

This metric was used for the primary analyses (Figure 3). It quantifies the overall scale of proteomic remodeling without regard to direction.

### Signed Log\u2082 fold change (Log2FC)
Signed Log2 fold changes were calculated as:

```
Log2FC = mean(log\u2082 Intensity_Group1) − mean(log\u2082 Intensity_Group2)
```

Distributions of signed Log2FC values were generated as supplementary analyses (Supplementary Figures S1–S4) to more directly address the reviewer’s suggestion.

All divergence analyses were performed on the full QC-filtered proteome (n = 3404 proteins). Day-mean profiles were used to compute both metrics. Individual biological replicates are shown as separate points in the raincloud plots.

## Statistical analysis
Differences between distributions of proteomic divergence were evaluated using two-sided Kolmogorov-Smirnov tests. All pairwise comparisons of interest were performed. Because the number of tests was moderate and the analyses were exploratory, *P*-values were adjusted for multiple testing using the Benjamini-Hochberg procedure to control the false discovery rate (FDR). Raw and FDR-adjusted *P*-values are reported in Supplementary Table X.

## Visualization
Distributions of proteomic changes were visualized as raincloud plots comprising violin plots (showing the full density), overlaid boxplots (indicating median and interquartile range), and individual data points (representing per-protein values). Plots were generated in R using `ggplot2` and `ggsignif`. Significance brackets on figures indicate raw *P*-values from Kolmogorov-Smirnov tests.
