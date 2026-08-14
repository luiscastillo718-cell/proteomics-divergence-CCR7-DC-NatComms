# Response to Reviewer #2 – Proteomic Divergence Analysis

We thank the reviewer for raising this important point. We agree that calculating a single median value after significance-based pre-filtering can introduce bias due to differential statistical power across conditions and timepoints. We have completely revised this analysis to address the concern.

In the revised manuscript, proteomic divergence is now evaluated at the **distribution level** using the **complete QC-filtered proteome** (n = 3404 proteins) **without any pre-filtering** on statistical significance. For each comparison, we calculated the absolute difference in log\u2082 intensity (|\Delta log\u2082|) between group means for every protein and compared the resulting distributions using Kolmogorov-Smirnov tests. This approach avoids the statistical power bias highlighted by the reviewer while focusing on the overall magnitude of proteomic remodeling.

To more directly follow the reviewer’s suggestion to analyze Log2FC data at the distribution level, we have also included distributions of **signed Log2 fold changes** as supplementary material (Supplementary Figures S1–S4). These analyses support the conclusions drawn from the primary |\Delta log\u2082| approach.

All pairwise comparisons were adjusted for multiple testing using the Benjamini-Hochberg method to control the false discovery rate (FDR). Raw and FDR-adjusted *P*-values are provided in Supplementary Table X.

Replicate structure is transparently represented: day-mean profiles were used to generate the distributions, while individual biological replicates are shown as separate points in the raincloud plots (Figure 3 and related panels).

These revisions are presented in **Figure 3** (including the new panel comparing temporal changes between GFP+ and GFP–) and **Supplementary Figures S1–S4**, with full statistical details in Supplementary Table X. We believe this revised framework is more transparent, statistically robust, and directly responsive to the reviewer’s concern.
