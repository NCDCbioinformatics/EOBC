# Methods Summary

## EOBC Biomarker Set

The analysis begins from a curated EOBC biomarker panel assigned to six biological families: immune, repair, glycolysis/TCA, fatty acid, kinase signaling, and hormone signaling. RNA expression and raw TSS methylation beta values are processed as matched tumor-level matrices and summarized as z-scores.

EOBC states are ordered as:

1. `G1 | H`: hypermethylated EOBC state
2. `G2 | I`: intermediate EOBC state
3. `G3 | L-like`: luminal-like EOBC state
4. `G4 | L`: luminal EOBC state

For group-level marker maps, each marker is assigned to the EOBC group with the highest group mean z-score. Methylation features retain raw TSS beta-value direction in the final paper figures.

## RNA-TSS Coupling

RNA/TSS coupling is computed as Spearman correlation between matched RNA expression and raw TSS methylation values. The final label classes are:

- `Meth-RNA inverse`: negative Spearman rho, supporting promoter/TSS-linked repression.
- `Meth-RNA positive`: positive Spearman rho.
- `Meth-RNA weak`: weak or non-classified coupling.

The compact heatmaps and Sankey biomarker labels use the same coupling source table to avoid figure-to-figure drift.

## RF Layer Scores

Random-forest layer scores are read from `gene_driver_summary_all_targets.csv`. For each gene, the final Sankey uses the maximum RF importance observed in the methylation layer and the maximum RF importance observed in the transcriptome layer. The displayed label records whether the marker is methylation-dominant, RNA-dominant, mixed, or weak by these layer-wise RF maxima.

The final Sankey audit table verifies that the plotted RF scores match the source table.

## OS Evidence

OS evidence is generated from RNA and methylation marker strata using Kaplan-Meier/log-rank workflows and consensus summaries. Final Sankey OS routes are restricted to KM-significant marker evidence:

- `OS | Protective`
- `OS | Adverse`

Individual OS plots are generated as supplementary Kaplan-Meier panels with risk tables.

## DepMap/GDSC Drug Evidence

DepMap/GDSC drug evidence is summarized at the marker-drug level using AUC/IC50 association statistics and dose-response curves. Final Sankey routing keeps curated therapeutic classes only and removes other/exploratory and experimental/natural drug classes from the terminal evidence routes.

Dose-response supplementary panels compare marker-high and marker-low cell-line strata across log10 drug dose. Dense dose grids are thinned for readability while preserving the response-curve shape used for AUC interpretation.

## Immune TIL/TMB Evidence

Immune evidence is computed as continuous marker expression versus TIL score and marker expression versus TMB score Spearman correlations. TIL and TMB evidence are harmonized into:

- `Immune | TIL/TMB positive`
- `Immune | TIL/TMB negative`
- `Immune | TIL/TMB weak`

TMB-specific correlations exclude the six zero/outlier TMB-log1p samples identified during final figure QC. TIL analyses keep the full matched immune cohort. The same TMB exclusion rule is applied to final Sankey routing and supplementary immune correlation panels.

## Final Sankey

The final Sankey is built from harmonized route evidence only:

- KM-significant OS routes.
- Curated DepMap/GDSC drug-response routes.
- Significant continuous TIL/TMB immune routes after TMB outlier exclusion.

The column order is:

1. EOBC methylation state
2. Biological program
3. Terminal biomarker
4. Functional evidence domain

The ribbon color is intentionally progressive: the first flow segment emphasizes EOBC state, the second biological program, the third biomarker coupling/RF-layer dominance, and the final segment functional evidence. This keeps colors aligned with the earlier figure palette while making the route interpretation visually explicit.

