# EOBC TSS-Methylation Biomarker Analysis

This repository contains the analysis code used to generate the EOBC raw TSS-methylation biomarker figures and supporting evidence panels for manuscript preparation.

The code links EOBC methylation states, RNA/TSS-methylation marker behavior, biological programs, OS evidence, curated DepMap/GDSC drug-response evidence, and TIL/TMB immune correlations. Raw patient-level and controlled-access input data are not included.

## Repository Contents

- `scripts/R/`: R scripts for biomarker landscape, contextual-domain figures, survival/drug/immune evidence panels, and final Sankey plots.
- `scripts/python/`: Python scripts for OS subset screening, consensus evidence summaries, domain evidence plots, and Sankey-aligned supplement panels.
- `scripts/workflow/`: Convenience workflow script for running the publication pipeline after local paths are configured.
- `config/paths_template.yml`: Template for local input/output paths.
- `metadata/`: Non-sensitive gene-level manifests/audits used to verify the final figure set.
- `docs/`: Methods, data availability, reproducibility notes, and manuscript code-availability text.

## Quick Start

1. Install R and Python dependencies listed in `requirements-r.txt` and `requirements-python.txt`.
2. Copy `.env.example` to `.env` or set the environment variables directly.
3. Point `EOBC_BIOMARKER_ROOT` to the local biomarker analysis root that contains the expected input and intermediate table folders.
4. Run the workflow:

```powershell
.\scripts\workflow\run_publication_pipeline.ps1 -BiomarkerRoot "D:\path\to\biomarker"
```

The scripts write results under `EOBC_FINAL_ANALYSIS_DIR`, defaulting to `<EOBC_BIOMARKER_ROOT>/final_analysis`.

## Key Final Analysis Choices

- Raw TSS methylation beta-value z-scores are used for the final EOBC methylation-state figures.
- RNA/TSS coupling is summarized as Spearman rho; negative rho supports promoter/TSS-linked repression.
- TMB correlations exclude the six zero/outlier TMB-log1p samples for TMB-only analyses; TIL analyses retain the full matched cohort.
- DepMap routes use curated therapeutic classes and exclude other/exploratory or experimental/natural classes from final Sankey evidence routing.
- RF layer scores are read from `gene_driver_summary_all_targets.csv` and summarized as gene-wise methylation/transcriptome RF maxima.
- The final Sankey uses harmonized evidence routes only: KM-significant OS, curated DepMap drug-response classes, and continuous TIL/TMB correlations.

## Verification Metadata

The final Sankey audit table in `metadata/meth_group_pathway_evidence_biomarker_right_sankey_final_audit_v12.csv` checks that:

- EOBC group labels match the source group-dominance table.
- Meth-RNA labels match the RNA/TSS coupling source.
- RF methylation and RNA values match the source RF table.
- Supplementary evidence panels agree with Sankey evidence routes.

The reconciliation table in `metadata/supplement_vs_sankey_evidence_reconciliation.csv` records the final supplement-versus-Sankey evidence match status.

## Data Availability

This repository intentionally excludes raw patient-level matrices, clinical metadata, controlled-access data, and bulky figure outputs. See `docs/data_availability.md` for expected input tables and external data sources.

## Citation

Please cite the associated manuscript when available. A provisional citation template is provided in `CITATION.cff`.

