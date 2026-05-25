# Reproducibility Notes

## Environment

The final development environment used R 4.3.x and Python 3.10+ on Windows. Equivalent Linux/macOS environments should work if paths and fonts are adjusted.

Install Python dependencies:

```powershell
python -m pip install -r requirements-python.txt
```

Install R dependencies:

```r
install.packages(readLines("requirements-r.txt"))
```

Some survival plots use `survminer`, which may require system libraries depending on the R platform.

## Path Configuration

Set these environment variables before running:

- `EOBC_BIOMARKER_ROOT`: local biomarker root directory.
- `EOBC_FINAL_ANALYSIS_DIR`: final analysis directory; defaults conceptually to `<EOBC_BIOMARKER_ROOT>/final_analysis`.
- `EOBC_PROJECT_ROOT`: parent research root, required by scripts that access DepMap source folders outside the biomarker folder.
- `EOBC_VALIDATION_TABLES_DIR`: optional override for the curated validation table folder.

Do not commit `.env` or `config/paths.local.yml`.

## Suggested Run Order

1. `scripts/R/make_group_biomarker_landscape_r_v1.R`
2. `scripts/R/make_os_elasticnet_full_candidate_r_v1.R`
3. `scripts/R/make_full_layer_candidate_validation_r_v1.R`
4. `scripts/python/run_os_phase1_subset_screen.py`
5. `scripts/python/run_os_consensus_analysis_v3.py`
6. `scripts/python/run_biomarker_publication_suite_v2.py`
7. `scripts/R/make_group_to_outcome_bridge_r_v1.R`
8. `scripts/R/make_group_marker_domain_specific_r_v1.R`
9. `scripts/R/make_group_biomarker_contextual_domain_raw_meth_r_v1.R`
10. `scripts/R/make_group_biomarker_trajectory_figures_r_v1.R`
11. `scripts/R/make_eobc_biomarker_publication_final_figures_raw_meth_r_v1.R`
12. `scripts/R/make_final_individual_biomarker_evidence_plots_r_v1.R`
13. `scripts/R/make_meth_drugclass_biomarker_right_sankey.R`
14. `scripts/python/make_sankey_aligned_supplement_evidence_figures.py`

The convenience wrapper `scripts/workflow/run_publication_pipeline.ps1` executes the same high-level order.

## Final Figure Consistency Checks

After the final Sankey run, inspect:

- `metadata/meth_group_pathway_evidence_biomarker_right_sankey_final_audit_v12.csv`
- `metadata/supplement_vs_sankey_evidence_reconciliation.csv`

Expected final checks:

- RF METH mismatch: 0
- RF RNA mismatch: 0
- Meth-RNA label mismatch: 0
- EOBC group mismatch: 0
- Supplement/Sankey evidence mismatch: 0

