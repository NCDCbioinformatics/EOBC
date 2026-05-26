#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
})

env_path <- function(name, default = NULL, must_work = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) {
    return(normalizePath(value, winslash = "/", mustWork = must_work))
  }
  if (!is.null(default)) {
    return(normalizePath(default, winslash = "/", mustWork = must_work))
  }
  stop("Set environment variable ", name, " (see config/paths_template.yml).")
}

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
current_script <- if (length(script_file)) {
  normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("make_eobc_biomarker_publication_final_figures_raw_meth_r_v1.R", winslash = "/", mustWork = FALSE)
}
script_dir <- dirname(current_script)

biomarker_root <- env_path("EOBC_BIOMARKER_ROOT")
root_dir <- env_path("EOBC_FINAL_ANALYSIS_DIR", file.path(biomarker_root, "final_analysis"))
src_run_dir <- file.path(root_dir, "group_biomarker_contextual_domain_raw_meth_r_v1")
src_plot_dir <- file.path(src_run_dir, "plots")
src_table_dir <- file.path(src_run_dir, "tables")

final_dir <- file.path(root_dir, "eobc_biomarker_publication_final_figures_r_v1")
final_plot_dir <- file.path(final_dir, "plots")
final_table_dir <- file.path(final_dir, "tables")
final_source_dir <- file.path(final_dir, "source")
final_fig_dir <- file.path(root_dir, "final_fig")
final_revision_dir <- file.path(final_fig_dir, "final_paper_figures_20260522_v2")
final_plot_dirs <- unique(c(final_plot_dir, final_fig_dir, final_revision_dir))

invisible(lapply(final_plot_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(final_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_source_dir, recursive = TRUE, showWarnings = FALSE)

source_script <- file.path(script_dir, "make_group_biomarker_contextual_domain_raw_meth_r_v1.R")

message("Sourcing contextual-domain plotting script...")
source(source_script, local = FALSE)

save_final <- function(plot, basename, width, height, dpi = 450) {
  saved <- unlist(lapply(final_plot_dirs, function(dest_dir) {
    png_path <- file.path(dest_dir, paste0(basename, ".png"))
    pdf_path <- file.path(dest_dir, paste0(basename, ".pdf"))
    ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
    ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
    c(png_path, pdf_path)
  }))
  invisible(saved)
}

copy_panel <- function(src_basename, final_basename) {
  for (ext in c("png", "pdf")) {
    src <- file.path(src_plot_dir, paste0(src_basename, ".", ext))
    if (file.exists(src)) {
      for (dest_dir in final_plot_dirs) {
        dst <- file.path(dest_dir, paste0(final_basename, ".", ext))
        file.copy(src, dst, overwrite = TRUE)
      }
    } else {
      warning("Missing source panel: ", src)
    }
  }
}

publication_theme <- theme(
  legend.position = "bottom",
  legend.box = "vertical",
  legend.box.just = "center",
  legend.margin = margin(t = 3, r = 2, b = 3, l = 2),
  legend.box.margin = margin(t = 4, r = 2, b = 6, l = 2),
  plot.margin = margin(8, 12, 24, 12),
  plot.background = element_rect(fill = "white", color = NA)
)

state_guides <- guides(
  fill = guide_legend(
    nrow = 2, byrow = TRUE, order = 1,
    override.aes = list(shape = 21, color = "#111827", size = 3.0, alpha = 1, stroke = 0.7)
  ),
  color = guide_legend(
    nrow = 1, byrow = TRUE, order = 2,
    override.aes = list(shape = 21, fill = "white", size = 3.0, alpha = 1, stroke = 1.0)
  ),
  size = guide_legend(
    nrow = 1, byrow = TRUE, order = 3,
    override.aes = list(shape = 21, color = "#111827", fill = "white", alpha = 1)
  )
)

immune_state_guides <- guides(
  fill = guide_legend(
    nrow = 2, byrow = TRUE, order = 1,
    override.aes = list(shape = 21, color = "#111827", size = 3.0, alpha = 1, stroke = 0.7)
  ),
  color = guide_legend(
    nrow = 2, byrow = TRUE, order = 2,
    override.aes = list(shape = 21, fill = "white", size = 3.0, alpha = 1, stroke = 1.0)
  ),
  size = guide_legend(
    nrow = 1, byrow = TRUE, order = 3,
    override.aes = list(shape = 21, color = "#111827", fill = "white", alpha = 1)
  )
)

trajectory_guides <- guides(
  fill = guide_legend(nrow = 1, byrow = TRUE, order = 1),
  size = guide_legend(nrow = 1, byrow = TRUE, order = 2)
)

# Main-text panel A is the sample-level heatmap because it carries the full
# tumor-level EOBC/pseudotime ordering. The compact group-level heatmap is
# retained as a supplemental panel below.
copy_panel(
  "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1",
  "Figure_1A_EOBC_biomarker_activity_heatmap_with_RNA_TSS_coupling"
)

supplemental_source_panels <- c(
  "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1",
  "Figure_01A_EOBC_group_activity_coupling_heatmap_R_v3",
  "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v3",
  "Figure_02A_OS_state_resolved_group_marker_evidence_R_v3",
  "Figure_03A_DepMap_state_resolved_group_marker_evidence_R_v3",
  "Figure_04A_Immune_state_resolved_group_marker_evidence_R_v3",
  "Figure_05A_OS_annotated_family_trajectory_R_v2",
  "Figure_05B_DepMap_annotated_family_trajectory_R_v2",
  "Figure_05C_Immune_annotated_family_trajectory_R_v2"
)

invisible(lapply(supplemental_source_panels, function(panel) copy_panel(panel, panel)))

if (exists("p_cor_matrix_v3")) {
  save_final(
    p_cor_matrix_v3 +
      theme(
        plot.margin = margin(10, 16, 12, 10),
        legend.position = "right",
        legend.box = "vertical",
        plot.background = element_rect(fill = "white", color = NA)
      ),
    "Figure_1B_RNA_TSS_methylation_coupling_matrix",
    width = 9.4,
    height = 7.8
  )
}

if (exists("p_os_state_v2")) {
  save_final(
    p_os_state_v2 + publication_theme + state_guides,
    "Figure_1C_OS_evidence_projected_on_EOBC_group_markers",
    width = 15.8,
    height = 10.4
  )
}

if (exists("p_depmap_state_v2")) {
  save_final(
    p_depmap_state_v2 + publication_theme + state_guides,
    "Figure_1D_DepMap_evidence_projected_on_EOBC_group_markers",
    width = 15.8,
    height = 10.4
  )
}

if (exists("p_immune_state_v2")) {
  save_final(
    p_immune_state_v2 + publication_theme + immune_state_guides,
    "Figure_1E_Immune_evidence_projected_on_EOBC_group_markers",
    width = 15.8,
    height = 10.4
  )
}

if (exists("p_os_biomarker_summary_v1")) {
  save_final(
    p_os_biomarker_summary_v1,
    "Figure_S2A_OS_significant_biomarker_summary_no_group",
    width = 11.6,
    height = 6.5
  )
}

if (exists("p_depmap_biomarker_summary_v1")) {
  save_final(
    p_depmap_biomarker_summary_v1,
    "Figure_S2B_DepMap_significant_biomarker_summary_no_group",
    width = 12.2,
    height = 7.4
  )
}

if (exists("p_immune_biomarker_summary_v1")) {
  save_final(
    p_immune_biomarker_summary_v1,
    "Figure_S2C_Immune_significant_biomarker_summary_no_group",
    width = 12.2,
    height = 8.0
  )
}

trajectory_theme <- theme(
  legend.position = "bottom",
  legend.box = "vertical",
  legend.box.just = "center",
  legend.box.margin = margin(t = 4, r = 2, b = 8, l = 2),
  legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
  plot.margin = margin(10, 12, 24, 12),
  plot.background = element_rect(fill = "white", color = NA)
)

if (exists("p_traj_os_v2")) {
  save_final(
    p_traj_os_v2 + trajectory_theme + trajectory_guides,
    "Figure_1F_OS_evidence_on_EOBC_family_trajectories",
    width = 14.2,
    height = 7.4
  )
}

if (exists("p_traj_depmap_v2")) {
  save_final(
    p_traj_depmap_v2 + trajectory_theme + trajectory_guides,
    "Figure_1G_DepMap_evidence_on_EOBC_family_trajectories",
    width = 14.2,
    height = 7.4
  )
}

if (exists("p_traj_immune_v2")) {
  save_final(
    p_traj_immune_v2 + trajectory_theme + trajectory_guides,
    "Figure_1H_Immune_evidence_on_EOBC_family_trajectories",
    width = 14.2,
    height = 7.4
  )
}

if (exists("fig_group_activity_heat_v3")) {
  save_final(
    fig_group_activity_heat_v3 +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = margin(8, 10, 12, 10),
        plot.background = element_rect(fill = "white", color = NA)
      ),
    "Figure_S1A_EOBC_group_mean_activity_coupling_heatmap",
    width = 12.8,
    height = 10.8
  )
}

panel_manifest <- tribble(
  ~final_panel, ~final_basename, ~source_or_object, ~recommended_use, ~interpretation_note,
  "1A", "Figure_1A_EOBC_biomarker_activity_heatmap_with_RNA_TSS_coupling", "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1", "Main figure", "Tumors are ordered by EOBC group/pseudotime; methylation rows are raw TSS beta-value z-scores.",
  "1B", "Figure_1B_RNA_TSS_methylation_coupling_matrix", "p_cor_matrix_v3", "Main figure", "Raw RNA versus raw TSS methylation Spearman rho; negative diagonal correlations support promoter/TSS-linked repression.",
  "1C", "Figure_1C_OS_evidence_projected_on_EOBC_group_markers", "p_os_state_v2", "Main figure", "Nominal Cox regression p < 0.05 OS-linked markers are projected onto the group-defining contrast/Kruskal-Wallis coordinate system; point size encodes best-endpoint Cox regression p-value.",
  "1D", "Figure_1D_DepMap_evidence_projected_on_EOBC_group_markers", "p_depmap_state_v2", "Main figure", "DepMap-linked markers are projected onto the group-defining contrast/Kruskal-Wallis coordinate system; point size encodes the AUC/IC50 dose-response curve p-value.",
  "1E", "Figure_1E_Immune_evidence_projected_on_EOBC_group_markers", "p_immune_state_v2", "Main figure", "TIL/TMB immune-linked markers are projected onto the group-defining contrast/Kruskal-Wallis coordinate system; rows require either TIL or TMB Spearman p < 0.05 plus concordant non-zero TIL/TMB rho directions, with six zero-TMB outlier samples excluded from TMB correlations.",
  "S2A", "Figure_S2A_OS_significant_biomarker_summary_no_group", "p_os_biomarker_summary_v1", "Supplement", "Domain-only nominal Cox regression p < 0.05 OS biomarker checklist without EOBC G1-G4 paneling.",
  "S2B", "Figure_S2B_DepMap_significant_biomarker_summary_no_group", "p_depmap_biomarker_summary_v1", "Supplement", "Domain-only DepMap AUC/IC50 dose-response curve p < 0.05 biomarker checklist without EOBC G1-G4 paneling.",
  "S2C", "Figure_S2C_Immune_significant_biomarker_summary_no_group", "p_immune_biomarker_summary_v1", "Supplement", "Domain-only immune checklist requiring either TIL or TMB Spearman p < 0.05 plus concordant non-zero TIL/TMB rho directions; six zero-TMB outlier samples are excluded from TMB correlations.",
  "1F", "Figure_1F_OS_evidence_on_EOBC_family_trajectories", "p_traj_os_v2", "Main or supporting figure", "Family-level EOBC trajectories annotated with counts/direction of OS-linked markers.",
  "1G", "Figure_1G_DepMap_evidence_on_EOBC_family_trajectories", "p_traj_depmap_v2", "Main or supporting figure", "Family-level EOBC trajectories annotated with drug-response-linked marker counts.",
  "1H", "Figure_1H_Immune_evidence_on_EOBC_family_trajectories", "p_traj_immune_v2", "Main or supporting figure", "Family-level EOBC trajectories annotated with immune-linked marker counts.",
  "S1A", "Figure_S1A_EOBC_group_mean_activity_coupling_heatmap", "fig_group_activity_heat_v3", "Supplement", "Compact group-level version of the activity/coupling heatmap."
)

write_csv(panel_manifest, file.path(final_dir, "publication_final_figure_manifest.csv"))
write_csv(panel_manifest, file.path(final_fig_dir, "publication_final_figure_manifest.csv"))
write_csv(panel_manifest, file.path(final_revision_dir, "publication_final_figure_manifest.csv"))

if (dir.exists(src_table_dir)) {
  table_files <- list.files(src_table_dir, pattern = "\\.csv$", full.names = TRUE)
  file.copy(table_files, final_table_dir, overwrite = TRUE)
}

file.copy(source_script, file.path(final_source_dir, basename(source_script)), overwrite = TRUE)
file.copy(current_script,
          file.path(final_source_dir, "make_eobc_biomarker_publication_final_figures_raw_meth_r_v1.R"),
          overwrite = TRUE)

readme <- c(
  "# EOBC biomarker publication final figure set",
  "",
  "This folder contains the recommended publication-ready contextual EOBC biomarker figure set.",
  "",
  "## Figure order",
  "",
  "1. Figure 1A: tumor-level RNA expression/raw TSS methylation heatmap with RNA-TSS coupling and domain evidence strips.",
  "2. Figure 1B: raw RNA versus raw TSS methylation coupling matrix.",
  "3. Figure 1C: OS evidence plotted for nominal Cox regression p < 0.05 best-endpoint hits within EOBC group-defining marker programs; point size encodes Cox regression p-value.",
  "4. Figure 1D: DepMap drug-response evidence plotted for AUC/IC50 dose-response curve p < 0.05 hits within EOBC group-defining marker programs; point size encodes the same curve p-value.",
  "5. Figure 1E: immune evidence requires either TIL or TMB Spearman p < 0.05 plus concordant non-zero TIL/TMB rho directions within EOBC group-defining marker programs; point size encodes the best of the two p-values and TMB correlations exclude six zero-TMB outlier samples.",
  "6. Figure S2A-C: domain-only significant biomarker summaries for OS, DepMap, and immune evidence without EOBC G1-G4 paneling.",
  "7. Figure 1F-H: family-trajectory summaries annotated by OS, DepMap, and immune evidence.",
  "",
  "## Interpretation notes",
  "",
  "- EOBC group order is fixed as G1 | H, G2 | I, G3 | L-like, G4 | L.",
  "- TSS/promoter methylation panels use raw beta-value z-scores; high methylation is shown directly without direction flipping.",
  "- In heatmap coupling strips, sample rho is the matched-tumor Spearman correlation; * marks nominal Spearman p < 0.05.",
  "- RNA-TSS correlation panels retain raw methylation values; negative Spearman rho is the expected direction for promoter/TSS methylation-linked transcriptional repression.",
  "- DepMap evidence uses the available gene-level AUC/IC50 dose-response curve statistics; exploratory drugs are retained when they are the prioritized biomarker-drug curve hit.",
  "- Immune evidence is based on continuous TIL score and TMB log1p Spearman association tests; no externally assigned immune hot/cold subtype is used.",
  "- Both PNG and editable PDF files are exported for each final panel.",
  "",
  "## Reproducibility",
  "",
  "- The source R scripts are copied into the source/ folder.",
  "- Supporting CSV tables are copied into the tables/ folder.",
  "- Panel-level provenance is recorded in publication_final_figure_manifest.csv."
)
writeLines(readme, file.path(final_dir, "README_publication_final_figure_set.md"))
writeLines(readme, file.path(final_fig_dir, "README_publication_final_figure_set.md"))
writeLines(readme, file.path(final_revision_dir, "README_publication_final_figure_set.md"))

message("Final publication figure set saved to: ", final_dir)
message("Final figure copies saved to: ", final_fig_dir)
message("Final revision copies saved to: ", final_revision_dir)
