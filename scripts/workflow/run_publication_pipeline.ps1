param(
    [Parameter(Mandatory = $true)]
    [string]$BiomarkerRoot,

    [string]$ProjectRoot = "",
    [string]$FinalAnalysisDir = "",
    [string]$RscriptPath = "",
    [string]$PythonPath = "python"
)

$ErrorActionPreference = "Stop"

if (-not $ProjectRoot) {
    $ProjectRoot = Split-Path -Parent $BiomarkerRoot
}

if (-not $FinalAnalysisDir) {
    $FinalAnalysisDir = Join-Path $BiomarkerRoot "final_analysis"
}

if (-not $RscriptPath) {
    if ($env:R_SCRIPT_PATH) {
        $RscriptPath = $env:R_SCRIPT_PATH
    } else {
        $RscriptPath = "Rscript"
    }
}

$env:EOBC_PROJECT_ROOT = $ProjectRoot
$env:EOBC_BIOMARKER_ROOT = $BiomarkerRoot
$env:EOBC_FINAL_ANALYSIS_DIR = $FinalAnalysisDir

Write-Host "EOBC_BIOMARKER_ROOT=$env:EOBC_BIOMARKER_ROOT"
Write-Host "EOBC_FINAL_ANALYSIS_DIR=$env:EOBC_FINAL_ANALYSIS_DIR"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Push-Location $repoRoot

try {
    & $RscriptPath "scripts/R/make_group_biomarker_landscape_r_v1.R"
    & $RscriptPath "scripts/R/make_os_elasticnet_full_candidate_r_v1.R"
    & $RscriptPath "scripts/R/make_full_layer_candidate_validation_r_v1.R"
    & $PythonPath "scripts/python/run_os_phase1_subset_screen.py"
    & $PythonPath "scripts/python/run_os_consensus_analysis_v3.py"
    & $PythonPath "scripts/python/run_biomarker_publication_suite_v2.py"
    & $RscriptPath "scripts/R/make_group_to_outcome_bridge_r_v1.R"
    & $RscriptPath "scripts/R/make_group_marker_domain_specific_r_v1.R"
    & $RscriptPath "scripts/R/make_group_biomarker_contextual_domain_raw_meth_r_v1.R"
    & $RscriptPath "scripts/R/make_group_biomarker_trajectory_figures_r_v1.R"
    & $RscriptPath "scripts/R/make_eobc_biomarker_publication_final_figures_raw_meth_r_v1.R"
    & $RscriptPath "scripts/R/make_final_individual_biomarker_evidence_plots_r_v1.R"
    & $RscriptPath "scripts/R/make_meth_drugclass_biomarker_right_sankey.R"
    & $PythonPath "scripts/python/make_sankey_aligned_supplement_evidence_figures.py"
}
finally {
    Pop-Location
}

