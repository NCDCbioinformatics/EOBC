from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
PUB_SCRIPT = SCRIPT_DIR / "run_biomarker_publication_suite_v2.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pub = load_module(PUB_SCRIPT, "publication_suite_v2")
phase1 = pub.phase1

BIOMARKER_ROOT = pub.BIOMARKER_ROOT
FINAL_ANALYSIS_ROOT = pub.FINAL_ANALYSIS_ROOT
FINAL_FIG_DIR = FINAL_ANALYSIS_ROOT / "final_fig"
FINAL_REVISION_DIR = FINAL_FIG_DIR / "final_paper_figures_20260522_v2"
OUT_ROOTS = [
    FINAL_FIG_DIR / "individual_biomarker_evidence_20260522",
    FINAL_REVISION_DIR / "individual_biomarker_evidence_20260522",
]

SANKY_TABLE = (
    FINAL_ANALYSIS_ROOT
    / "eobc_integrated_story_sankey_r_v5"
    / "tables"
    / "meth_group_pathway_evidence_biomarker_right_sankey_routes_full.csv"
)
OS_TABLE = (
    FINAL_ANALYSIS_ROOT
    / "group_marker_domain_specific_r_v1"
    / "tables"
    / "OS_group_defined_marker_endpoint_results.csv"
)
DEP_PROBE_TABLE = pub.VALIDATION_TABLES / "depmap_methylation_probe_selection.csv"
DEP_METH_TABLE = (
    BIOMARKER_ROOT.parent
    / "depmap"
    / "meth"
    / "Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv"
)
DEP_DOSE_TABLE = (
    BIOMARKER_ROOT.parent
    / "depmap"
    / "drug"
    / "Drug_sensitivity_replicate-level_dose_(Sanger_GDSC2)_subsetted_NAsdropped.csv"
)
EXPR_TABLE = (
    FINAL_ANALYSIS_ROOT
    / "group_biomarker_landscape_r_v1"
    / "tables"
    / "group_biomarker_long_values.csv"
)
IMMUNE_TABLE = (
    BIOMARKER_ROOT
    / "10_external_validation_TIL_TMB_DepMap_final_v7"
    / "tables"
    / "immune_signature_tcga_til_tmb_merged.csv"
)


TEXT = "#111827"
MUTED = "#64748B"
GRID = "#E2E8F0"
BLUE = "#2563EB"
ORANGE = "#D95F02"

FAMILY_COLORS = {
    "Immune": "#4EA5F0",
    "Repair": "#47C56B",
    "Glycolysis / TCA": "#E6BC18",
    "Fatty acid": "#F4A259",
    "Kinase signaling": "#9B6AE8",
    "Hormone signaling": "#9AAABC",
}
IMMUNE_CLASS_COLORS = {
    "Immune | TIL/TMB positive": "#0F9F9A",
    "Immune | TIL/TMB weak": "#94A3B8",
    "Immune | TIL/TMB negative": "#475569",
}


def ensure_dirs() -> None:
    for root in OUT_ROOTS:
        (root / "supplement_onepage").mkdir(parents=True, exist_ok=True)
        (root / "tables").mkdir(parents=True, exist_ok=True)


def save_to_all(fig: plt.Figure, stub: str) -> None:
    for root in OUT_ROOTS:
        out = root / "supplement_onepage" / stub
        fig.savefig(out.with_suffix(".png"), dpi=450, bbox_inches="tight", facecolor="white")
        fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")


def write_table_all(df: pd.DataFrame, name: str) -> None:
    for root in OUT_ROOTS:
        df.to_csv(root / "tables" / name, index=False, encoding="utf-8-sig")


def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    vals = np.asarray([np.nan if pd.isna(v) else float(v) for v in pvals], dtype=float)
    q = np.full(vals.shape, np.nan, dtype=float)
    ok = np.isfinite(vals)
    if not ok.any():
        return q
    p = vals[ok]
    order = np.argsort(p)
    ranked = p[order]
    m = ranked.size
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(p)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def fmt_p(x: float) -> str:
    if x is None or not np.isfinite(float(x)):
        return "NA"
    x = float(x)
    if x < 0.001:
        return "<0.001"
    return f"{x:.3f}"


def fmt_rho(x: float) -> str:
    if x is None or not np.isfinite(float(x)):
        return "NA"
    return f"{float(x):.2f}"


def clean_evidence_label(x: str) -> str:
    return str(x).replace("\n", " ").replace("  ", " ").strip()


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def apply_axis_style(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#1F2937")
    ax.spines["bottom"].set_color("#1F2937")
    ax.tick_params(colors="#475569", labelsize=7.5)
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.set_axisbelow(True)


def plot_os_km(route: pd.DataFrame) -> pd.DataFrame:
    os_routes = (
        route.loc[route["route_type"].eq("OS")]
        .copy()
        .sort_values(["evidence_node", "gene_clean"])
    )
    genes = os_routes["gene_clean"].drop_duplicates().tolist()
    os_meta = pd.read_csv(OS_TABLE)
    selected = (
        os_meta.loc[
            os_meta["gene"].isin(genes)
            & os_meta["modality"].eq("METH")
            & os_meta["selected"].astype(bool)
        ]
        .sort_values(["gene", "p", "logrank_p"], na_position="last")
        .groupby("gene", as_index=False)
        .head(1)
    )
    if selected.empty:
        raise RuntimeError("No selected METH OS rows found for Sankey OS genes")

    clinical = phase1.load_clinical(phase1.CLINICAL_PATH)
    meth = phase1.zscore_df(phase1.load_omics_matrix(phase1.METH_PATH, genes))
    route_by_gene = os_routes.set_index("gene_clean").to_dict("index")

    ncols = 3
    nrows = math.ceil(len(selected) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.8, 4.15 * nrows), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    manifest_rows: List[Dict[str, object]] = []

    for ax, (_, row) in zip(axes, selected.iterrows()):
        gene = str(row["gene"])
        endpoint = str(row["endpoint"])
        endpoint_df = phase1.build_endpoint_df(clinical, endpoint)
        endpoint_df = endpoint_df[endpoint_df["Sample"].isin(meth.index)].copy()
        endpoint_df["score"] = endpoint_df["Sample"].map(meth[gene])
        endpoint_df["time_months"] = endpoint_df["time_days"] / phase1.DAYS_PER_MONTH
        split = float(endpoint_df["score"].median())
        endpoint_df["group"] = np.where(endpoint_df["score"] >= split, "High TSS-METH", "Low TSS-METH")
        endpoint_df = endpoint_df.dropna(subset=["time_months", "event", "group"]).copy()

        p_logrank = np.nan
        try:
            hi = endpoint_df[endpoint_df["group"].eq("High TSS-METH")]
            lo = endpoint_df[endpoint_df["group"].eq("Low TSS-METH")]
            p_logrank = float(
                logrank_test(
                    hi["time_months"],
                    lo["time_months"],
                    event_observed_A=hi["event"],
                    event_observed_B=lo["event"],
                ).p_value
            )
        except Exception:
            pass

        for group_name, color in [("Low TSS-METH", BLUE), ("High TSS-METH", ORANGE)]:
            sub = endpoint_df[endpoint_df["group"].eq(group_name)]
            kmf = KaplanMeierFitter(label=group_name)
            kmf.fit(sub["time_months"], sub["event"])
            surv = kmf.survival_function_.reset_index()
            ax.step(
                surv.iloc[:, 0],
                surv.iloc[:, 1],
                where="post",
                color=color,
                linewidth=2.0,
                label=group_name,
            )

        route_row = route_by_gene.get(gene, {})
        family = str(route_row.get("Family6", ""))
        border = FAMILY_COLORS.get(family, "#94A3B8")
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
        ax.spines["left"].set_color(border)
        ax.spines["bottom"].set_color(border)
        apply_axis_style(ax)
        ax.set_ylim(0, 1.03)
        ax.set_xlim(left=0)
        ax.set_title(
            f"{gene} | {clean_evidence_label(route_row.get('evidence_node', 'OS'))}",
            loc="left",
            fontsize=10.2,
            fontweight="bold",
            color=TEXT,
            pad=7,
        )
        ax.text(
            0.02,
            0.08,
            f"{row['endpoint_label']}  p={fmt_p(p_logrank)}\nHR={fmt_rho(float(row['hr']))}; selected in OS model",
            transform=ax.transAxes,
            fontsize=7.6,
            color=MUTED,
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#E5E7EB", "alpha": 0.9},
        )
        manifest_rows.append(
            {
                "domain": "OS",
                "gene": gene,
                "endpoint": endpoint,
                "endpoint_label": row["endpoint_label"],
                "route_evidence": route_row.get("evidence_node", ""),
                "os_model_selected": True,
                "individual_logrank_p": p_logrank,
                "os_table_p": row["p"],
                "os_table_q": row["q"],
                "hr": row["hr"],
                "n": int(endpoint_df.shape[0]),
                "n_high": int(endpoint_df["group"].eq("High TSS-METH").sum()),
                "n_low": int(endpoint_df["group"].eq("Low TSS-METH").sum()),
            }
        )

    for ax in axes[len(selected) :]:
        ax.axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=9.2)
    fig.suptitle(
        "Sankey-aligned OS Kaplan-Meier evidence for terminal biomarkers",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.02,
        0.958,
        "Panels are restricted to OS routes shown in the final Sankey. METH-high/low groups are median splits of raw TSS-methylation z-score.",
        ha="left",
        fontsize=10,
        color=MUTED,
    )
    fig.text(0.5, 0.04, "Months", ha="center", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.text(0.026, 0.52, "Survival probability", va="center", rotation="vertical", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.115, top=0.89, wspace=0.24, hspace=0.38)
    save_to_all(fig, "Figure_S_OS_KM_sankey_aligned_biomarkers")
    plt.close(fig)
    return pd.DataFrame(manifest_rows)


def parse_drug_from_route(top_drugs: str) -> str:
    text = str(top_drugs)
    text = text.split(" (")[0]
    return text.strip()


def selected_dose_columns(drugs: Iterable[str]) -> Tuple[List[str], pd.DataFrame]:
    header = pd.read_csv(DEP_DOSE_TABLE, nrows=0)
    drug_set = set(drugs)
    rows = []
    columns = []
    pat = re.compile(r"^(?P<drug>.+? \(GDSC2:\d+\)) (?P<dose>[0-9.eE+-]+).+ rep(?P<rep>\d+)$")
    for col in header.columns:
        match = pat.match(col)
        if not match:
            continue
        drug_raw = match.group("drug")
        drug_clean = drug_raw.split(" (GDSC2:")[0].strip()
        if drug_clean not in drug_set:
            continue
        rows.append(
            {
                "column": col,
                "drug_raw": drug_raw,
                "drug_clean": drug_clean,
                "dose_uM": float(match.group("dose")),
                "rep": int(match.group("rep")),
            }
        )
        columns.append(col)
    return columns, pd.DataFrame(rows)


def thin_curve(df: pd.DataFrame, max_points: int = 8) -> pd.DataFrame:
    out = []
    for (group, drug), sub in df.groupby(["marker_group", "drug_clean"], sort=False):
        sub = sub.sort_values("log10_uM").copy()
        unique = sub["log10_uM"].nunique()
        if unique <= max_points:
            out.append(sub)
            continue
        dose_order = (
            sub[["log10_uM"]]
            .drop_duplicates()
            .sort_values("log10_uM")
            .assign(bin=lambda x: pd.qcut(np.arange(x.shape[0]), q=max_points, labels=False, duplicates="drop"))
        )
        binned = sub.merge(dose_order, on="log10_uM", how="left")
        agg = (
            binned.groupby("bin", as_index=False)
            .agg(
                dose_uM=("dose_uM", "mean"),
                log10_uM=("log10_uM", "mean"),
                mean_viability=("mean_viability", "mean"),
                se_viability=("se_viability", "mean"),
                n=("n", "sum"),
            )
            .assign(marker_group=group, drug_clean=drug)
        )
        out.append(agg)
    return pd.concat(out, ignore_index=True) if out else df


def plot_depmap(route: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    dep_routes = (
        route.loc[route["route_type"].eq("DepMap")]
        .copy()
        .sort_values(["gene_clean"])
        .reset_index(drop=True)
    )
    dep_routes["drug_clean"] = dep_routes["top_drugs"].map(parse_drug_from_route)
    genes = dep_routes["gene_clean"].drop_duplicates().tolist()
    drugs = dep_routes["drug_clean"].drop_duplicates().tolist()

    probe = pd.read_csv(DEP_PROBE_TABLE)
    selected_features = probe.loc[probe["gene"].isin(genes), ["gene", "selected_feature"]]
    feature_map = selected_features.set_index("gene")["selected_feature"].to_dict()
    meth_cols = ["depmap_id", "cell_line_display_name"] + sorted(set(feature_map.values()))
    meth = pd.read_csv(DEP_METH_TABLE, usecols=lambda c: c in set(meth_cols))

    dose_cols, dose_meta = selected_dose_columns(drugs)
    if dose_meta.empty:
        raise RuntimeError("No replicate-level dose columns matched Sankey DepMap drugs")
    dose = pd.read_csv(DEP_DOSE_TABLE, usecols=["depmap_id"] + dose_cols)

    manifest_rows: List[Dict[str, object]] = []
    curve_rows: List[pd.DataFrame] = []
    ncols = 3
    nrows = math.ceil(dep_routes.shape[0] / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.8, 4.15 * nrows), sharey=True)
    axes = np.atleast_1d(axes).ravel()

    long_cache = {}
    for drug in drugs:
        cols = dose_meta.loc[dose_meta["drug_clean"].eq(drug), "column"].tolist()
        dm = dose[["depmap_id"] + cols].melt("depmap_id", var_name="column", value_name="viability")
        dm = dm.merge(dose_meta, on="column", how="left")
        dm["viability"] = pd.to_numeric(dm["viability"], errors="coerce")
        dm = dm.dropna(subset=["viability", "dose_uM"]).copy()
        dm["log10_uM"] = np.log10(dm["dose_uM"].astype(float))
        long_cache[drug] = dm

    for ax, (_, row) in zip(axes, dep_routes.iterrows()):
        gene = str(row["gene_clean"])
        drug = str(row["drug_clean"])
        feature = feature_map.get(gene)
        if feature is None or feature not in meth.columns:
            ax.axis("off")
            continue
        features = meth[["depmap_id", feature]].copy()
        features[feature] = pd.to_numeric(features[feature], errors="coerce")
        dlong = long_cache[drug].merge(features, on="depmap_id", how="inner")
        dlong = dlong.dropna(subset=[feature, "viability"]).copy()
        split = float(dlong.drop_duplicates("depmap_id")[feature].median())
        dlong["marker_group"] = np.where(dlong[feature] >= split, "High TSS-METH", "Low TSS-METH")
        curve = (
            dlong.groupby(["drug_clean", "marker_group", "dose_uM", "log10_uM"], as_index=False)
            .agg(
                mean_viability=("viability", "mean"),
                se_viability=("viability", lambda x: float(x.std(ddof=1) / math.sqrt(max(x.notna().sum(), 1))) if x.notna().sum() > 1 else 0.0),
                n=("viability", "count"),
            )
            .sort_values(["marker_group", "log10_uM"])
        )
        curve = thin_curve(curve, max_points=8)
        curve["gene"] = gene
        curve["selected_feature"] = feature
        curve_rows.append(curve)

        for group_name, color in [("Low TSS-METH", BLUE), ("High TSS-METH", ORANGE)]:
            sub = curve[curve["marker_group"].eq(group_name)].sort_values("log10_uM")
            if sub.empty:
                continue
            x = sub["log10_uM"].to_numpy(dtype=float)
            y = sub["mean_viability"].to_numpy(dtype=float)
            se = sub["se_viability"].fillna(0).to_numpy(dtype=float)
            ax.fill_between(x, np.maximum(0, y - se), np.minimum(1.12, y + se), color=color, alpha=0.15, linewidth=0)
            ax.plot(x, y, color=color, linewidth=2.0, marker="o", markersize=4.5, markeredgecolor="white", markeredgewidth=0.6, label=group_name)

        family = str(row["Family6"])
        border = FAMILY_COLORS.get(family, "#94A3B8")
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
        ax.spines["left"].set_color(border)
        ax.spines["bottom"].set_color(border)
        apply_axis_style(ax)
        ax.set_ylim(0, 1.1)
        ax.set_title(f"{gene} | {drug}", loc="left", fontsize=10.2, fontweight="bold", color=TEXT, pad=7)
        ax.text(
            0.02,
            0.08,
            f"{clean_evidence_label(row['evidence_node'])}\nFDR={fmt_p(float(row['best_fdr']))}; feature={feature}",
            transform=ax.transAxes,
            fontsize=7.4,
            color=MUTED,
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#E5E7EB", "alpha": 0.88},
        )
        manifest_rows.append(
            {
                "domain": "DepMap",
                "gene": gene,
                "drug_clean": drug,
                "drug_class": clean_evidence_label(str(row["evidence_node"]).replace("Drug |", "")),
                "route_evidence": row["evidence_node"],
                "best_fdr": row["best_fdr"],
                "selected_feature": feature,
                "n_curve_points_after_thinning": int(curve.shape[0]),
                "n_cell_line_dose_values": int(dlong.shape[0]),
            }
        )

    for ax in axes[dep_routes.shape[0] :]:
        ax.axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=9.2)
    fig.suptitle(
        "Sankey-aligned DepMap/GDSC dose-response evidence",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.02,
        0.958,
        "Panels are restricted to DepMap routes shown in the final Sankey. Replicate dose points are averaged and thinned to avoid over-jagged curves.",
        ha="left",
        fontsize=10,
        color=MUTED,
    )
    fig.text(0.5, 0.04, "Drug dose (log10 uM)", ha="center", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.text(0.026, 0.52, "Mean viability", va="center", rotation="vertical", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.115, top=0.89, wspace=0.24, hspace=0.38)
    save_to_all(fig, "Figure_S_DepMap_dose_response_sankey_aligned_biomarkers")
    plt.close(fig)
    return pd.DataFrame(manifest_rows), pd.concat(curve_rows, ignore_index=True)


def linear_fit_line(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.array([]), np.array([])
    coef = np.polyfit(x[ok], y[ok], deg=1)
    xs = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 60)
    ys = coef[0] * xs + coef[1]
    return xs, ys


def plot_immune(route: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    imm_routes = (
        route.loc[route["route_type"].eq("Immune")]
        .copy()
        .sort_values(["evidence_node", "gene_clean"])
        .reset_index(drop=True)
    )
    genes = imm_routes["gene_clean"].drop_duplicates().tolist()
    route_by_gene = imm_routes.set_index("gene_clean").to_dict("index")

    expr = pd.read_csv(EXPR_TABLE)
    expr = (
        expr.loc[expr["modality"].eq("RNA") & expr["gene"].isin(genes)]
        .copy()
        .rename(columns={"value_z": "rna_z"})
    )
    expr["rna_z"] = pd.to_numeric(expr["rna_z"], errors="coerce")
    immune = pd.read_csv(IMMUNE_TABLE)
    immune = immune[["Sample", "TIL_score", "tmb_log1p"]].copy()
    immune["TIL_score"] = pd.to_numeric(immune["TIL_score"], errors="coerce")
    immune["tmb_log1p"] = pd.to_numeric(immune["tmb_log1p"], errors="coerce")
    tmb_outliers = sorted(immune.loc[immune["tmb_log1p"].le(0) & immune["tmb_log1p"].notna(), "Sample"].unique().tolist())

    til = (
        expr.merge(immune[["Sample", "TIL_score"]], on="Sample", how="inner")
        .assign(metric="TIL", metric_raw=lambda x: x["TIL_score"])
        .dropna(subset=["rna_z", "metric_raw"])
    )
    tmb = (
        expr.merge(immune[["Sample", "tmb_log1p"]], on="Sample", how="inner")
        .loc[lambda x: ~x["Sample"].isin(tmb_outliers)]
        .assign(metric="TMB", metric_raw=lambda x: x["tmb_log1p"])
        .dropna(subset=["rna_z", "metric_raw"])
    )
    plot_df = pd.concat([til, tmb], ignore_index=True)
    plot_df["metric_z"] = plot_df.groupby("metric")["metric_raw"].transform(zscore)

    stats_rows = []
    for (gene, metric), sub in plot_df.groupby(["gene", "metric"], sort=False):
        rho, p = spearmanr(sub["rna_z"], sub["metric_raw"], nan_policy="omit")
        stats_rows.append(
            {
                "domain": "Immune",
                "gene": gene,
                "metric": metric,
                "n": int(sub[["rna_z", "metric_raw"]].dropna().shape[0]),
                "rho": float(rho) if np.isfinite(rho) else np.nan,
                "p": float(p) if np.isfinite(p) else np.nan,
                "route_evidence": route_by_gene.get(gene, {}).get("evidence_node", ""),
            }
        )
    stats = pd.DataFrame(stats_rows)
    stats["q"] = stats.groupby("metric")["p"].transform(lambda s: bh_fdr(s.to_numpy()))
    for idx, row in stats.iterrows():
        route_row = route_by_gene.get(row["gene"], {})
        metric = str(row["metric"])
        for field in ["rho", "p", "q", "n"]:
            route_key = f"{metric}_{field}"
            if route_key in route_row and pd.notna(route_row.get(route_key)):
                stats.at[idx, field] = route_row.get(route_key)

    ncols = 4
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15.2, 3.25 * nrows), sharex=False, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    metric_colors = {"TIL": "#16A34A", "TMB": "#7C3AED"}

    for ax, gene in zip(axes, genes):
        sub_gene = plot_df.loc[plot_df["gene"].eq(gene)].copy()
        route_row = route_by_gene.get(gene, {})
        evidence = str(route_row.get("evidence_node", ""))
        for metric in ["TIL", "TMB"]:
            sub = sub_gene.loc[sub_gene["metric"].eq(metric)]
            ax.scatter(
                sub["rna_z"],
                sub["metric_z"],
                s=21,
                color=metric_colors[metric],
                alpha=0.72 if metric == "TIL" else 0.56,
                edgecolor="white",
                linewidth=0.35,
                label=metric,
            )
            xs, ys = linear_fit_line(sub["rna_z"].to_numpy(dtype=float), sub["metric_z"].to_numpy(dtype=float))
            if xs.size:
                ax.plot(xs, ys, color=metric_colors[metric], linewidth=1.75, alpha=0.95)
        family = str(route_row.get("Family6", ""))
        class_color = IMMUNE_CLASS_COLORS.get(evidence, FAMILY_COLORS.get(family, "#64748B"))
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
        ax.spines["left"].set_color(class_color)
        ax.spines["bottom"].set_color(class_color)
        apply_axis_style(ax)
        ax.axhline(0, color="#CBD5E1", linewidth=0.7, linestyle="--")
        st = stats.loc[stats["gene"].eq(gene)]
        til_s = st.loc[st["metric"].eq("TIL")].iloc[0]
        tmb_s = st.loc[st["metric"].eq("TMB")].iloc[0]
        short_evidence = evidence.replace("Immune | ", "")
        ax.set_title(f"{gene} | {short_evidence}", loc="left", fontsize=9.4, fontweight="bold", color=TEXT, pad=6)
        ax.text(
            0.02,
            0.04,
            f"TIL rho={fmt_rho(til_s['rho'])}, q={fmt_p(til_s['q'])}\nTMB rho={fmt_rho(tmb_s['rho'])}, q={fmt_p(tmb_s['q'])}",
            transform=ax.transAxes,
            fontsize=6.9,
            color=MUTED,
            ha="left",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#E5E7EB", "alpha": 0.86},
        )

    for ax in axes[len(genes) :]:
        ax.axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc="lower center", ncol=2, frameon=False, fontsize=9.2)
    fig.suptitle(
        "Sankey-aligned immune TIL/TMB correlation evidence",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.02,
        0.958,
        f"Only immune routes retained in the final Sankey are shown. TMB correlations exclude {len(tmb_outliers)} recurrent zero/outlier samples; TIL uses all matched samples.",
        ha="left",
        fontsize=10,
        color=MUTED,
    )
    fig.text(0.5, 0.04, "RNA expression z-score", ha="center", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.text(0.024, 0.52, "Immune metric z-score", va="center", rotation="vertical", fontsize=10.5, fontweight="bold", color=TEXT)
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.10, top=0.90, wspace=0.22, hspace=0.34)
    save_to_all(fig, "Figure_S_Immune_TIL_TMB_sankey_aligned_biomarkers")
    plt.close(fig)
    outlier_df = pd.DataFrame({"Sample": tmb_outliers, "exclusion_scope": "TMB only", "reason": "tmb_log1p <= 0"})
    return stats, outlier_df


def reconciliation_table(route: pd.DataFrame) -> pd.DataFrame:
    sankey = (
        route[["route_type", "gene_clean"]]
        .drop_duplicates()
        .rename(columns={"route_type": "domain", "gene_clean": "gene"})
    )
    manifest_path = OUT_ROOTS[-1] / "tables" / "individual_biomarker_evidence_manifest.csv"
    if not manifest_path.exists():
        return sankey.assign(status="sankey_only")
    old = pd.read_csv(manifest_path)
    old = old[["domain", "gene"]].drop_duplicates()
    rows = []
    for domain in sorted(set(sankey["domain"]).union(old["domain"])):
        sgenes = set(sankey.loc[sankey["domain"].eq(domain), "gene"])
        ogenes = set(old.loc[old["domain"].eq(domain), "gene"])
        for gene in sorted(sgenes & ogenes):
            rows.append({"domain": domain, "gene": gene, "status": "matched_old_individual_and_sankey"})
        for gene in sorted(sgenes - ogenes):
            rows.append({"domain": domain, "gene": gene, "status": "sankey_only_in_old_outputs"})
        for gene in sorted(ogenes - sgenes):
            rows.append({"domain": domain, "gene": gene, "status": "old_individual_only_not_in_sankey"})
    return pd.DataFrame(rows)


def supplement_vs_sankey_table(route: pd.DataFrame, supplement_manifest: pd.DataFrame) -> pd.DataFrame:
    sankey = (
        route[["route_type", "gene_clean"]]
        .drop_duplicates()
        .rename(columns={"route_type": "domain", "gene_clean": "gene"})
    )
    supp = supplement_manifest[["domain", "gene"]].drop_duplicates()
    rows = []
    for domain in sorted(set(sankey["domain"]).union(supp["domain"])):
        sgenes = set(sankey.loc[sankey["domain"].eq(domain), "gene"])
        pgenes = set(supp.loc[supp["domain"].eq(domain), "gene"])
        for gene in sorted(sgenes & pgenes):
            rows.append({"domain": domain, "gene": gene, "status": "matched_supplement_and_sankey"})
        for gene in sorted(sgenes - pgenes):
            rows.append({"domain": domain, "gene": gene, "status": "sankey_only"})
        for gene in sorted(pgenes - sgenes):
            rows.append({"domain": domain, "gene": gene, "status": "supplement_only"})
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    route = pd.read_csv(SANKY_TABLE)
    os_manifest = plot_os_km(route)
    dep_manifest, dep_curve = plot_depmap(route)
    immune_stats, tmb_outliers = plot_immune(route)
    recon = reconciliation_table(route)

    supplement_manifest = pd.concat(
        [
            os_manifest[["domain", "gene", "route_evidence"]],
            dep_manifest[["domain", "gene", "route_evidence"]],
            immune_stats[["domain", "gene", "route_evidence"]].drop_duplicates(),
        ],
        ignore_index=True,
    ).drop_duplicates()
    supplement_recon = supplement_vs_sankey_table(route, supplement_manifest)
    write_table_all(os_manifest, "supplement_OS_KM_sankey_aligned_source_table.csv")
    write_table_all(dep_manifest, "supplement_DepMap_sankey_aligned_manifest.csv")
    write_table_all(dep_curve, "supplement_DepMap_sankey_aligned_dose_curve_table.csv")
    write_table_all(immune_stats, "supplement_Immune_TIL_TMB_sankey_aligned_correlation_table.csv")
    write_table_all(tmb_outliers, "supplement_Immune_TMB_outlier_samples_excluded.csv")
    write_table_all(recon, "individual_vs_sankey_evidence_reconciliation.csv")
    write_table_all(supplement_recon, "supplement_vs_sankey_evidence_reconciliation.csv")
    write_table_all(supplement_manifest, "sankey_aligned_supplement_manifest.csv")

    print("Saved Sankey-aligned one-page supplementary evidence figures:")
    for root in OUT_ROOTS:
        print(f"  {root / 'supplement_onepage'}")


if __name__ == "__main__":
    main()
