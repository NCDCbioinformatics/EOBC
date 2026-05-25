from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


def env_path(name: str, default: Path | None = None) -> Path:
    value = os.environ.get(name)
    if value:
        return Path(value).expanduser()
    if default is not None:
        return Path(default).expanduser()
    raise RuntimeError(f"Set {name}; see config/paths_template.yml.")


BIOMARKER_ROOT = env_path("EOBC_BIOMARKER_ROOT")
FINAL_ANALYSIS_ROOT = env_path("EOBC_FINAL_ANALYSIS_DIR", BIOMARKER_ROOT / "final_analysis")
PHASE1_ROOT = env_path("EOBC_OS_PHASE1_ROOT", FINAL_ANALYSIS_ROOT / "os_phase1_subset_screen_v1")
TABLE_DIR = PHASE1_ROOT / "tables"
INT_SUMMARY = env_path(
    "EOBC_INTEGRATED_EVIDENCE_SUMMARY",
    FINAL_ANALYSIS_ROOT / "int_os_drug_imm_v1" / "tables" / "gene_integrated_evidence_summary.csv",
)
OUTPUT_ROOT = env_path("EOBC_OS_CONSENSUS_OUT", FINAL_ANALYSIS_ROOT / "os_consensus_analysis_v3")
PLOTS_DIR = OUTPUT_ROOT / "plots"
TABLES_DIR = OUTPUT_ROOT / "tables"

TOP_N = 30
MODALITIES = ["RNA", "METH"]
ENDPOINTS = ["overall_os", "os_5y", "os_10y"]
SCREEN_ORDER = [
    ("RNA", "overall_os"),
    ("RNA", "os_5y"),
    ("RNA", "os_10y"),
    ("METH", "overall_os"),
    ("METH", "os_5y"),
    ("METH", "os_10y"),
]
SCREEN_LABELS = {
    ("RNA", "overall_os"): "RNA\nOverall",
    ("RNA", "os_5y"): "RNA\n5y",
    ("RNA", "os_10y"): "RNA\n10y",
    ("METH", "overall_os"): "METH\nOverall",
    ("METH", "os_5y"): "METH\n5y",
    ("METH", "os_10y"): "METH\n10y",
}
ROW_LABELS = {
    ("RNA", "overall_os"): "RNA | Overall OS",
    ("RNA", "os_5y"): "RNA | 5-year OS",
    ("RNA", "os_10y"): "RNA | 10-year OS",
    ("METH", "overall_os"): "METH | Overall OS",
    ("METH", "os_5y"): "METH | 5-year OS",
    ("METH", "os_10y"): "METH | 10-year OS",
}

FAMILY_COLORS = {
    "Immune": "#2a9d8f",
    "Repair": "#d62828",
    "Glycolysis / TCA": "#f77f00",
    "Hormone signaling": "#577590",
    "Fatty acid": "#f4a261",
    "Kinase signaling": "#7b2cbf",
}
LAYER_EDGE = {
    "Transcriptome": "#0b6e4f",
    "Methylation": "#8d5524",
}
LAYER_MARKERS = {"Transcriptome": "o", "Methylation": "s"}
CMAP = LinearSegmentedColormap.from_list(
    "os_consensus", ["#69b3ff", "#f7f7f7", "#f4a261"], N=256
)


def ensure_dirs() -> None:
    for path in [OUTPUT_ROOT, PLOTS_DIR, TABLES_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def save_figure_bundle(fig: plt.Figure, output_path: Path, dpi: int = 320) -> None:
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def family_chip(label: str, family: str) -> str:
    short = {
        "Immune": "Imm",
        "Repair": "Rep",
        "Glycolysis / TCA": "Gly/TCA",
        "Hormone signaling": "Horm",
        "Fatty acid": "FA",
        "Kinase signaling": "Kin",
    }.get(family, family)
    return f"{label}  [{short}]"


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(INT_SUMMARY)
    meta = meta[["gene", "gene_label", "Layer", "Family6", "target_label"]].copy()
    meta["gene_label"] = meta["gene_label"].fillna(meta["gene"])
    return meta.drop_duplicates(subset=["gene"]).reset_index(drop=True)


def load_screen_tables() -> tuple[Dict[Tuple[str, str], pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame]]:
    volcano = {}
    top_hits = {}
    km_tables = {}
    for modality, endpoint in SCREEN_ORDER:
        key = (modality, endpoint)
        volcano[key] = pd.read_csv(TABLE_DIR / f"volcano_allsets_{modality}_{endpoint}.csv")
        volcano[key]["gene_list"] = volcano[key]["genes"].str.split("|", regex=False)
        top_hits[key] = pd.read_csv(TABLE_DIR / f"top_subsets_{modality}_{endpoint}.csv")
        km_tables[key] = pd.read_csv(TABLE_DIR / f"km_table_{modality}_{endpoint}.csv")
    return volcano, top_hits, km_tables


def gene_in_subset(df: pd.DataFrame, gene: str) -> pd.Series:
    return df["gene_list"].apply(lambda xs: gene in xs)


def build_screen_gene_metrics(meta: pd.DataFrame, volcano: Dict[Tuple[str, str], pd.DataFrame], top_hits: Dict[Tuple[str, str], pd.DataFrame], km_tables: Dict[Tuple[str, str], pd.DataFrame]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for modality, endpoint in SCREEN_ORDER:
        key = (modality, endpoint)
        df = volcano[key].sort_values("p_value_raw").reset_index(drop=True)
        topn = df.head(TOP_N).copy()
        top1_label = str(top_hits[key].iloc[0]["subset_label"])
        top1_genes = set(str(top_hits[key].iloc[0]["genes"]).split("|"))
        top1_hr_adj = float(top_hits[key].iloc[0]["hr_adj"])
        top1_p_adj = float(top_hits[key].iloc[0]["p_value_adj"])
        km = km_tables[key].copy()
        n_high = int((km["group"] == "High score").sum())
        n_low = int((km["group"] == "Low score").sum())

        for _, meta_row in meta.iterrows():
            gene = str(meta_row["gene"])
            gene_all = df[gene_in_subset(df, gene)].copy()
            gene_top = topn[gene_in_subset(topn, gene)].copy()

            if gene_all.empty:
                continue

            best = gene_all.sort_values("p_value_raw").iloc[0]
            if gene_top.empty:
                weighted_beta = 0.0
                mean_beta = 0.0
                support = 0.0
                topn_neglogp_mean = 0.0
            else:
                weights = np.maximum(gene_top["neglog10_p_raw"].to_numpy(dtype=float), 1e-6)
                weighted_beta = float(np.average(gene_top["beta_raw"], weights=weights))
                mean_beta = float(gene_top["beta_raw"].mean())
                support = float(gene_top.shape[0] / TOP_N)
                topn_neglogp_mean = float(gene_top["neglog10_p_raw"].mean())

            screen_contrib = weighted_beta * support
            rows.append(
                {
                    "gene": gene,
                    "gene_label": meta_row["gene_label"],
                    "Layer": meta_row["Layer"],
                    "Family6": meta_row["Family6"],
                    "target_label": meta_row["target_label"],
                    "modality": modality,
                    "endpoint": endpoint,
                    "screen": f"{modality}_{endpoint}",
                    "support_top30": support,
                    "weighted_beta_top30": weighted_beta,
                    "mean_beta_top30": mean_beta,
                    "top30_mean_neglog10p": topn_neglogp_mean,
                    "screen_contrib": screen_contrib,
                    "best_subset": best["subset_label"],
                    "best_beta_raw": float(best["beta_raw"]),
                    "best_neglog10p_raw": float(best["neglog10_p_raw"]),
                    "best_p_raw": float(best["p_value_raw"]),
                    "top1_hit": gene in top1_genes,
                    "top1_signature": top1_label,
                    "top1_hr_adj": top1_hr_adj,
                    "top1_p_adj": top1_p_adj,
                    "top1_n_high": n_high,
                    "top1_n_low": n_low,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(TABLES_DIR / "os_gene_screen_metrics.csv", index=False)
    return out


def classify_gene(row: pd.Series) -> str:
    strength = float(row["consensus_strength"])
    support = float(row["support_sum"])
    early = float(row["early_score"])
    late = float(row["late_score"])
    overall = float(row["overall_score"])
    consensus = float(row["consensus_score"])

    if strength < 0.12 or support < 0.20:
        return "Weak / context-dependent"

    if consensus >= 0:
        if early > 0 and late > 0:
            if abs(early - late) < 0.08:
                return "Persistent adverse"
            if early > late:
                return "Early-dominant adverse"
            return "Late-dominant adverse"
        if early > 0 and late <= 0:
            return "Early-specific adverse"
        if late > 0 and early <= 0:
            return "Late-specific adverse"
        if overall > 0:
            return "Overall-adverse / mixed"
        return "Adverse-leaning mixed"

    if early < 0 and late < 0:
        if abs(early - late) < 0.08:
            return "Persistent protective"
        if abs(early) > abs(late):
            return "Early-dominant protective"
        return "Late-dominant protective"
    if early < 0 and late >= 0:
        return "Early-specific protective"
    if late < 0 and early >= 0:
        return "Late-specific protective"
    if overall < 0:
        return "Overall-protective / mixed"
    return "Protective-leaning mixed"


def build_gene_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    pivot_contrib = metrics.pivot(index="gene", columns="screen", values="screen_contrib").fillna(0.0)
    pivot_support = metrics.pivot(index="gene", columns="screen", values="support_top30").fillna(0.0)

    meta = metrics[["gene", "gene_label", "Layer", "Family6", "target_label"]].drop_duplicates("gene").set_index("gene")
    summary = meta.copy()
    summary["consensus_score"] = pivot_contrib.sum(axis=1)
    summary["consensus_strength"] = pivot_contrib.abs().sum(axis=1)
    summary["support_sum"] = pivot_support.sum(axis=1)
    summary["support_mean"] = pivot_support.mean(axis=1)
    summary["early_score"] = pivot_contrib[["RNA_os_5y", "METH_os_5y"]].sum(axis=1)
    summary["late_score"] = pivot_contrib[["RNA_os_10y", "METH_os_10y"]].sum(axis=1)
    summary["overall_score"] = pivot_contrib[["RNA_overall_os", "METH_overall_os"]].sum(axis=1)
    summary["rna_score"] = pivot_contrib[["RNA_overall_os", "RNA_os_5y", "RNA_os_10y"]].sum(axis=1)
    summary["meth_score"] = pivot_contrib[["METH_overall_os", "METH_os_5y", "METH_os_10y"]].sum(axis=1)
    summary["top1_hits"] = metrics.groupby("gene")["top1_hit"].sum()
    summary["best_screen"] = metrics.sort_values("best_p_raw").groupby("gene").first()["screen"]
    summary["best_subset"] = metrics.sort_values("best_p_raw").groupby("gene").first()["best_subset"]
    summary["best_neglog10p"] = metrics.groupby("gene")["best_neglog10p_raw"].max()
    summary["consensus_class"] = summary.apply(classify_gene, axis=1)
    summary = summary.sort_values(
        ["consensus_strength", "support_sum", "best_neglog10p"],
        ascending=[False, False, False],
    ).reset_index()
    summary.to_csv(TABLES_DIR / "os_gene_consensus_summary.csv", index=False)
    return summary


def build_top_signature_matrix(summary: pd.DataFrame, top_hits: Dict[Tuple[str, str], pd.DataFrame], metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    signatures = []
    for modality, endpoint in SCREEN_ORDER:
        row = top_hits[(modality, endpoint)].iloc[0]
        genes = str(row["genes"]).split("|")
        signatures.append((modality, endpoint, str(row["subset_label"]), genes, float(row["hr_adj"]), float(row["p_value_adj"])))
    union_genes = []
    for _, _, _, genes, _, _ in signatures:
        for gene in genes:
            if gene not in union_genes:
                union_genes.append(gene)
    gene_order = summary.set_index("gene").loc[[g for g in union_genes if g in set(summary["gene"])]].sort_values(
        "consensus_strength", ascending=False
    ).index.tolist()

    for modality, endpoint, subset_label, genes, hr_adj, p_adj in signatures:
        km_row = metrics[
            (metrics["modality"] == modality)
            & (metrics["endpoint"] == endpoint)
        ].iloc[0]
        row = {
            "screen": f"{modality}_{endpoint}",
            "row_label": ROW_LABELS[(modality, endpoint)],
            "subset_label": subset_label,
            "hr_adj": hr_adj,
            "p_adj": p_adj,
            "n_high": int(km_row["top1_n_high"]),
            "n_low": int(km_row["top1_n_low"]),
        }
        for gene in gene_order:
            row[gene] = int(gene in genes)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(TABLES_DIR / "os_top_signature_chip_matrix.csv", index=False)
    return out


def plot_bubble_matrix(ax: plt.Axes, metrics: pd.DataFrame, summary: pd.DataFrame) -> None:
    gene_order = summary["gene"].tolist()
    metric_map = metrics.set_index(["gene", "screen"])
    screen_names = [f"{m}_{e}" for m, e in SCREEN_ORDER]
    max_abs = max(0.25, float(np.abs(metrics["screen_contrib"]).max()))
    norm = Normalize(vmin=-max_abs, vmax=max_abs)

    ax.set_facecolor("#fbfcfe")
    ax.axvspan(-0.5, 2.5, color="#edf6ff", alpha=0.75, zorder=0)
    ax.axvspan(2.5, 5.5, color="#fff4ea", alpha=0.75, zorder=0)

    for xi in range(len(screen_names) + 1):
        ax.axvline(xi - 0.5, color="#e5e7eb", linewidth=0.7, zorder=0)
    for yi in range(len(gene_order) + 1):
        ax.axhline(yi - 0.5, color="#eef2f7", linewidth=0.55, zorder=0)

    for y, gene in enumerate(gene_order):
        layer = summary.loc[summary["gene"] == gene, "Layer"].iloc[0]
        family = summary.loc[summary["gene"] == gene, "Family6"].iloc[0]
        marker = LAYER_MARKERS.get(layer, "o")
        for x, screen in enumerate(screen_names):
            row = metric_map.loc[(gene, screen)]
            support = float(row["support_top30"])
            contrib = float(row["screen_contrib"])
            top1_hit = bool(row["top1_hit"])

            if support <= 0:
                ax.scatter(
                    x,
                    y,
                    s=24,
                    facecolors="none",
                    edgecolors="#d1d5db",
                    linewidths=0.5,
                    marker=marker,
                    zorder=2,
                )
                continue

            ax.scatter(
                x,
                y,
                s=80 + support * 820,
                c=[CMAP(norm(contrib))],
                edgecolors=LAYER_EDGE.get(layer, "#4b5563"),
                linewidths=1.0,
                marker=marker,
                zorder=3,
            )
            if top1_hit:
                ax.scatter(
                    x,
                    y,
                    s=46,
                    c="#ffd166",
                    edgecolors="white",
                    linewidths=0.5,
                    marker="*",
                    zorder=4,
                )

        label = summary.loc[summary["gene"] == gene, "gene_label"].iloc[0]
        ax.text(
            -0.78,
            y,
            family_chip(label, family),
            ha="right",
            va="center",
            fontsize=8.8,
            color="#111827",
        )

    ax.text(1.0, -1.12, "RNA screens", ha="center", va="center", fontsize=10, weight="bold", color="#2563eb")
    ax.text(4.0, -1.12, "Methylation screens", ha="center", va="center", fontsize=10, weight="bold", color="#ea580c")

    ax.set_xlim(-1.25, len(screen_names) - 0.4)
    ax.set_ylim(len(gene_order) - 0.5, -1.4)
    ax.set_xticks(range(len(screen_names)))
    ax.set_xticklabels([SCREEN_LABELS[k] for k in SCREEN_ORDER], fontsize=9, weight="bold")
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "A. Gene-level consensus extracted from exhaustive subset screening",
        loc="left",
        fontsize=13.5,
        weight="bold",
        pad=14,
    )

    size_vals = [0.10, 0.30, 0.60]
    size_handles = [
        plt.scatter([], [], s=80 + s * 820, color="#9ca3af", edgecolors="#4b5563", linewidths=0.8, label=f"{int(s*100)}% of top {TOP_N}")
        for s in size_vals
    ]
    shape_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Transcriptome"], markeredgewidth=1, label="Transcriptome", markersize=8),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Methylation"], markeredgewidth=1, label="Methylation", markersize=8),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="#ffd166", markeredgecolor="white", markeredgewidth=0.6, label="Present in screen-leading signature", markersize=10),
    ]
    leg1 = ax.legend(handles=size_handles, title="Top-subset support", loc="upper left", bbox_to_anchor=(0.00, -0.08), ncol=3, frameon=False, fontsize=8, title_fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=shape_handles, title="Gene layer", loc="upper left", bbox_to_anchor=(0.00, -0.16), ncol=3, frameon=False, fontsize=8, title_fontsize=9)

    cax = ax.inset_axes([0.74, -0.20, 0.22, 0.03])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_label("Signed OS consensus within top subsets  (protective  ←  →  adverse)", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)


def plot_temporal_scatter(ax: plt.Axes, summary: pd.DataFrame) -> None:
    ax.set_facecolor("#ffffff")
    ax.axhspan(0, max(0.05, summary["late_score"].max() * 1.15), xmin=0.5, xmax=1.0, color="#fff0e6", alpha=0.55, zorder=0)
    ax.axhspan(summary["late_score"].min() * 1.15, 0, xmin=0.0, xmax=0.5, color="#ecf6ff", alpha=0.55, zorder=0)
    ax.axvspan(0, max(0.05, summary["early_score"].max() * 1.15), ymin=0.0, ymax=1.0, color="#fff7ed", alpha=0.28, zorder=0)
    ax.axvspan(summary["early_score"].min() * 1.15, 0, ymin=0.0, ymax=1.0, color="#eff6ff", alpha=0.28, zorder=0)
    ax.axhline(0, color="#6b7280", linestyle="--", linewidth=1.0)
    ax.axvline(0, color="#6b7280", linestyle="--", linewidth=1.0)

    size_scale = 1200 / max(summary["consensus_strength"].max(), 0.25)
    top_labels = summary.head(10).copy()

    for _, row in summary.iterrows():
        family = row["Family6"]
        layer = row["Layer"]
        ax.scatter(
            row["early_score"],
            row["late_score"],
            s=90 + row["consensus_strength"] * size_scale,
            c=FAMILY_COLORS.get(family, "#9ca3af"),
            marker=LAYER_MARKERS.get(layer, "o"),
            edgecolors=LAYER_EDGE.get(layer, "#374151"),
            linewidths=1.0,
            alpha=0.92,
            zorder=2,
        )

    placed: List[Tuple[float, float]] = []
    for _, row in top_labels.iterrows():
        x = float(row["early_score"])
        y = float(row["late_score"])
        dx = 0.015 if x >= 0 else -0.015
        dy = 0.012 if y >= 0 else -0.012
        lx = x + dx
        ly = y + dy
        for px, py in placed:
            if abs(lx - px) < 0.12 and abs(ly - py) < 0.08:
                ly += 0.05 if dy >= 0 else -0.05
        ax.text(
            lx,
            ly,
            row["gene_label"],
            fontsize=8.5,
            color="#111827",
            ha="left" if dx > 0 else "right",
            va="bottom" if dy > 0 else "top",
            zorder=3,
        )
        placed.append((lx, ly))

    ax.text(0.02, 0.97, "10y adverse", transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color="#b45309")
    ax.text(0.98, 0.03, "5y adverse", transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5, color="#c2410c")
    ax.text(0.02, 0.03, "5y/10y protective", transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color="#2563eb")

    ax.set_xlabel("Short-term OS consensus (5-year screens)", fontsize=10.5, weight="bold")
    ax.set_ylabel("Long-term OS consensus (10-year screens)", fontsize=10.5, weight="bold")
    ax.grid(True, color="#e5e7eb", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_title(
        "B. Temporal OS behavior of each biomarker",
        loc="left",
        fontsize=13.5,
        weight="bold",
        pad=14,
    )

    fam_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none", label=family, markersize=8)
        for family, color in FAMILY_COLORS.items()
    ]
    layer_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Transcriptome"], markeredgewidth=1, label="Transcriptome", markersize=8),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Methylation"], markeredgewidth=1, label="Methylation", markersize=8),
    ]
    leg1 = ax.legend(handles=fam_handles, title="Biological class", loc="upper left", bbox_to_anchor=(0.00, -0.20), ncol=2, frameon=False, fontsize=8, title_fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=layer_handles, title="Layer", loc="upper left", bbox_to_anchor=(0.60, -0.20), ncol=2, frameon=False, fontsize=8, title_fontsize=9)


def plot_signature_matrix(ax: plt.Axes, signature_df: pd.DataFrame, summary: pd.DataFrame) -> None:
    gene_cols = [c for c in signature_df.columns if c not in {"screen", "row_label", "subset_label", "hr_adj", "p_adj", "n_high", "n_low"}]
    summary_idx = summary.set_index("gene")
    ax.set_facecolor("#fbfcfe")

    for yi in range(signature_df.shape[0] + 1):
        ax.axhline(yi - 0.5, color="#e5e7eb", linewidth=0.7, zorder=0)
    for xi in range(len(gene_cols) + 1):
        ax.axvline(xi - 0.5, color="#eef2f7", linewidth=0.7, zorder=0)

    for y, (_, row) in enumerate(signature_df.iterrows()):
        ax.text(-0.85, y, row["row_label"], ha="right", va="center", fontsize=9.0, color="#111827", weight="bold")
        for x, gene in enumerate(gene_cols):
            val = int(row[gene])
            family = summary_idx.loc[gene, "Family6"]
            layer = summary_idx.loc[gene, "Layer"]
            edge = LAYER_EDGE.get(layer, "#4b5563")
            face = FAMILY_COLORS.get(family, "#d1d5db") if val == 1 else "#ffffff"
            chip = FancyBboxPatch(
                (x - 0.33, y - 0.22),
                0.66,
                0.44,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                facecolor=face,
                edgecolor=edge if val == 1 else "#d1d5db",
                linewidth=1.2 if val == 1 else 0.7,
                alpha=0.98,
            )
            ax.add_patch(chip)
        ax.text(
            len(gene_cols) + 0.55,
            y,
            f"Adj.HR {row['hr_adj']:.2f}  |  adj.P {format_p(row['p_adj'])}",
            ha="left",
            va="center",
            fontsize=8.6,
            color="#374151",
        )

    for x, gene in enumerate(gene_cols):
        ax.text(
            x,
            -0.82,
            summary_idx.loc[gene, "gene_label"],
            ha="right",
            va="center",
            rotation=45,
            fontsize=8.5,
            color="#111827",
            weight="bold",
        )

    ax.set_xlim(-1.35, len(gene_cols) + 2.35)
    ax.set_ylim(signature_df.shape[0] - 0.5, -1.15)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "C. Reusable gene architecture of the leading OS signatures",
        loc="left",
        fontsize=13.5,
        weight="bold",
        pad=14,
    )


def make_main_figure(metrics: pd.DataFrame, summary: pd.DataFrame, signature_df: pd.DataFrame) -> Path:
    fig = plt.figure(figsize=(18.4, 12.4))
    gs = GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[1.45, 1.0],
        height_ratios=[1.0, 0.72],
        wspace=0.24,
        hspace=0.34,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    plot_bubble_matrix(ax_a, metrics, summary)
    plot_temporal_scatter(ax_b, summary)
    plot_signature_matrix(ax_c, signature_df, summary)

    fig.suptitle(
        "EOBC OS biomarker consensus: integrating exhaustive subset screening into gene-level survival programs",
        fontsize=19,
        weight="bold",
        y=0.992,
    )
    fig.text(
        0.5,
        0.965,
        "New summary framework: each biomarker is scored by how repeatedly it appears in the strongest gene sets, which direction it drives raw OS association, and whether its effect is early, late, or persistent across RNA and methylation screens.",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#4b5563",
    )
    output = PLOTS_DIR / "Figure_01_OS_gene_consensus_framework.png"
    save_figure_bundle(fig, output)
    plt.close(fig)
    return output


def save_individual_panels(metrics: pd.DataFrame, summary: pd.DataFrame, signature_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10.6, 11.2))
    plot_bubble_matrix(ax, metrics, summary)
    fig.tight_layout()
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01A_OS_gene_consensus_matrix.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    plot_temporal_scatter(ax, summary)
    fig.tight_layout()
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01B_OS_temporal_scatter.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    plot_signature_matrix(ax, signature_df, summary)
    fig.tight_layout()
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01C_OS_signature_architecture.png")
    plt.close(fig)


def make_supplementary_heatmap(metrics: pd.DataFrame, summary: pd.DataFrame) -> Path:
    matrix = metrics.pivot(index="gene", columns="screen", values="screen_contrib")
    matrix = matrix[[f"{m}_{e}" for m, e in SCREEN_ORDER]]
    matrix = matrix.loc[summary["gene"]]
    display = matrix.rename(columns={f"{m}_{e}": SCREEN_LABELS[(m, e)].replace("\n", " | ") for m, e in SCREEN_ORDER})
    max_abs = max(0.25, float(np.abs(display.values).max()))

    fig, ax = plt.subplots(figsize=(10.6, 10.8))
    im = ax.imshow(display.values, aspect="auto", cmap=CMAP, vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(np.arange(display.shape[1]))
    ax.set_xticklabels(display.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(display.shape[0]))
    ax.set_yticklabels(summary["gene_label"], fontsize=8.5)
    ax.set_title("OS gene consensus heatmap (signed support-weighted effect)", fontsize=14, weight="bold", pad=12)
    for i in range(display.shape[0]):
        ax.text(display.shape[1] + 0.12, i, summary.iloc[i]["consensus_class"], va="center", fontsize=7.6, color="#374151")
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("Protective  ←  signed screen contribution  →  Adverse")
    fig.tight_layout(rect=[0, 0, 0.92, 1])
    output = PLOTS_DIR / "Figure_S01_OS_consensus_heatmap.png"
    save_figure_bundle(fig, output)
    plt.close(fig)
    return output


def write_summary(summary: pd.DataFrame, figure_path: Path) -> None:
    top = summary.head(12)
    lines = [
        "# OS consensus analysis v3",
        "",
        "This analysis converts exhaustive subset screening into gene-level OS consensus metrics.",
        "",
        "## Strongest consensus genes",
        "",
    ]
    for _, row in top.iterrows():
        lines.append(
            f"- {row['gene_label']} ({row['Family6']}, {row['Layer']}): "
            f"{row['consensus_class']} | score={row['consensus_score']:.3f}, "
            f"strength={row['consensus_strength']:.3f}, support_sum={row['support_sum']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Main figure",
            "",
            f"- {figure_path}",
        ]
    )
    (OUTPUT_ROOT / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    meta = load_meta()
    volcano, top_hits, km_tables = load_screen_tables()
    metrics = build_screen_gene_metrics(meta, volcano, top_hits, km_tables)
    summary = build_gene_summary(metrics)
    signature_df = build_top_signature_matrix(summary, top_hits, metrics)
    figure_path = make_main_figure(metrics, summary, signature_df)
    save_individual_panels(metrics, summary, signature_df)
    make_supplementary_heatmap(metrics, summary)
    write_summary(summary, figure_path)
    print("OS consensus analysis v3 complete.")
    print(f"Saved to {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
