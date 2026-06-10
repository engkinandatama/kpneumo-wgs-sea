import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Headless backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist

# ============================================================
# Publication-Quality AMR Genotypic Heatmap
# Input : results/amr/summary_abricate.tab
#         config/samples.tsv
#         results/downstream/metadata/metadata_summary.tsv
# Output: results/downstream/figures/amr_heatmap.pdf
#         results/downstream/figures/amr_heatmap.png
# ============================================================

def load_abricate_summary(path: str) -> pd.DataFrame:
    """Parse ABRicate summary table into a clean presence/absence matrix."""
    df = pd.read_csv(path, sep="\t")
    # First column is file path, extract sample_id from filename
    df.index = df.iloc[:, 0].apply(
        lambda x: os.path.basename(x).replace("_abricate.tab", "")
    )
    df = df.drop(df.columns[0], axis=1)
    # Keep only gene columns (drop NUM_FOUND)
    df = df.drop(columns=["NUM_FOUND"], errors="ignore")
    # Convert presence/absence: "." = 0, numeric = 1
    df = df.apply(lambda col: col.map(lambda x: 0 if x == "." or pd.isna(x) else 1))
    return df


def filter_genes(matrix: pd.DataFrame, min_presence: int = 1) -> pd.DataFrame:
    """Keep only genes present in at least min_presence samples."""
    return matrix.loc[:, matrix.sum() >= min_presence]


def main():
    global snakemake
    if 'snakemake' not in globals():
        class MockSnakemake:
            input = type('Input', (), {
                'abricate': 'results/amr/summary_abricate.tab',
                'samples': 'config/samples.tsv',
                'metadata': 'results/downstream/metadata/metadata_summary.tsv'
            })()
            output = type('Output', (), {
                'pdf': 'results/downstream/figures/amr_heatmap.pdf',
                'png': 'results/downstream/figures/amr_heatmap.png'
            })()
        snakemake = MockSnakemake()

    abricate_path = snakemake.input.abricate
    samples_path  = snakemake.input.samples
    metadata_path = snakemake.input.metadata
    out_pdf       = snakemake.output.pdf
    out_png       = snakemake.output.png

    # Load data
    matrix = load_abricate_summary(abricate_path)
    # Simplify gene names by removing prefix
    matrix.columns = [c.replace("Klebsiella_pneumoniae_", "") for c in matrix.columns]
    matrix = filter_genes(matrix, min_presence=1)

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")
    meta_df    = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")

    # Muted, professional color palettes
    country_colors = {
        "Indonesia": "#2c7bb6",  # Elegant Blue
        "Malaysia":  "#1a9641",  # Elegant Green
        "Vietnam":   "#d7191c",  # Elegant Red
        "Thailand":  "#fdae61",  # Elegant Orange
    }
    species_colors = {
        "Klebsiella pneumoniae": "#252525",  # Dark Charcoal
        "Klebsiella quasipneumoniae subsp. quasipneumoniae": "#41b6c4",  # Muted Teal
    }
    species_labels = {
        "Klebsiella pneumoniae": "K. pneumoniae",
        "Klebsiella quasipneumoniae subsp. quasipneumoniae": "K. quasipneumoniae",
    }

    # Sort samples initially
    common_samples = [s for s in matrix.index if s in samples_df.index and s in meta_df.index]
    matrix = matrix.loc[common_samples]

    # Group genes by category for logical column ordering
    carb_genes = sorted([c for c in matrix.columns
                  if any(kw in c for kw in ["NDM", "KPC", "OXA", "GES", "IMP", "VIM", "BRP"])])
    esbl_genes = sorted([c for c in matrix.columns
                  if any(kw in c for kw in ["CTX-M", "TEM", "SHV", "VEB", "OKP"])])
    flq_genes  = sorted([c for c in matrix.columns
                  if any(kw in c for kw in ["Qnr", "qnr", "oqx"])])
    mcr_genes  = sorted([c for c in matrix.columns
                  if c.lower().startswith("mcr")])
    other_genes = sorted([c for c in matrix.columns
                   if c not in carb_genes + esbl_genes + flq_genes + mcr_genes])
    
    col_order = carb_genes + esbl_genes + flq_genes + mcr_genes + other_genes
    matrix = matrix[col_order]

    # Cluster rows using Jaccard distance
    Z = linkage(pdist(matrix.values, "jaccard"), method="average")
    order = leaves_list(Z)
    matrix = matrix.iloc[order]
    meta_df = meta_df.loc[matrix.index]

    n_samples, n_genes = matrix.shape

    # Set publication-quality font
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 8,
        "figure.dpi": 300,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "patch.linewidth": 0.8,
    })

    # Layout constants (inches) - expanded margins to prevent label overlap
    MARGIN_L = 0.15
    MARGIN_R = 3.0
    MARGIN_T = 0.60
    MARGIN_B = 2.0
    DEND_W   = 1.0
    ANNOT_W  = 0.22
    GAP      = 0.06
    HEAT_W   = n_genes * 0.165
    CAT_H    = 0.30
    HEAT_H   = n_samples * 0.30

    FIG_W = MARGIN_L + DEND_W + GAP + ANNOT_W + GAP + ANNOT_W + GAP + HEAT_W + MARGIN_R
    FIG_H = MARGIN_B + HEAT_H + CAT_H + MARGIN_T

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    # Define GridSpec
    gs = GridSpec(
        2, 4,
        left=MARGIN_L / FIG_W,
        right=1.0 - MARGIN_R / FIG_W,
        bottom=MARGIN_B / FIG_H,
        top=1.0 - MARGIN_T / FIG_H,
        width_ratios=[DEND_W, ANNOT_W, ANNOT_W, HEAT_W],
        height_ratios=[CAT_H, HEAT_H],
        wspace=GAP / np.mean([DEND_W, ANNOT_W, HEAT_W]),
        hspace=0.0
    )

    ax_dend   = fig.add_subplot(gs[1, 0])
    ax_sp     = fig.add_subplot(gs[1, 1])
    ax_co     = fig.add_subplot(gs[1, 2])
    ax_heat   = fig.add_subplot(gs[1, 3])
    ax_cat    = fig.add_subplot(gs[0, 3])

    # Row dendrogram
    dendrogram(
        Z, orientation="left", ax=ax_dend,
        color_threshold=0, above_threshold_color="#555555",
        link_color_func=lambda k: "#555555",
        no_labels=True
    )
    ax_dend.set_axis_off()

    # Annotation bar drawing helper
    def draw_annotation_bar(ax, values, palette, title):
        colors = [palette.get(v, "#cccccc") for v in values]
        for i, c in enumerate(colors):
            ax.barh(i + 0.5, 1.0, color=c, height=0.9, edgecolor="white", linewidth=0.5)
        ax.set_xlim(0, 1.0)
        ax.set_ylim(0, n_samples)
        ax.set_axis_off()
        ax.text(0.5, 1.02, title, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=7.5, fontweight="bold", rotation=90)

    draw_annotation_bar(ax_sp, [meta_df.loc[s, "species"] for s in matrix.index], species_colors, "Species")
    draw_annotation_bar(ax_co, [samples_df.loc[s, "country"] for s in matrix.index], country_colors, "Country")

    # Main Heatmap
    cmap = ListedColormap(["#f1f5f9", "#1e3a8a"])  # slate-white (absent) -> deep navy (present)
    im = ax_heat.imshow(
        matrix.values,
        aspect="auto",
        cmap=cmap,
        vmin=0, vmax=1,
        interpolation="nearest"
    )

    # Grid lines
    ax_heat.set_xticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_samples, 1), minor=True)
    ax_heat.grid(which="minor", color="#cbd5e1", linewidth=0.4)
    ax_heat.tick_params(which="minor", length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(True)
        spine.set_color("#475569")
        spine.set_linewidth(0.8)

    # Tick labels
    ax_heat.set_xticks(range(n_genes))
    ax_heat.set_xticklabels(col_order, rotation=90, fontsize=7.5, fontweight="bold")
    ax_heat.set_yticks(range(n_samples))
    ax_heat.set_yticklabels(matrix.index, fontsize=7.5)
    ax_heat.yaxis.set_ticks_position("right")
    ax_heat.tick_params(axis="both", which="both", length=2, color="#475569", pad=3)

    # Category split lines
    offsets = [len(carb_genes), len(carb_genes) + len(esbl_genes),
               len(carb_genes) + len(esbl_genes) + len(flq_genes),
               len(carb_genes) + len(esbl_genes) + len(flq_genes) + len(mcr_genes)]
    for x in offsets:
        if 0 < x < n_genes:
            ax_heat.axvline(x=x - 0.5, color="#1e293b", linewidth=1.2, linestyle="--", alpha=0.8)

    # Category labels above heatmap
    ax_cat.set_xlim(ax_heat.get_xlim())
    ax_cat.set_ylim(0, 1)
    ax_cat.set_axis_off()

    def add_category_header(start, end, text):
        if end > start:
            mid = (start + end - 1) / 2
            ax_cat.plot([start - 0.35, end - 0.65], [0.15, 0.15], color="#334155", linewidth=1.0)
            ax_cat.plot([start - 0.35, start - 0.35], [0.15, 0.35], color="#334155", linewidth=1.0)
            ax_cat.plot([end - 0.65, end - 0.65], [0.15, 0.35], color="#334155", linewidth=1.0)
            ax_cat.text(mid, 0.45, text, ha="center", va="bottom", fontsize=7.5,
                        fontweight="bold", color="#1e293b")

    add_category_header(0, len(carb_genes), "Carbapenems")
    add_category_header(len(carb_genes), len(carb_genes) + len(esbl_genes), "ESBLs")
    add_category_header(len(carb_genes) + len(esbl_genes), len(carb_genes) + len(esbl_genes) + len(flq_genes), "Quinolones")
    add_category_header(len(carb_genes) + len(esbl_genes) + len(flq_genes), len(carb_genes) + len(esbl_genes) + len(flq_genes) + len(mcr_genes), "MCR")
    add_category_header(len(carb_genes) + len(esbl_genes) + len(flq_genes) + len(mcr_genes), n_genes, "Other AMR")

    # Legend Panel - Shifted further right to prevent overlap with y-tick labels
    # Gap for y-tick labels is ~1.2 inches. So we place the legend at 1.4 inches from the heatmap.
    LEGEND_L = 1.4  # Inches from right edge of heatmap to left edge of legend axis
    ax_leg = fig.add_axes([
        1.0 - (MARGIN_R - LEGEND_L) / FIG_W,
        MARGIN_B / FIG_H,
        (MARGIN_R - LEGEND_L - 0.15) / FIG_W,
        HEAT_H / FIG_H
    ])
    ax_leg.set_axis_off()

    # Build legend handles
    legend_patches = []
    # Species
    legend_patches.append(mpatches.Patch(facecolor="white", edgecolor="white", label="Species"))
    for sp, c in species_colors.items():
        legend_patches.append(mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.5,
                                             label=species_labels.get(sp, sp)))
    legend_patches.append(mpatches.Patch(facecolor="white", edgecolor="white", label=""))
    # Country
    legend_patches.append(mpatches.Patch(facecolor="white", edgecolor="white", label="Country"))
    for co, c in country_colors.items():
        legend_patches.append(mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.5, label=co))
    legend_patches.append(mpatches.Patch(facecolor="white", edgecolor="white", label=""))
    # Presence
    legend_patches.append(mpatches.Patch(facecolor="white", edgecolor="white", label="AMR Presence"))
    legend_patches.append(mpatches.Patch(facecolor="#1e3a8a", edgecolor="black", linewidth=0.5, label="Present"))
    legend_patches.append(mpatches.Patch(facecolor="#f1f5f9", edgecolor="#cbd5e1", linewidth=0.5, label="Absent"))

    ax_leg.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        frameon=False,
        fontsize=7.5,
        borderpad=0.0,
        handlelength=1.1,
        labelspacing=0.35
    )

    # Title
    fig.suptitle(
        "Resistome Profile of Klebsiella pneumoniae Complex Clinical Isolates\n"
        "from Southeast Asia (n=20)",
        fontsize=11, fontweight="bold",
        x=(MARGIN_L + DEND_W + GAP + ANNOT_W + GAP + ANNOT_W + GAP + HEAT_W/2) / FIG_W,
        y=1.0 - 0.20 / FIG_H,
        ha="center"
    )

    # Save to outputs
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[amr_heatmap] Heatmap saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
