import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Headless backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import SymLogNorm
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
import os

# ============================================================
# Publication-Quality SNP Pairwise Distance Clustermap
# Input : results/phylogeny/core.tab  (snippy-core output)
#         config/samples.tsv
#         results/downstream/metadata/metadata_summary.tsv
# Output: results/downstream/snp_distance/snp_matrix.tsv
#         results/downstream/snp_distance/snp_heatmap.pdf
#         results/downstream/snp_distance/snp_heatmap.png
# ============================================================

def parse_core_tab(path: str) -> pd.DataFrame:
    """
    Compute pairwise SNP distance matrix from snippy-core core.tab file.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Rename columns to remove '_snippy' suffix
    df.columns = [c.replace("_snippy", "") for c in df.columns]
    
    # Drop positional metadata columns
    meta_cols = ["CHR", "POS", "REF"]
    sample_cols = [c for c in df.columns if c not in meta_cols]
    genotype_matrix = df[sample_cols]

    n = len(sample_cols)
    dist = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            diff = (genotype_matrix.iloc[:, i] != genotype_matrix.iloc[:, j]).sum()
            dist[i, j] = diff
            dist[j, i] = diff

    return pd.DataFrame(dist, index=sample_cols, columns=sample_cols)


def main():
    global snakemake
    if 'snakemake' not in globals():
        class MockSnakemake:
            input = type('Input', (), {
                'core_tab': 'results/phylogeny/core.tab',
                'samples': 'config/samples.tsv',
                'metadata': 'results/downstream/metadata/metadata_summary.tsv'
            })()
            output = type('Output', (), {
                'tsv': 'results/downstream/snp_distance/snp_matrix.tsv',
                'pdf': 'results/downstream/snp_distance/snp_heatmap.pdf',
                'png': 'results/downstream/snp_distance/snp_heatmap.png'
            })()
        snakemake = MockSnakemake()

    core_tab = snakemake.input.core_tab
    out_tsv  = snakemake.output.tsv
    out_pdf  = snakemake.output.pdf
    out_png  = snakemake.output.png

    samples_path  = snakemake.input.samples
    metadata_path = snakemake.input.metadata

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")
    meta_df    = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")

    print("[snp_distance] Parsing core.tab ...")
    dist_matrix = parse_core_tab(core_tab)
    dist_matrix.to_csv(out_tsv, sep="\t")
    print(f"[snp_distance] Matrix saved → {out_tsv}")

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

    # Align to common samples
    common  = [s for s in dist_matrix.index if s in meta_df.index and s in samples_df.index]
    dist_matrix = dist_matrix.loc[common, common].astype(float)
    meta_df = meta_df.loc[dist_matrix.index]

    n = len(dist_matrix)

    # Hierarchical clustering (average linkage on condensed distance matrix)
    condensed = squareform(dist_matrix.values, checks=False)
    Z         = linkage(condensed, method="average")
    order     = leaves_list(Z)
    dist_matrix = dist_matrix.iloc[order, order]
    meta_df = meta_df.loc[dist_matrix.index]

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

    # Layout constants (inches) - expanded margins to prevent label overlap
    MARGIN_L = 0.15
    MARGIN_R = 3.2
    MARGIN_T = 0.80
    MARGIN_B = 1.60
    DEND_W   = 0.90
    ANNOT_W  = 0.22
    GAP      = 0.05
    HEAT_SZ  = 5.4

    # Calculate figure dimensions
    FIG_W = MARGIN_L + DEND_W + GAP + ANNOT_W + GAP + ANNOT_W + GAP + HEAT_SZ + MARGIN_R
    COL_HDR_H = ANNOT_W + GAP + ANNOT_W + GAP + DEND_W
    FIG_H = MARGIN_B + HEAT_SZ + COL_HDR_H + MARGIN_T

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    # Define gridspec for symmetric matrix plotting
    gs = GridSpec(
        4, 4,
        left=MARGIN_L / FIG_W,
        right=1.0 - MARGIN_R / FIG_W,
        bottom=MARGIN_B / FIG_H,
        top=1.0 - MARGIN_T / FIG_H,
        width_ratios=[DEND_W, ANNOT_W, ANNOT_W, HEAT_SZ],
        height_ratios=[DEND_W, ANNOT_W, ANNOT_W, HEAT_SZ],
        wspace=GAP / np.mean([DEND_W, ANNOT_W, HEAT_SZ]),
        hspace=GAP / np.mean([DEND_W, ANNOT_W, HEAT_SZ])
    )

    ax_row_dend = fig.add_subplot(gs[3, 0])
    ax_row_sp   = fig.add_subplot(gs[3, 1])
    ax_row_co   = fig.add_subplot(gs[3, 2])
    
    ax_col_dend = fig.add_subplot(gs[0, 3])
    ax_col_sp   = fig.add_subplot(gs[1, 3])
    ax_col_co   = fig.add_subplot(gs[2, 3])
    
    ax_heat     = fig.add_subplot(gs[3, 3])

    # Row & Column Dendrograms
    dendrogram(
        Z, orientation="left", ax=ax_row_dend,
        color_threshold=0, above_threshold_color="#555555",
        link_color_func=lambda k: "#555555",
        no_labels=True
    )
    ax_row_dend.set_axis_off()

    dendrogram(
        Z, orientation="top", ax=ax_col_dend,
        color_threshold=0, above_threshold_color="#555555",
        link_color_func=lambda k: "#555555",
        no_labels=True
    )
    ax_col_dend.set_axis_off()

    # Row Annotation bars
    def draw_row_annot(ax, values, palette):
        colors = [palette.get(v, "#cccccc") for v in values]
        for i, c in enumerate(colors):
            ax.barh(i + 0.5, 1.0, color=c, height=0.9, edgecolor="white", linewidth=0.5)
        ax.set_xlim(0, 1.0)
        ax.set_ylim(0, n)
        ax.set_axis_off()

    draw_row_annot(ax_row_sp, [meta_df.loc[s, "species"] for s in dist_matrix.index], species_colors)
    draw_row_annot(ax_row_co, [samples_df.loc[s, "country"] for s in dist_matrix.index], country_colors)

    # Column Annotation bars
    def draw_col_annot(ax, values, palette):
        colors = [palette.get(v, "#cccccc") for v in values]
        for i, c in enumerate(colors):
            ax.bar(i + 0.5, 1.0, color=c, width=0.9, edgecolor="white", linewidth=0.5)
        ax.set_ylim(0, 1.0)
        ax.set_xlim(0, n)
        ax.set_axis_off()

    draw_col_annot(ax_col_sp, [meta_df.loc[s, "species"] for s in dist_matrix.columns], species_colors)
    draw_col_annot(ax_col_co, [samples_df.loc[s, "country"] for s in dist_matrix.columns], country_colors)

    # Main Heatmap
    max_val = dist_matrix.values.max()
    norm = SymLogNorm(linthresh=50, linscale=0.8, vmin=0, vmax=max_val, base=10)
    
    im = ax_heat.imshow(
        dist_matrix.values,
        aspect="auto",
        cmap="magma_r",  # Muted sequential colormap
        norm=norm,
        interpolation="nearest"
    )

    # Grid lines and border
    ax_heat.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax_heat.grid(which="minor", color="#cbd5e1", linewidth=0.4)
    ax_heat.tick_params(which="minor", length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(True)
        spine.set_color("#475569")
        spine.set_linewidth(0.8)

    # Plain black tick labels
    tick_labels = dist_matrix.index.tolist()
    ax_heat.set_xticks(range(n))
    ax_heat.set_xticklabels(tick_labels, rotation=90, fontsize=7.5, fontweight="semibold")
    ax_heat.set_yticks(range(n))
    ax_heat.set_yticklabels(tick_labels, fontsize=7.5, fontweight="semibold")
    ax_heat.yaxis.set_ticks_position("right")
    ax_heat.tick_params(axis="both", which="both", length=2, color="#475569", pad=3)

    # Lineage overlay rectangles
    # K. quasipneumoniae group
    qp_idx = [i for i, s in enumerate(dist_matrix.index) if "quasipneumoniae" in str(meta_df.loc[s, "species"])]
    if qp_idx:
        lo, hi = min(qp_idx) - 0.5, max(qp_idx) + 0.5
        ax_heat.add_patch(plt.Rectangle((lo, lo), hi - lo, hi - lo,
                                         fill=False, edgecolor="#14b8a6",
                                         linewidth=1.8, linestyle="--", zorder=5))
        # Label text box - placed neatly at the center of the cluster diagonal
        mid = lo + (hi - lo) / 2
        ax_heat.text(mid, mid, "K. quasipneumoniae",
                     ha="center", va="center", fontsize=7.5,
                     color="#004d40", fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#14b8a6", lw=0.8, alpha=0.85))

    # Vietnamese ST15 outbreak cohort
    vn_idx = [i for i, s in enumerate(dist_matrix.index) if samples_df.loc[s, "country"] == "Vietnam" if s in samples_df.index]
    if vn_idx:
        lo, hi = min(vn_idx) - 0.5, max(vn_idx) + 0.5
        ax_heat.add_patch(plt.Rectangle((lo, lo), hi - lo, hi - lo,
                                         fill=False, edgecolor="#ef4444",
                                         linewidth=1.8, linestyle=":", zorder=5))
        # Label text box - placed neatly at the center of the cluster diagonal
        mid = lo + (hi - lo) / 2
        ax_heat.text(mid, mid, "ST15 Clonal Cluster",
                     ha="center", va="center", fontsize=7.5,
                     color="#7f1d1d", fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#ef4444", lw=0.8, alpha=0.85))

    # Colorbar and Legend Panel - Shifted right to prevent overlap with y-tick labels
    # Gap for y-tick labels is ~1.2 inches. Colorbar and Legend start at 1.4 inches from heatmap.
    LEGEND_L = 1.4  # Inches from right edge of heatmap to left edge of colorbar/legend axes

    # Add Colorbar manually inside the legend axis area to keep it perfectly aligned
    cb_ax = fig.add_axes([
        1.0 - (MARGIN_R - LEGEND_L) / FIG_W,
        (MARGIN_B + HEAT_SZ - 1.8) / FIG_H,
        0.18 / FIG_W,
        1.5 / FIG_H
    ])
    cbar = fig.colorbar(im, cax=cb_ax, orientation="vertical")
    cbar.set_label("Pairwise SNP Distance (log scale)", fontsize=7.5, fontweight="bold", labelpad=10)
    cb_ax.yaxis.set_ticks_position("right")
    cb_ax.yaxis.set_label_position("right")  # Label on the right side of the colorbar, away from sample names
    cb_ax.tick_params(labelsize=7)

    # Legend Panel Axis
    ax_leg = fig.add_axes([
        1.0 - (MARGIN_R - LEGEND_L) / FIG_W,
        MARGIN_B / FIG_H,
        (MARGIN_R - LEGEND_L - 0.15) / FIG_W,
        (HEAT_SZ - 2.0) / FIG_H
    ])
    ax_leg.set_axis_off()

    # Custom legends
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

    ax_leg.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.52),
        frameon=False,
        fontsize=7.5,
        borderpad=0.0,
        handlelength=1.1,
        labelspacing=0.35
    )

    # Title
    fig.suptitle(
        "Clustered Core-genome Pairwise SNP Distance\n"
        "Klebsiella pneumoniae Complex — Southeast Asia (n=20)",
        fontsize=11, fontweight="bold",
        x=(MARGIN_L + DEND_W + GAP + ANNOT_W + GAP + ANNOT_W + GAP + HEAT_SZ/2) / FIG_W,
        y=1.0 - 0.25 / FIG_H,
        ha="center"
    )

    # Save to outputs
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[snp_distance] Heatmap saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
