import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import ListedColormap
import os

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
    matrix = filter_genes(matrix, min_presence=1)

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")
    meta_df    = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")

    # Muted, professional color palettes
    country_colors = {
        "Indonesia": "#3c8dbc",  # Steel blue
        "Malaysia":  "#00a65a",  # Emerald green
        "Vietnam":   "#f39c12",  # Amber
        "Thailand":  "#dd4b39",  # Crimson
    }
    species_colors = {
        "Klebsiella pneumoniae": "#2c3e50",  # Muted slate/navy
        "Klebsiella quasipneumoniae subsp. quasipneumoniae": "#16a085",  # Muted teal
    }

    # Sort samples initially, clustermap will reorder them hierarchically
    common_samples = [s for s in matrix.index if s in samples_df.index]
    matrix = matrix.loc[common_samples]

    # Group genes by category for logical column ordering
    carb_genes = [c for c in matrix.columns
                  if any(kw in c for kw in ["NDM", "KPC", "OXA", "GES", "IMP", "VIM", "MBL"])]
    esbl_genes = [c for c in matrix.columns
                  if any(kw in c for kw in ["CTX-M", "TEM", "SHV", "VEB"])]
    flq_genes  = [c for c in matrix.columns
                  if any(kw in c for kw in ["Qnr", "qnr", "aac(6')", "Gyra", "ParC"])]
    mcr_genes  = [c for c in matrix.columns
                  if c.lower().startswith("mcr")]
    other_genes = [c for c in matrix.columns
                   if c not in carb_genes + esbl_genes + flq_genes + mcr_genes]
    
    col_order = carb_genes + esbl_genes + flq_genes + mcr_genes + other_genes
    matrix = matrix[col_order]

    # Map row annotations (Species and Country)
    annot_df = pd.DataFrame(index=matrix.index)
    annot_df["Country"] = matrix.index.map(
        lambda s: country_colors.get(samples_df.loc[s, "country"] if s in samples_df.index else "Unknown", "#999999")
    )
    annot_df["Species"] = matrix.index.map(
        lambda s: species_colors.get(meta_df.loc[s, "species"] if s in meta_df.index else "Unknown", "#7f8c8d")
    )

    # Set publication-quality font
    sns.set_theme(style="white", rc={"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"]})

    # Plot dimensions
    fig_height = max(8, len(matrix) * 0.45)
    fig_width  = max(14, len(col_order) * 0.28)
    
    # Custom binary color map using ListedColormap: slate-white (absent) -> deep navy (present)
    cmap = ListedColormap(["#f8fafc", "#1e3a8a"])

    print("[amr_heatmap] Plotting clustered resistome heatmap ...")
    
    # Generate clustered heatmap (cluster rows/isolates, keep columns in functional order)
    g = sns.clustermap(
        matrix,
        row_cluster=True,
        col_cluster=False,
        row_colors=annot_df[["Species", "Country"]],
        cmap=cmap,
        linewidths=0.4,
        linecolor="#cbd5e1",
        figsize=(fig_width, fig_height),
        cbar=False,
        tree_kws={"linewidths": 0.8},
    )

    # Style axes
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8, fontweight="bold")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
    g.ax_heatmap.set_xlabel("Resistance Gene", fontsize=10, fontweight="bold", labelpad=10)
    g.ax_heatmap.set_ylabel("Isolate", fontsize=10, fontweight="bold", labelpad=10)
    
    # Add vertical separator lines between gene categories
    offsets = [len(carb_genes), len(carb_genes) + len(esbl_genes),
               len(carb_genes) + len(esbl_genes) + len(flq_genes),
               len(carb_genes) + len(esbl_genes) + len(flq_genes) + len(mcr_genes)]
    for x in offsets:
        if 0 < x < len(col_order):
            g.ax_heatmap.axvline(x=x, color="black", linewidth=1.5, linestyle="--", alpha=0.5)

    g.fig.suptitle(
        "Resistome Profile of Klebsiella pneumoniae Complex Clinical Isolates\n"
        "from Southeast Asia (n=20)",
        fontsize=13, fontweight="bold", y=0.98
    )

    # Adjust layout to prevent clipping and leave space for legends
    g.fig.subplots_adjust(top=0.90, right=0.84)

    # Add legends for row annotations
    # 1. Country legend
    country_patches = [
        mpatches.Patch(facecolor=c, edgecolor=c, label=k)
        for k, c in country_colors.items()
    ]
    g.fig.legend(
        handles=country_patches,
        title="Country of Origin",
        loc="upper right",
        bbox_to_anchor=(0.99, 0.85),
        frameon=True,
        fontsize=8,
        title_fontsize=9
    )

    # 2. Species legend (simplified labels)
    species_labels = {
        "Klebsiella pneumoniae": "K. pneumoniae",
        "Klebsiella quasipneumoniae subsp. quasipneumoniae": "K. quasipneumoniae"
    }
    species_patches = [
        mpatches.Patch(facecolor=c, edgecolor=c, label=species_labels.get(k, k))
        for k, c in species_colors.items()
    ]
    g.fig.legend(
        handles=species_patches,
        title="Species Identification",
        loc="upper right",
        bbox_to_anchor=(0.99, 0.70),
        frameon=True,
        fontsize=8,
        title_fontsize=9
    )

    # 3. Gene presence/absence legend
    amr_patches = [
        mpatches.Patch(facecolor="#1e3a8a", edgecolor="#1e3a8a", label="Present"),
        mpatches.Patch(facecolor="#f8fafc", edgecolor="#cbd5e1", label="Absent")
    ]
    g.fig.legend(
        handles=amr_patches,
        title="AMR Gene Presence",
        loc="upper right",
        bbox_to_anchor=(0.99, 0.55),
        frameon=True,
        fontsize=8,
        title_fontsize=9
    )

    # Save to outputs
    g.fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    g.fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[amr_heatmap] Heatmap saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
