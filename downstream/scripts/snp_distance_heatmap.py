import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import SymLogNorm
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
    sns.set_theme(style="white", rc={"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"]})

    # Muted, professional palette for metadata annotations
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

    # Map row and column annotations
    annot_df = pd.DataFrame(index=dist_matrix.index)
    annot_df["Country"] = dist_matrix.index.map(
        lambda s: country_colors.get(samples_df.loc[s, "country"] if s in samples_df.index else "Unknown", "#999999")
    )
    annot_df["Species"] = dist_matrix.index.map(
        lambda s: species_colors.get(meta_df.loc[s, "species"] if s in meta_df.index else "Unknown", "#7f8c8d")
    )

    # Use SymLogNorm to display both tight outbreak links (<10 SNPs)
    # and species-level differences (>180,000 SNPs) in the same color space.
    max_val = dist_matrix.values.max()
    norm = SymLogNorm(linthresh=10, linscale=1.0, vmin=0, vmax=max_val, base=10)

    print("[snp_distance] Plotting clustered heatmap with logarithmic normalization ...")
    
    # Generate hierarchical clustermap
    g = sns.clustermap(
        dist_matrix,
        row_colors=annot_df[["Species", "Country"]],
        col_colors=annot_df[["Species", "Country"]],
        cmap="rocket_r",       # Light colors for low distance, dark purple/black for high distance
        norm=norm,
        linewidths=0.2,
        linecolor="#f1f5f9",
        figsize=(11, 10),
        cbar_kws={"label": "Pairwise SNP distance (log scale)", "shrink": 0.5},
        tree_kws={"linewidths": 0.8},
    )

    # Rotate tick labels and adjust font sizes
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    # Adjust layout to make space for legends
    g.fig.subplots_adjust(top=0.92, right=0.82)
    g.fig.suptitle(
        "Clustered Core-genome Pairwise SNP Distance Heatmap\n"
        "Klebsiella pneumoniae Complex (n=20)",
        fontsize=12, fontweight="bold", y=0.98
    )

    # Add legends for row colors
    # 1. Country legend
    country_patches = [
        mpatches.Patch(color=c, label=k)
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
        mpatches.Patch(color=c, label=species_labels.get(k, k))
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

    # Save to outputs
    g.fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    g.fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[snp_distance] Heatmap saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
