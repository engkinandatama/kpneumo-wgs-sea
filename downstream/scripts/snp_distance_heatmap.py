import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# ============================================================
# SNP Pairwise Distance Heatmap
# Input : results/phylogeny/core.tab  (snippy-core output)
# Output: results/downstream/snp_distance/snp_matrix.tsv
#         results/downstream/snp_distance/snp_heatmap.pdf
#         results/downstream/snp_distance/snp_heatmap.png
# ============================================================

def parse_core_tab(path: str) -> pd.DataFrame:
    """
    Compute pairwise SNP distance matrix from snippy-core core.tab file.
    core.tab columns: CHR  POS  REF  Sample1_snippy  Sample2_snippy ...
    Each cell is the nucleotide call at that position.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Rename columns to remove '_snippy' suffix to match sample_id
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

    samples_path = snakemake.input.samples
    samples_df   = pd.read_csv(samples_path, sep="\t").set_index("sample_id")

    print("[snp_distance] Parsing core.tab ...")
    dist_matrix = parse_core_tab(core_tab)
    dist_matrix.to_csv(out_tsv, sep="\t")
    print(f"[snp_distance] Matrix saved → {out_tsv}")

    # Sort by country for visual grouping
    common = [s for s in dist_matrix.index if s in samples_df.index]
    order = sorted(common,
                   key=lambda s: (samples_df.loc[s, "country"], s))
    dist_matrix = dist_matrix.loc[order, order]

    country_colors = {
        "Indonesia": "#4285F4",
        "Malaysia":  "#34A853",
        "Vietnam":   "#FBBC05",
        "Thailand":  "#EA4335",
    }
    row_colors = dist_matrix.index.map(
        lambda s: country_colors.get(samples_df.loc[s, "country"]
                                     if s in samples_df.index else "Unknown", "#999999")
    )

    n = len(order)
    fig, ax = plt.subplots(figsize=(max(10, n * 0.55), max(8, n * 0.5)))

    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    heatmap = sns.heatmap(
        dist_matrix, ax=ax,
        cmap=cmap,
        annot=(n <= 25),          # show numbers if matrix is small enough
        fmt="d",
        linewidths=0.3,
        linecolor="#eeeeee",
        cbar_kws={"label": "Pairwise SNP distance", "shrink": 0.6},
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,  fontsize=7)
    ax.set_title(
        "Pairwise Core-genome SNP Distance Matrix\n"
        "K. pneumoniae Complex (n=20)",
        fontsize=11, fontweight="bold", pad=10
    )

    # Colour y-tick labels by country
    for tick_label, color in zip(ax.get_yticklabels(), row_colors):
        tick_label.set_color(color)

    # Country legend
    import matplotlib.patches as mpatches
    legend_patches = [
        mpatches.Patch(color=c, label=k)
        for k, c in country_colors.items()
    ]
    ax.legend(handles=legend_patches, title="Country",
              loc="upper left", bbox_to_anchor=(1.15, 1),
              frameon=True, fontsize=8, title_fontsize=9)

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[snp_distance] Heatmap saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
