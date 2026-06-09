import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import sys
import os

# ============================================================
# AMR Heatmap - Distribution of Resistance Genes across Isolates
# Input : results/amr/summary_abricate.tab  +  config/samples.tsv
# Output: results/downstream/figures/amr_heatmap.pdf  +  .png
# ============================================================

def load_abricate_summary(path: str) -> pd.DataFrame:
    """Parse ABRicate summary table into a clean presence/absence matrix."""
    df = pd.read_csv(path, sep="\t", comment="#", header=0)
    # First column is file path, extract sample_id from filename
    df.index = df.iloc[:, 0].apply(
        lambda x: os.path.basename(x).replace("_abricate.tab", "")
    )
    df = df.drop(df.columns[0], axis=1)
    # Keep only gene columns (drop NUM_FOUND)
    df = df.drop(columns=["NUM_FOUND"], errors="ignore")
    # Convert presence/absence: "." = 0, numeric = 1
    df = df.applymap(lambda x: 0 if x == "." or pd.isna(x) else 1)
    return df


def filter_genes(matrix: pd.DataFrame, min_presence: int = 1) -> pd.DataFrame:
    """Keep only genes present in at least min_presence samples."""
    return matrix.loc[:, matrix.sum() >= min_presence]


def add_country_annotation(ax, samples_df: pd.DataFrame, order: list,
                            country_colors: dict, fig_width: float) -> None:
    """Add a country color strip to the left side of the heatmap."""
    country_col = samples_df.loc[order, "country"]
    colors = [country_colors[c] for c in country_col]
    for i, color in enumerate(colors):
        ax.add_patch(
            mpatches.Rectangle(
                xy=(-1.5, i), width=1, height=1,
                color=color, transform=ax.transData, clip_on=False
            )
        )


def main():
    # --- Paths ---
    abricate_path = snakemake.input.abricate
    samples_path  = snakemake.input.samples
    out_pdf       = snakemake.output.pdf
    out_png       = snakemake.output.png

    # --- Load data ---
    matrix = load_abricate_summary(abricate_path)
    matrix = filter_genes(matrix, min_presence=1)

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")

    # --- Country colour palette ---
    country_colors = {
        "Indonesia": "#4285F4",
        "Malaysia":  "#34A853",
        "Vietnam":   "#FBBC05",
        "Thailand":  "#EA4335",
    }

    # --- Sort samples by country then sample_id ---
    common_samples = [s for s in matrix.index if s in samples_df.index]
    order = sorted(common_samples,
                   key=lambda s: (samples_df.loc[s, "country"], s))
    matrix = matrix.loc[order]

    # --- Group important gene categories for column ordering ---
    carb_genes = [c for c in matrix.columns
                  if any(kw in c for kw in
                         ["NDM", "KPC", "OXA", "GES", "IMP", "VIM", "MBL"])]
    esbl_genes = [c for c in matrix.columns
                  if any(kw in c for kw in ["CTX-M", "TEM", "SHV", "VEB"])]
    flq_genes  = [c for c in matrix.columns
                  if any(kw in c for kw in ["Qnr", "qnr", "aac(6')", "Gyra", "ParC"])]
    mcr_genes  = [c for c in matrix.columns
                  if c.startswith("MCR") or c.startswith("mcr")]
    other_genes = [c for c in matrix.columns
                   if c not in carb_genes + esbl_genes + flq_genes + mcr_genes]
    col_order = carb_genes + esbl_genes + flq_genes + mcr_genes + other_genes
    matrix = matrix[col_order]

    # --- Plot ---
    fig_height = max(8, len(order) * 0.45)
    fig_width  = max(14, len(col_order) * 0.28)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Custom colour map: white (absent) → dark red (present)
    cmap = sns.color_palette(["#f7f7f7", "#c0392b"], as_cmap=True)

    sns.heatmap(
        matrix, ax=ax,
        cmap=cmap,
        linewidths=0.4,
        linecolor="#dddddd",
        cbar=False,
        xticklabels=True,
        yticklabels=True,
    )

    # Style axes
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90,
                       fontsize=7, fontfamily="DejaVu Sans")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0,
                       fontsize=8, fontfamily="DejaVu Sans")
    ax.set_xlabel("Resistance Gene", fontsize=10, labelpad=10)
    ax.set_ylabel("Isolate", fontsize=10, labelpad=10)
    ax.set_title(
        "Resistome Profile of K. pneumoniae Complex Isolates\n"
        "from Southeast Asia (n=20)",
        fontsize=12, fontweight="bold", pad=12
    )

    # Country colour strip
    add_country_annotation(ax, samples_df, order, country_colors, fig_width)

    # Gene category separator lines
    offsets = [len(carb_genes), len(carb_genes) + len(esbl_genes),
               len(carb_genes) + len(esbl_genes) + len(flq_genes)]
    for x in offsets:
        if 0 < x < len(col_order):
            ax.axvline(x=x, color="black", linewidth=1.5, linestyle="--", alpha=0.5)

    # Country legend
    legend_patches = [
        mpatches.Patch(color=c, label=k)
        for k, c in country_colors.items()
    ]
    ax.legend(handles=legend_patches, title="Country of Origin",
              loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=True, fontsize=8, title_fontsize=9)

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[amr_heatmap] Saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
