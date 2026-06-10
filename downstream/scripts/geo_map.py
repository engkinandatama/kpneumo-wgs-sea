import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Headless backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import geopandas as gpd
import os

# ============================================================
# Geographic Distribution Map – Southeast Asia
# Shows isolate count per country with pie charts of species
# Input : config/samples.tsv
#         results/downstream/fastani/ani_classified.tsv  (species labels)
# Output: results/downstream/figures/geo_map.pdf
#         results/downstream/figures/geo_map.png
# ============================================================

# Southeast Asian countries of interest (ISO codes)
SEA_COUNTRIES = {
    "Indonesia": "IDN",
    "Malaysia":  "MYS",
    "Vietnam":   "VNM",
    "Thailand":  "THA",
}

# Country centroid for bubble placement (lon, lat)
COUNTRY_CENTROIDS = {
    "Indonesia": (118.0, -2.5),
    "Malaysia":  (109.5, 4.5),
    "Vietnam":   (106.0, 16.0),
    "Thailand":  (101.0, 15.5),
}

# Unified publication-grade country palette (matching heatmaps)
COUNTRY_COLORS = {
    "Indonesia": "#2c7bb6",  # Elegant Blue
    "Malaysia":  "#1a9641",  # Elegant Green
    "Vietnam":   "#d7191c",  # Elegant Red
    "Thailand":  "#fdae61",  # Elegant Orange
}

# Unified species palette
SPECIES_COLORS = {
    "Klebsiella pneumoniae":       "#252525",  # Dark Charcoal
    "Klebsiella quasipneumoniae":  "#41b6c4",  # Muted Teal
    "Unknown":                     "#94a3b8",  # Slate
}


def load_species_labels(fastani_tsv: str, samples_df: pd.DataFrame) -> pd.Series:
    """
    Parse FastANI result to assign species labels.
    If FastANI TSV not available, default to metadata column.
    """
    try:
        ani = pd.read_csv(fastani_tsv, sep="\t",
                          names=["query", "reference", "ANI", "mapped", "total"])
        # query is the sample contigs path → extract sample_id
        ani["sample_id"] = ani["query"].apply(
            lambda x: x.split("/")[-2] if "/" in x else x.replace("_contigs.fasta", "")
        )
        # reference path → extract species label
        ani["species_label"] = ani["reference"].apply(
            lambda x: (
                "Klebsiella quasipneumoniae" if "quasipneumoniae" in x.lower()
                else "Klebsiella pneumoniae"
            )
        )
        # Best ANI hit per sample (highest ANI)
        best = (
            ani.sort_values("ANI", ascending=False)
               .drop_duplicates("sample_id")
               .set_index("sample_id")["species_label"]
         )
        return best
    except Exception:
        # Fallback: assume Kp for all
        return pd.Series("Klebsiella pneumoniae",
                         index=samples_df.index, name="species_label")


def main():
    global snakemake
    if 'snakemake' not in globals():
        class MockSnakemake:
            input = type('Input', (), {
                'samples': 'config/samples.tsv',
                'fastani': 'results/downstream/fastani/ani_classified.tsv'
            })()
            output = type('Output', (), {
                'pdf': 'results/downstream/figures/geo_map.pdf',
                'png': 'results/downstream/figures/geo_map.png'
            })()
        snakemake = MockSnakemake()

    samples_path  = snakemake.input.samples
    fastani_path  = snakemake.input.fastani
    out_pdf       = snakemake.output.pdf
    out_png       = snakemake.output.png

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")
    species    = load_species_labels(fastani_path, samples_df)
    samples_df["species"] = samples_df.index.map(species).fillna("Unknown")

    # Count by country × species
    summary = (
        samples_df.groupby(["country", "species"])
                  .size()
                  .reset_index(name="count")
    )
    total_per_country = samples_df.groupby("country").size()

    # Load world map
    try:
        world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    except (AttributeError, ValueError):
        # Geopandas >= 1.0 fallback
        url = "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip"
        world = gpd.read_file(url)
        # Rename columns to match old naturalearth_lowres
        rename_mapper = {
            "POP_EST": "pop_est", 
            "CONTINENT": "continent", 
            "ADMIN": "name", 
            "ADM0_A3": "iso_a3", 
            "GDP_MD": "gdp_md_est"
        }
        world.rename(columns=rename_mapper, inplace=True)

    sea_iso = list(SEA_COUNTRIES.values())
    sea_map = world[world["iso_a3"].isin(sea_iso)].copy()

    # Set publication-quality font
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "figure.dpi": 300,
        "axes.linewidth": 0.8,
        "patch.linewidth": 0.8,
    })

    fig, ax = plt.subplots(1, 1, figsize=(9, 7), facecolor="white")
    
    # Plot global background map
    world.plot(ax=ax, color="#f8fafc", edgecolor="#e2e8f0", linewidth=0.5)
    # Plot Southeast Asia region background
    sea_map.plot(ax=ax, color="#f1f5f9", edgecolor="#cbd5e1", linewidth=0.8)

    # Highlight target countries
    for country, iso in SEA_COUNTRIES.items():
        country_shape = sea_map[sea_map["iso_a3"] == iso]
        if country_shape.empty:
            continue
        # Subtly shade target countries using their assigned country palette color
        country_shape.plot(
            ax=ax,
            color=COUNTRY_COLORS.get(country, "#cccccc"),
            edgecolor="#334155",
            linewidth=0.8,
            alpha=0.25,
        )

    # Draw pie charts at centroids
    max_count = total_per_country.max()
    for country, (lon, lat) in COUNTRY_CENTROIDS.items():
        if country not in total_per_country.index:
            continue
        total = total_per_country[country]

        sub = summary[summary["country"] == country]
        sp_counts = dict(zip(sub["species"], sub["count"]))

        # Pie chart inset slices
        sizes  = [sp_counts.get(sp, 0) for sp in ["Klebsiella pneumoniae", "Klebsiella quasipneumoniae", "Unknown"]]
        colors = [SPECIES_COLORS[sp] for sp in ["Klebsiella pneumoniae", "Klebsiella quasipneumoniae", "Unknown"]]
        valid  = [(s, c) for s, c in zip(sizes, colors) if s > 0]
        if valid:
            s_vals, c_vals = zip(*valid)
        else:
            s_vals, c_vals = [1], ["#cccccc"]

        # Place inset pie chart using inset_axes (automatically handles projections/transforms)
        ax_inset = inset_axes(
            ax, 
            width=0.45,  # Width in inches
            height=0.45, # Height in inches
            loc="center",
            bbox_to_anchor=(lon, lat),
            bbox_transform=ax.transData,
            borderpad=0
        )
        
        ax_inset.pie(
            s_vals, 
            colors=c_vals, 
            startangle=90,
            wedgeprops={"linewidth": 0.5, "edgecolor": "white"}
        )
        # Draw thin border around the pie chart
        circle = plt.Circle((0,0), 1.0, fill=False, edgecolor="#334155", linewidth=0.6)
        ax_inset.add_artist(circle)
        ax_inset.set_aspect("equal")

        # Label country and sample size below the bubble
        ax.annotate(
            f"n={total}",
            xy=(lon, lat - 3.2),
            ha="center", va="center",
            fontsize=8.5, fontweight="bold",
            color=COUNTRY_COLORS.get(country, "#1e293b"),
        )
        ax.annotate(
            country,
            xy=(lon, lat - 5.0),
            ha="center", va="center",
            fontsize=8, color="#475569",
            fontweight="semibold"
        )

    # Set map extent for Southeast Asia focus
    ax.set_xlim(95, 135)
    ax.set_ylim(-12, 28)
    ax.set_axis_off()

    ax.set_title(
        "Geographic Distribution of\n"
        "K. pneumoniae Complex Isolates in Southeast Asia (n=20)",
        fontsize=11, fontweight="bold", pad=12, color="#1e293b"
    )

    # Neatly styled legends
    species_patches = [
        mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.5, label=sp.replace("Klebsiella ", "K. "))
        for sp, c in SPECIES_COLORS.items()
        if sp != "Unknown"
    ]
    country_patches = [
        mpatches.Patch(facecolor=c, edgecolor="black", linewidth=0.5, label=k)
        for k, c in COUNTRY_COLORS.items()
    ]

    legend1 = ax.legend(
        handles=species_patches, 
        title="Species", 
        loc="lower left", 
        fontsize=7.5, 
        title_fontsize=8,
        frameon=True,
        facecolor="white",
        edgecolor="#cbd5e1"
    )
    ax.add_artist(legend1)
    
    ax.legend(
        handles=country_patches, 
        title="Country of Origin", 
        loc="lower right", 
        fontsize=7.5, 
        title_fontsize=8,
        frameon=True,
        facecolor="white",
        edgecolor="#cbd5e1"
    )

    # Save output figures
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[geo_map] Saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
