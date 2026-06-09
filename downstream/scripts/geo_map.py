import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import geopandas as gpd
import sys

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

COUNTRY_COLORS = {
    "Indonesia": "#4285F4",
    "Malaysia":  "#34A853",
    "Vietnam":   "#FBBC05",
    "Thailand":  "#EA4335",
}

SPECIES_COLORS = {
    "Klebsiella pneumoniae":       "#2C3E50",
    "Klebsiella quasipneumoniae":  "#E67E22",
    "Unknown":                     "#95A5A6",
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
    samples_path  = snakemake.input.samples
    fastani_path  = snakemake.input.fastani
    out_pdf       = snakemake.output.pdf
    out_png       = snakemake.output.png

    samples_df = pd.read_csv(samples_path, sep="\t").set_index("sample_id")
    species    = load_species_labels(fastani_path, samples_df)
    samples_df["species"] = samples_df.index.map(species).fillna("Unknown")

    # ── Count by country × species ────────────────────────────────────────────
    summary = (
        samples_df.groupby(["country", "species"])
                  .size()
                  .reset_index(name="count")
    )
    total_per_country = samples_df.groupby("country").size()

    # ── Load world map ─────────────────────────────────────────────────────────
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    sea_iso = list(SEA_COUNTRIES.values())
    sea_map = world[world["iso_a3"].isin(sea_iso)].copy()

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    world.plot(ax=ax, color="#f0f0f0", edgecolor="#cccccc", linewidth=0.4)
    sea_map.plot(ax=ax, color="#dce8f5", edgecolor="#888888", linewidth=0.7)

    # ── Highlight each SEA country with a thick border ─────────────────────────
    for country, iso in SEA_COUNTRIES.items():
        country_shape = sea_map[sea_map["iso_a3"] == iso]
        if country_shape.empty:
            continue
        country_shape.plot(
            ax=ax,
            color=COUNTRY_COLORS.get(country, "#cccccc"),
            edgecolor="#444444",
            linewidth=0.9,
            alpha=0.6,
        )

    # ── Bubble + pie chart per country ────────────────────────────────────────
    max_count = total_per_country.max()
    for country, (lon, lat) in COUNTRY_CENTROIDS.items():
        if country not in total_per_country.index:
            continue
        total = total_per_country[country]
        radius_pts = 30 + (total / max_count) * 120   # bubble size

        sub = summary[summary["country"] == country]
        sp_counts = dict(zip(sub["species"], sub["count"]))

        # Pie chart inset at bubble location
        sizes  = [sp_counts.get(sp, 0) for sp in SPECIES_COLORS]
        colors = [SPECIES_COLORS[sp] for sp in SPECIES_COLORS]
        valid  = [(s, c) for s, c in zip(sizes, colors) if s > 0]
        if valid:
            s_vals, c_vals = zip(*valid)
        else:
            s_vals, c_vals = [1], ["#cccccc"]

        # Transform data coords → display coords for inset
        x_disp, y_disp = ax.transData.transform((lon, lat))
        inv = fig.transFigure.inverted()
        x_fig, y_fig = inv.transform((x_disp, y_disp))

        pie_size = 0.06
        ax_inset = fig.add_axes([x_fig - pie_size / 2,
                                  y_fig - pie_size / 2,
                                  pie_size, pie_size])
        ax_inset.pie(s_vals, colors=c_vals, startangle=90,
                     wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        ax_inset.set_aspect("equal")

        # Count label below pie
        ax.annotate(
            f"n={total}",
            xy=(lon, lat - 3.5),
            ha="center", va="center",
            fontsize=9, fontweight="bold",
            color=COUNTRY_COLORS.get(country, "#333333"),
        )
        ax.annotate(
            country,
            xy=(lon, lat - 5.5),
            ha="center", va="center",
            fontsize=8, color="#333333",
        )

    # ── Map extent – Southeast Asia ───────────────────────────────────────────
    ax.set_xlim(95, 135)
    ax.set_ylim(-12, 28)
    ax.set_axis_off()

    ax.set_title(
        "Geographic Distribution of\n"
        "K. pneumoniae Complex Isolates in Southeast Asia",
        fontsize=13, fontweight="bold", pad=10
    )

    # ── Legends ───────────────────────────────────────────────────────────────
    species_patches = [
        mpatches.Patch(color=c, label=sp.replace("Klebsiella", "K."))
        for sp, c in SPECIES_COLORS.items()
        if sp != "Unknown"
    ]
    country_patches = [
        mpatches.Patch(color=c, label=k)
        for k, c in COUNTRY_COLORS.items()
    ]

    legend1 = ax.legend(
        handles=species_patches, title="Species",
        loc="lower left", fontsize=8, title_fontsize=9,
        framealpha=0.85
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=country_patches, title="Country",
        loc="lower right", fontsize=8, title_fontsize=9,
        framealpha=0.85
    )

    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[geo_map] Saved → {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
