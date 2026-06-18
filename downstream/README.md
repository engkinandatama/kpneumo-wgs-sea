# Downstream Genomic and Statistical Analysis Pipeline

This downstream analysis sub-workflow performs functional annotation, taxonomic validation, plasmid typing, pan-genome comparison, genomic distance profiling, and statistical modeling on the processed isolates.

## Analyses Performed

| # | Rule | Output | Scientific Purpose |
|---|---|---|---|
| 1 | `prokka` | `results/downstream/prokka/` | Functional annotation of genes across draft assemblies |
| 2 | `fastani_compute` | `results/downstream/fastani/ani_matrix.tsv` | Speciation verification using average nucleotide identity (ANI) |
| 3 | `plasmidfinder` | `results/downstream/plasmid/` | Identifies plasmid replicons carrying resistance genes |
| 4 | `roary` | `results/downstream/roary/` | Pan-genome analysis: core vs. accessory gene determination |
| 5 | `amr_heatmap` | `results/downstream/figures/amr_heatmap.pdf` | Generates clustered resistome heatmaps for the publication |
| 6 | `snp_distance` | `results/downstream/snp_distance/` | Computes pairwise core SNP distances between strains |
| 7 | `statistical_tests` | `results/downstream/statistics/stats_report.txt` | Calculates Fisher's Exact, Chi-Square, and Kruskal-Wallis tests |
| 8 | `geo_map` | `results/downstream/figures/geo_map.pdf` | Plots regional distribution and metadata maps |

## Prerequisites

Ensure that the main Snakemake workflow (defined in the root directory) has finished and generated:
- `results/amr/{sample}_assembly/contigs.fasta` for all 20 samples.
- `results/phylogeny/core.tab` (from Snippy-core).
- `results/amr/summary_abricate.tab`.
- `results/ready_to_download/typing/metadata_summary.tsv`.

## Usage

Run the downstream Snakemake sub-workflow by specifying the downstream Snakefile path:

### 1. Dry-run Validation
```bash
snakemake -s downstream/Snakefile --use-conda -n --cores 16
```

### 2. Generate DAG Diagram (Optional)
```bash
snakemake -s downstream/Snakefile --dag | dot -Tpdf > downstream/dag.pdf
```

### 3. Run Pipeline
```bash
snakemake -s downstream/Snakefile --use-conda --cores 16
```

## Output Structure

```text
results/downstream/
├── prokka/
│   └── {sample}/                # Prokka GFF, GBK, FAA, and FNA annotations
├── fastani/
│   ├── ani_matrix.tsv           # Pairwise FastANI matrix
│   └── ani_classified.tsv       # Filtered ANI species designations
├── plasmid/
│   ├── {sample}/                # PlasmidFinder outputs per sample
│   └── plasmid_summary.tsv      # Aggregated plasmid replicon profiles
├── roary/
│   ├── gene_presence_absence.csv
│   └── summary_statistics.txt   # Core genome vs. pan-genome stats
├── snp_distance/
│   ├── snp_matrix.tsv           # Core genome pairwise SNP distance matrix
│   ├── snp_heatmap.pdf          # Log-scaled SNP heatmap (Publication quality)
│   └── snp_heatmap.png          # Log-scaled SNP heatmap (Preview)
├── statistics/
│   └── stats_report.txt         # Statistical analysis outputs (scipy.stats)
└── figures/
    ├── amr_heatmap.pdf          # Publication Jaccard-clustered resistome heatmap
    ├── amr_heatmap.png          # Preview Jaccard-clustered resistome heatmap
    ├── geo_map.pdf              # Regional map
    └── geo_map.png              # Regional map preview
```

## Important Notes
- **Taxonomic Diversity:** The four Indonesian isolates are confirmed as *K. quasipneumoniae subsp. quasipneumoniae* (ANI ~93.8% against *K. pneumoniae* HS11286 reference genome). They are included in pan-genome analysis to evaluate KpSC-wide gene diversity.
- **PlasmidFinder:** The workflow automatically handles database setup and checks for conserved plasmid types like IncFIB, IncFII, and IncX4 that are associated with carbapenemase mobilization.
- **Output Formats:** Figures are produced in **PDF** format (vector graphics for publication) and **PNG** format (for fast rendering/previews).
