# Downstream Analysis Pipeline

Pipeline Snakemake untuk analisis lanjutan (downstream) pasca pipeline WGS utama
selesai. Dirancang untuk menghasilkan data berkualitas publikasi pada jurnal Q1/Q2.

## Analisis yang Dilakukan

| # | Rule | Output | Tujuan Ilmiah |
|---|---|---|---|
| 1 | `prokka` | `results/downstream/prokka/` | Anotasi fungsional gen per isolat |
| 2 | `fastani_compute` | `results/downstream/fastani/ani_matrix.tsv` | Konfirmasi spesiasi ANI >95% |
| 3 | `plasmidfinder` | `results/downstream/plasmid/` | Identifikasi tipe plasmid pembawa AMR |
| 4 | `roary` | `results/downstream/roary/` | Pan-genome: core vs. accessory genome |
| 5 | `amr_heatmap` | `results/downstream/figures/amr_heatmap.pdf` | Visualisasi resistome untuk manuskrip |
| 6 | `snp_distance` | `results/downstream/snp_distance/` | Matriks jarak genomik antar isolat |
| 7 | `statistical_tests` | `results/downstream/statistics/stats_report.txt` | Uji Fisher's Exact, Chi-square, Kruskal-Wallis |
| 8 | `geo_map` | `results/downstream/figures/geo_map.pdf` | Peta distribusi geografis Asia Tenggara |

## Prasyarat

Pastikan **pipeline utama (`Snakefile` di root project)** sudah selesai dan menghasilkan:
- `results/amr/{sample}_assembly/contigs.fasta` untuk semua sampel
- `results/phylogeny/core.tab` (output Snippy-core)
- `results/amr/summary_abricate.tab`
- `results/ready_to_download/typing/metadata_summary.tsv`

## Cara Menjalankan di HPC

### 1. Push ke GitHub (dari lokal)
```bash
git add downstream/
git commit -m "feat(downstream): add downstream analysis pipeline (Prokka, Roary, FastANI, PlasmidFinder, statistics)"
git push
```

### 2. Pull di HPC
```bash
cd /path/to/project
git pull
```

### 3. Dry-run (validasi tanpa eksekusi)
```bash
snakemake -s downstream/Snakefile \
    --use-conda \
    --cores 32 \
    --resources mem_mb=30000 \
    --dry-run
```

### 4. Generate DAG diagram
```bash
snakemake -s downstream/Snakefile \
    --dag | dot -Tpdf > downstream/dag.pdf
```

### 5. Jalankan pipeline penuh
```bash
snakemake -s downstream/Snakefile \
    --use-conda \
    --conda-prefix ~/.conda-envs \
    --cores 32 \
    --resources mem_mb=30000 \
    --rerun-incomplete \
    2>&1 | tee logs/downstream_run.log
```

## Estimasi Waktu

| Analisis | Estimasi (32 cores) |
|---|---|
| Prokka × 20 sampel | 30–40 menit |
| FastANI all-vs-all | 5 menit |
| PlasmidFinder × 20 sampel | 20 menit |
| Roary (setelah Prokka selesai) | 30–60 menit |
| Statistik + Heatmap + Peta | 5–10 menit |
| **Total** | **~1.5 – 2.5 jam** |

## Struktur Output

```
results/downstream/
├── prokka/
│   └── {sample}/
│       ├── {sample}.gff
│       └── {sample}.faa
├── fastani/
│   ├── ani_matrix.tsv
│   └── ani_classified.tsv
├── plasmid/
│   ├── {sample}/results_tab.tsv
│   └── plasmid_summary.tsv
├── roary/
│   ├── gene_presence_absence.csv
│   └── summary_statistics.txt
├── snp_distance/
│   ├── snp_matrix.tsv
│   ├── snp_heatmap.pdf
│   └── snp_heatmap.png
├── statistics/
│   └── stats_report.txt
└── figures/
    ├── amr_heatmap.pdf
    ├── amr_heatmap.png
    ├── geo_map.pdf
    └── geo_map.png
```

## Catatan Penting

- **K. quasipneumoniae** termasuk dalam analisis pan-genome Roary karena merupakan anggota
  *Klebsiella pneumoniae* species complex (KpSC) dan perbandingannya secara ilmiah menarik.
- **PlasmidFinder** membutuhkan database yang akan otomatis diunduh saat pertama kali dijalankan.
- Semua figure dihasilkan dalam format **PDF** (untuk publikasi) dan **PNG** (untuk preview).
