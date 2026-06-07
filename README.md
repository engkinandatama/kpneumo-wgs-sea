# K. pneumoniae WGS Pipeline — SEA AMR Study

Pipeline Snakemake untuk analisis Whole Genome Sequencing (WGS) isolat *Klebsiella pneumoniae* dari Asia Tenggara (Indonesia, Malaysia, Vietnam).

## Fitur
- **Otomatisasi penuh:** Dari download SRA hingga report final dan phylogenetic tree.
- **Efisiensi Penyimpanan:** Menggunakan fitur `temp()` Snakemake untuk menghapus file FASTQ dan BAM mentah setelah diproses.
- **Reproduksibilitas:** Environment dikelola via Conda/Mamba dengan versi yang terpinned.
- **QC Komprehensif:** FastQC, fastp, Qualimap (BAM coverage), QUAST (assembly), MultiQC.
- **Dual Variant Calling:** FreeBayes (+ SnpEff + SnpSift HIGH/MODERATE filter) & Snippy (bacterial-optimized).
- **Phylogeny:** Core genome SNP alignment (snippy-core) → FastTree (GTR) → siap upload ke Microreact/iTOL.
- **AMR Profiling Komprehensif:** Kleborate, ABRicate (CARD), AMRFinderPlus (NCBI), ResFinder (plasmid-mediated).
- **Typing:** Kleborate (MLST + K/O-locus + virulence score) + standalone MLST (Pasteur scheme).

## Struktur Folder
```
kpneumo-wgs-sea/
├── config/
│   ├── config.yaml     # Konfigurasi pipeline
│   └── samples.tsv     # Daftar 15 sampel SEA
├── data/               # Data & referensi (diabaikan git)
├── results/            # Output pipeline (diabaikan git)
├── logs/               # Log setiap step (diabaikan git)
├── scripts/
│   └── setup_databases.sh  # Setup database awal
├── Snakefile           # Workflow utama
└── environment.yml     # Conda environment
```

## Dataset
15 isolat *K. pneumoniae* complex dari Asia Tenggara:
- 🇮🇩 **Indonesia (5):** SRR31897984, SRR31897983, SRR31897982, SRR31897981, SRR21679075
- 🇲🇾 **Malaysia (5):** SRR7964123–SRR7964127
- 🇻🇳 **Vietnam (5):** DRR076141–DRR076145

## Cara Penggunaan (HPC)

### 1. Clone Repository
```bash
git clone https://github.com/engkinandatama/kpneumo-wgs-sea.git
cd kpneumo-wgs-sea
```

### 2. Setup Environment
```bash
conda env create -f environment.yml
conda activate kpneumo_wgs
```

### 3. Setup Database & Reference
Jalankan sekali untuk mendownload semua database yang diperlukan.
```bash
bash scripts/setup_databases.sh
```

### 4. Run Pipeline
Lakukan *dry-run* terlebih dahulu untuk memastikan semua rule terhubung dengan benar:
```bash
snakemake -n --cores 64 --resources mem_mb=60000
```

Jalankan pipeline dengan resource limit (64 core, 60 GB RAM = 1/2 kapasitas HPC):
```bash
snakemake --cores 64 --resources mem_mb=60000
```

> **Resource Policy:** Pipeline dikonfigurasi untuk menggunakan maksimal **64 core** dan **60 GB RAM** sekaligus. Snakemake akan otomatis mengatur batching berdasarkan resource ini.
> - SPAdes (assembly): max **3 sampel paralel** (bottleneck RAM: 3×16GB=48GB)
> - Alignment/Snippy: max **7 sampel paralel** (7×8GB=56GB)
> - AMR tools: max **15+ sampel paralel** (ringan, 2-4GB)


### 5. Visualisasi Phylogenetic Tree (Opsional)
Upload `results/phylogeny/core.tree` ke [Microreact](https://microreact.org) atau [iTOL](https://itol.embl.de) bersama metadata `config/samples.tsv` untuk peta interaktif.

## Output Utama
| File | Deskripsi |
|---|---|
| `results/qc/multiqc_report.html` | Agregasi semua QC report |
| `results/variants/{sample}.highmod.vcf` | Variant HIGH/MODERATE impact (FreeBayes + SnpEff + SnpSift) |
| `results/variants/{sample}_snippy/snps.tab` | SNP summary per sampel (Snippy) |
| `results/phylogeny/core.tree` | Phylogenetic tree (FastTree GTR) |
| `results/typing/{sample}_kleborate.txt` | MLST, K/O-locus, virulence score |
| `results/amr/summary_abricate.tab` | Summary AMR genes semua sampel |

## Modifikasi Sampel
Daftar sampel dikelola di `config/samples.tsv`. Tambah atau kurangi ID aksesion SRA di sana, lalu jalankan ulang pipeline.
