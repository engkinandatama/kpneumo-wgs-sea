# K. pneumoniae WGS Pipeline (SEA Study)

Pipeline Snakemake untuk analisis Whole Genome Sequencing (WGS) isolat *Klebsiella pneumoniae* dari Asia Tenggara (Indonesia, Malaysia, Vietnam).

## Fitur
- **Otomatisasi:** Dari download SRA hingga report final.
- **Efisiensi Penyimpanan:** Menggunakan fitur `temp()` Snakemake untuk menghapus file FASTQ dan BAM mentah setelah diproses.
- **Reproduksibilitas:** Environment dikelola via Conda/Mamba.
- **Analisis Komprehensif:** QC, Variant Calling, Assembly, AMR Profiling (Kleborate, ABRicate), dan MultiQC.

## Struktur Folder
- `config/`: Konfigurasi pipeline dan daftar sampel.
- `data/`: Data mentah dan referensi (diabaikan oleh git).
- `results/`: Hasil analisis (diabaikan oleh git).
- `logs/`: Log file setiap step.
- `scripts/`: Helper scripts.

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
Jalankan script setup untuk mendownload database AMR dan referensi genome.
```bash
bash scripts/setup_databases.sh
```

### 4. Run Pipeline
Lakukan *dry-run* terlebih dahulu untuk memastikan semua rule terhubung:
```bash
snakemake -n
```

Jalankan pipeline dengan batas 32 core (1/4 kapasitas HPC) agar tetap ramah bagi pengguna lain:
```bash
snakemake --cores 32 --resources mem_mb=100000
```

## Modifikasi Sampel
Daftar sampel dikelola di `config/samples.tsv`. Anda bisa menambah atau mengurangi ID aksesion SRA di sana.
