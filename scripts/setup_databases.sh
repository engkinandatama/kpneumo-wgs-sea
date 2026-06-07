#!/bin/bash
# Script untuk setup database awal di HPC
# Jalankan SEKALI sebelum pipeline pertama kali dijalankan.
# Usage: bash scripts/setup_databases.sh
set -euo pipefail

echo "=============================================="
echo " Database Setup: K. pneumoniae WGS Pipeline  "
echo "=============================================="

mkdir -p data/reference

# ─── 1. Reference Genome ──────────────────────────────────────────────────────
REF="data/reference/Kpneumoniae_HS11286.fa"
if [ ! -f "$REF" ]; then
    echo "[1/6] Downloading K. pneumoniae HS11286 reference genome..."
    wget --tries=3 -q \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz" \
        -O ${REF}.gz
    gunzip ${REF}.gz
    bwa index $REF
    samtools faidx $REF
    echo "    Reference genome ready: $REF"
else
    echo "[1/6] Reference genome already exists, skipping."
fi

# ─── 2. SnpEff Database ───────────────────────────────────────────────────────
echo "[2/6] Downloading SnpEff database for K. pneumoniae HS11286..."
snpEff download -v Klebsiella_pneumoniae_subsp_pneumoniae_HS11286
echo "    SnpEff database ready."

# ─── 3. AMRFinderPlus Database ────────────────────────────────────────────────
echo "[3/6] Updating AMRFinderPlus database..."
amrfinder -u
echo "    AMRFinderPlus database ready."

# ─── 4. CARD Database (for ABRicate & RGI) ────────────────────────────────────
echo "[4/6] Updating ABRicate CARD database..."
abricate --setupdb
echo "    ABRicate databases ready."

# ─── 5. ResFinder Database ────────────────────────────────────────────────────
RESFINDER_DB="data/reference/resfinder_db"
if [ ! -d "$RESFINDER_DB" ]; then
    echo "[5/6] Cloning ResFinder database..."
    git clone https://git.sr.ht/~genomicepidemiology/resfinder_db "$RESFINDER_DB"
    cd "$RESFINDER_DB"
    python INSTALL.py
    cd -
    echo "    ResFinder database ready: $RESFINDER_DB"
else
    echo "[5/6] ResFinder database already exists. Updating..."
    git -C "$RESFINDER_DB" pull
    echo "    ResFinder database updated."
fi

# ─── 6. Kraken2 Database (reserved, download jika diperlukan) ─────────────────
# Uncomment block di bawah jika ingin enable contamination check dengan Kraken2
# KRAKEN_DB="data/reference/kraken2_db"
# if [ ! -d "$KRAKEN_DB" ]; then
#     echo "[6/6] Downloading Kraken2 Standard-8 database (~8GB)..."
#     mkdir -p "$KRAKEN_DB"
#     wget --tries=3 -q \
#         https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz \
#         -O "$KRAKEN_DB/k2_standard_08gb.tar.gz"
#     tar -xzvf "$KRAKEN_DB/k2_standard_08gb.tar.gz" -C "$KRAKEN_DB"
#     rm "$KRAKEN_DB/k2_standard_08gb.tar.gz"
#     echo "    Kraken2 database ready: $KRAKEN_DB"
# fi
echo "[6/6] Kraken2 skipped (uncomment in script to enable)."

echo ""
echo "=============================================="
echo " All databases setup complete! Ready to run. "
echo " Next: snakemake -n  (dry-run check)         "
echo "=============================================="
