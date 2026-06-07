#!/bin/bash
# Script untuk setup database awal di HPC
# Jalankan SEKALI sebelum pipeline pertama kali dijalankan.
# Usage: bash scripts/setup_databases.sh
set -euo pipefail

echo "=============================================="
echo " Database Setup: K. pneumoniae WGS Pipeline  "
echo "=============================================="

mkdir -p data/reference

# ─── 0. Create Conda Sub-environments ─────────────────────────────────────────
echo "[0/6] Checking and creating/updating Conda sub-environments..."
for env_file in envs/*.yaml; do
    env_name=$(grep "^name:" "$env_file" | cut -d' ' -f2)
    if ! conda env list | grep -q "^$env_name "; then
        echo "Creating environment $env_name..."
        mamba env create -f "$env_file"
    else
        echo "Environment $env_name already exists, updating..."
        mamba env update -f "$env_file"
    fi
done

# ─── 1. Reference Genome ──────────────────────────────────────────────────────
REF="data/reference/Kpneumoniae_HS11286.fa"
if [ ! -f "$REF" ]; then
    echo "[1/6] Downloading K. pneumoniae HS11286 reference genome..."
    wget --tries=3 -q \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz" \
        -O ${REF}.gz
    gunzip ${REF}.gz
    conda run -n kpneumo_align bwa index $REF
    conda run -n kpneumo_align samtools faidx $REF
    echo "    Reference genome ready: $REF"
else
    echo "[1/6] Reference genome already exists, skipping."
fi

# ─── 2. SnpEff Database ───────────────────────────────────────────────────────
SNPEFF_DB_DIR="data/snpeff_data/Klebsiella_pneumoniae_subsp_pneumoniae_HS11286"
if [ ! -f "$SNPEFF_DB_DIR/genes.snpEffectPredictor.bin" ]; then
    echo "[2/6] Building local SnpEff database for K. pneumoniae HS11286..."
    mkdir -p "$SNPEFF_DB_DIR"
    
    # Copy reference FASTA
    cp "$REF" "$SNPEFF_DB_DIR/sequences.fa"
    
    # Download GFF3
    wget --tries=3 -q \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.gff.gz" \
        -O "$SNPEFF_DB_DIR/genes.gff.gz"
        
    # Build database
    conda run -n kpneumo_variants snpEff build -gff3 -noCheckCds -noCheckProtein -v Klebsiella_pneumoniae_subsp_pneumoniae_HS11286
    echo "    SnpEff database built successfully."
else
    echo "[2/6] SnpEff database already exists, skipping."
fi

# ─── 3. AMRFinderPlus Database ────────────────────────────────────────────────
echo "[3/6] Updating AMRFinderPlus database..."
mkdir -p data/reference/amrfinder_db
conda run -n kpneumo_amr_typing amrfinder_update --database data/reference/amrfinder_db
echo "    AMRFinderPlus database ready."

# ─── 4. CARD Database (for ABRicate & RGI) ────────────────────────────────────
echo "[4/6] Updating ABRicate CARD database..."
conda run -n kpneumo_amr_typing abricate --setupdb
echo "    ABRicate databases ready."

# ─── 5. ResFinder Database ────────────────────────────────────────────────────
RESFINDER_DB="data/reference/resfinder_db"
if [ -d "$RESFINDER_DB" ]; then
    REMOTE_URL=$(git -C "$RESFINDER_DB" remote get-url origin 2>/dev/null || echo "")
    if [[ "$REMOTE_URL" != *"bitbucket.org"* ]]; then
        echo "Old ResFinder database remote detected. Re-cloning..."
        rm -rf "$RESFINDER_DB"
    fi
fi

if [ ! -d "$RESFINDER_DB" ]; then
    echo "[5/6] Cloning ResFinder database..."
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git "$RESFINDER_DB"
    cd "$RESFINDER_DB"
    conda run -n kpneumo_amr_typing python INSTALL.py
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
echo " Next: snakemake -n --use-conda              "
echo "=============================================="
