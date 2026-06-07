#!/bin/bash
# Script untuk setup database awal di HPC
set -e

echo "Starting database setup..."

# 1. AMRFinderPlus update
echo "Updating AMRFinderPlus database..."
amrfinder -u

# 2. SnpEff database download
echo "Downloading SnpEff database for K. pneumoniae..."
snpEff download -v Klebsiella_pneumoniae_subsp_pneumoniae_HS11286

# 3. CARD Database for RGI/ABRicate
echo "Downloading CARD database..."
mkdir -p data/reference/card
cd data/reference/card
wget https://card.mcmaster.ca/latest/data -O card-data.tar.bz2
tar -xjf card-data.tar.bz2
# Jika menggunakan RGI
# rgi load --card_json card.json --local
cd ../../../

# 4. Download Kraken2 Database (Standard-8, ~15GB)
# Ini mencakup Bacteria, Archaea, Virus, dan Human. Cocok untuk cek kontaminasi.
echo "Downloading Kraken2 Standard-8 database..."
mkdir -p data/reference/kraken2_db
cd data/reference/kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
tar -xzvf k2_standard_08gb_20240112.tar.gz
rm k2_standard_08gb_20240112.tar.gz
cd ../../../

# 5. Download Reference Genome
REF="data/reference/Kpneumoniae_HS11286.fa"
if [ ! -f "$REF" ]; then
    echo "Downloading Reference Genome..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz" -O ${REF}.gz
    gunzip ${REF}.gz
    bwa index $REF
    samtools faidx $REF
fi

echo "Setup complete!"
