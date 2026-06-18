# WGS and Downstream Analysis Pipeline for Southeast Asian KpSC Clinical Isolates

This repository contains a reproducible Snakemake workflow for processing raw whole-genome sequencing (WGS) data of *Klebsiella pneumoniae* species complex (KpSC) clinical isolates from Southeast Asia (Indonesia, Malaysia, Thailand, and Vietnam).

## Features
- **End-to-End Automation:** From raw SRA read retrieval to variant calling, assembly, typing, AMR profiling, phylogeny, and downstream statistical/visual analyses.
- **Storage Optimization:** Automatically deletes intermediate raw FASTQ and large BAM alignment files after processing using Snakemake's `temp()` directive.
- **Reproducibility:** Conda/Mamba environments manage all software versions.
- **Quality Control:** FastQC, fastp (read trimming), Qualimap (alignment stats), QUAST (assembly quality), and MultiQC.
- **Dual Variant Calling:** FreeBayes (+ SnpEff variant annotation) and Snippy (optimized for bacterial genomes).
- **Phylogeny:** Core-genome SNP alignment (via `snippy-core`) and phylogenetic tree construction (FastTree).
- **Genotyping & AMR Profiling:** Kleborate (MLST, K/O-locus typing, virulence scoring) and ABRicate (CARD, ResFinder, NCBI databases).

## Directory Structure
```text
kpneumo-wgs/
├── config/
│   ├── config.yaml          # Pipeline configurations
│   └── samples.tsv          # Metadata and SRA accessions for the 20 isolates
├── data/                    # Reference genome and databases (ignored by Git)
├── results/                 # Primary pipeline outputs (ignored by Git)
├── logs/                    # Rule logs (ignored by Git)
├── scripts/
│   └── setup_databases.sh   # Downloads and configures reference databases
├── envs/                    # Conda environment definition files
├── downstream/              # Downstream analysis pipeline
│   ├── Snakefile            # Downstream workflow definition
│   ├── config.yaml          # Downstream configurations
│   ├── envs/                # Conda environments for downstream tools
│   ├── scripts/             # Downstream R/Python analysis and plotting scripts
│   └── README.md            # Readme for downstream analysis
├── Snakefile                # Main variant calling and assembly workflow
├── environment.yml          # Core Conda environment specification
└── snpEff.config            # SnpEff database configuration
```

## Dataset
This study analyzes 20 clinical isolates from the *Klebsiella pneumoniae* species complex:
- 🇮🇩 **Indonesia (5):** SRR31897979, SRR31897980, SRR31897981, SRR31897982, SRR21679075
- 🇲🇾 **Malaysia (5):** ERR9538944, ERR9538945, ERR9538947, ERR9538952, ERR9538958
- 🇹🇭 **Thailand (5):** SRR26412217, SRR26412242, SRR26412243, SRR26412254, SRR26412293
- 🇻🇳 **Vietnam (5):** SRR6208298, SRR6208299, SRR6208300, SRR6208301, SRR6208302

## Installation & Setup

### 1. Clone the Repository
```bash
git clone https://github.com/engkinandatama/kpneumo-wgs.git
cd kpneumo-wgs
```

### 2. Configure the Conda Environment
Ensure Conda or Mamba is installed. Mamba is recommended for faster dependency resolution.
```bash
conda env create -f environment.yml
conda activate kpneumo_wgs
```

### 3. Setup Reference Databases
Execute the setup script to download reference databases (NCBI HS11286 reference genome and typing databases):
```bash
bash scripts/setup_databases.sh
```

## Execution

### Dry-run Validation
Verify rules and dependencies without executing the pipeline:
```bash
snakemake -n --cores 16
```

### Run Pipeline
Execute the full variant calling and assembly pipeline:
```bash
snakemake --use-conda --cores 16
```

## Main Outputs
- `results/qc/multiqc_report.html`: Interactive summary of quality metrics across all steps.
- `results/variants/{sample}.highmod.vcf`: High/moderate-impact SNPs and indels annotated using SnpEff.
- `results/variants/{sample}_snippy/snps.tab`: Snippy tabular SNP output.
- `results/phylogeny/core.tree`: FastTree core-genome SNP phylogenetic tree (ready for microreact.org or iTOL).
- `results/typing/{sample}_kleborate.txt`: Kleborate output (MLST, K/O antigens, virulence/resistance scores).
- `results/amr/summary_abricate.tab`: Matrix representing resistance genes across all samples.

## Downstream Analysis
Once the main variant calling and assembly steps are complete, proceed to the [downstream directory](file:///e:/Project/kpnumo-wgs/downstream/README.md) for functional annotation (Prokka), pan-genome profiling (Roary), plasmid replicon typing (PlasmidFinder), taxonomy validation (FastANI), and script-based statistical/visual analysis.
