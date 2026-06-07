import pandas as pd
import os

# Load config and samples
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample_id", drop=False)
SAMPLES = samples_df.index.tolist()

REF = config["reference"]

# ============================================================
# onstart: Ensure all directories exist before any rule runs.
# Snakemake creates dirs for declared outputs/logs automatically,
# but this guard prevents edge cases in run: blocks and nested paths.
# ============================================================
onstart:
    dirs = [
        "data/raw", "data/reference",
        "results/qc", "results/trimmed", "results/aligned",
        "results/variants", "results/amr", "results/typing", "results/phylogeny",
        "logs/download", "logs/qc", "logs/align", "logs/variants",
        "logs/assemble", "logs/amr", "logs/typing", "logs/phylogeny"
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    print("=" * 54)
    print(f"  K. pneumoniae WGS Pipeline — {len(SAMPLES)} samples")
    print(f"  Cores : {config['max_cores']} | RAM limit: {config['max_mem_mb']} MB")
    print(f"  Run   : snakemake --cores {config['max_cores']} --resources mem_mb={config['max_mem_mb']}")
    print("=" * 54)

rule all:
    input:
        # --- QC Reports ---
        "results/qc/multiqc_report.html",
        # --- AMR Summaries ---
        "results/amr/summary_abricate.tab",
        # --- Per-sample: BAM QC ---
        expand("results/qc/{sample}_flagstat.txt", sample=SAMPLES),
        expand("results/qc/{sample}_qualimap", sample=SAMPLES),
        # --- Per-sample: Assembly QC ---
        expand("results/qc/{sample}_quast", sample=SAMPLES),
        # --- Per-sample: Variant Calling (FreeBayes) ---
        expand("results/variants/{sample}.highmod.vcf", sample=SAMPLES),
        # --- Per-sample: Variant Calling (Snippy) ---
        expand("results/variants/{sample}_snippy/snps.vcf", sample=SAMPLES),
        # --- Phylogeny ---
        "results/phylogeny/core.tree",
        # --- Per-sample: Typing ---
        expand("results/typing/{sample}_kleborate.txt", sample=SAMPLES),
        expand("results/typing/{sample}_mlst.txt", sample=SAMPLES),
        # --- Per-sample: AMR Profiling ---
        expand("results/amr/{sample}_amrfinder.txt", sample=SAMPLES),
        expand("results/amr/{sample}_resfinder", sample=SAMPLES)

# ============================================================
# --- 1. DOWNLOAD DATA ---
# Threads: 2 (network-bound, not CPU-bound)
# RAM    : 4 GB | Max parallel: 15x (negligible)
# ============================================================
rule download_sra:
    output:
        r1 = temp("data/raw/{sample}_1.fastq.gz"),
        r2 = temp("data/raw/{sample}_2.fastq.gz")
    log: "logs/download/{sample}.log"
    threads: config["threads"]["download"]
    resources:
        mem_mb = 4000
    shell:
        """
        prefetch {wildcards.sample} -O data/raw/ &> {log}
        fasterq-dump data/raw/{wildcards.sample} -O data/raw/ --split-files --threads {threads} &>> {log}
        gzip data/raw/{wildcards.sample}_1.fastq &>> {log}
        gzip data/raw/{wildcards.sample}_2.fastq &>> {log}
        """

# ============================================================
# --- 2. QC & TRIMMING ---
# fastp  : 4 threads, 8 GB  | Max parallel: 7x (RAM bottleneck)
# fastqc : 2 threads, 2 GB  | Max parallel: 30x
# ============================================================
rule fastp:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"
    output:
        r1   = temp("results/trimmed/{sample}_R1.fastq.gz"),
        r2   = temp("results/trimmed/{sample}_R2.fastq.gz"),
        json = "results/qc/{sample}_fastp.json",
        html = "results/qc/{sample}_fastp.html"
    log: "logs/qc/{sample}_fastp.log"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb = 8000
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -j {output.json} -h {output.html} \
              -q 20 -l 36 --thread {threads} --detect_adapter_for_pe \
              2> {log}
        """

rule fastqc:
    input: "results/trimmed/{sample}_{read}.fastq.gz"
    output: "results/qc/{sample}_{read}_fastqc.html"
    log: "logs/qc/{sample}_{read}_fastqc.log"
    threads: 2
    resources:
        mem_mb = 2000
    shell: "fastqc {input} -o results/qc/ -t {threads} 2> {log}"

# ============================================================
# --- 3. ALIGNMENT & BAM QC ---
# align  : 8 threads, 8 GB  | Max parallel: 7x
# bam_qc : 8 threads, 8 GB  | Max parallel: 7x (combined with align: ~3-4x)
# Note: Qualimap java-mem-size=7G stays within 8GB resource budget
# ============================================================
rule align:
    input:
        r1  = "results/trimmed/{sample}_R1.fastq.gz",
        r2  = "results/trimmed/{sample}_R2.fastq.gz",
        ref = REF
    output:
        sorted_bam = "results/aligned/{sample}.sorted.bam"
    log: "logs/align/{sample}.log"
    threads: config["threads"]["align"]
    resources:
        mem_mb = config["mem"]["align"]
    shell:
        """
        bwa mem -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
            {input.ref} {input.r1} {input.r2} 2>> {log} \
            | samtools sort -@ {threads} -o {output.sorted_bam} 2>> {log}
        samtools index {output.sorted_bam} 2>> {log}
        """

rule bam_qc:
    """Quality gate before variant calling.
    flagstat: alignment rate | qualimap: per-base coverage depth."""
    input: "results/aligned/{sample}.sorted.bam"
    output:
        flagstat = "results/qc/{sample}_flagstat.txt",
        qualimap = directory("results/qc/{sample}_qualimap")
    log: "logs/qc/{sample}_bamqc.log"
    threads: config["threads"]["align"]
    resources:
        mem_mb = 8000
    shell:
        """
        samtools flagstat {input} > {output.flagstat} 2> {log}
        qualimap bamqc -bam {input} -outdir {output.qualimap} \
            --java-mem-size=7G -nt {threads} >> {log} 2>&1
        """

# ============================================================
# --- 4. VARIANT CALLING (FreeBayes + SnpEff + SnpSift) ---
# freebayes : 1 thread, 8 GB | Max parallel: 7x
# filter    : 1 thread, 2 GB | Max parallel: 30x
# snpeff    : 1 thread, 4 GB | Max parallel: 15x
# snpsift   : 1 thread, 4 GB | Max parallel: 15x
# ============================================================
rule variant_calling:
    input:
        bam = "results/aligned/{sample}.sorted.bam",
        ref = REF
    output:
        vcf = "results/variants/{sample}.raw.vcf"
    log: "logs/variants/{sample}_freebayes.log"
    resources:
        mem_mb = 8000
    shell: "freebayes -f {input.ref} {input.bam} > {output.vcf} 2> {log}"

rule filter_variants:
    input:
        vcf = "results/variants/{sample}.raw.vcf"
    output:
        vcf = "results/variants/{sample}.filtered.vcf"
    log: "logs/variants/{sample}_filter.log"
    resources:
        mem_mb = 2000
    shell: "vcffilter -f 'QUAL > 20 & DP > 10' {input.vcf} > {output.vcf} 2> {log}"

rule annotate_variants:
    input:
        vcf = "results/variants/{sample}.filtered.vcf"
    output:
        vcf = "results/variants/{sample}.annotated.vcf"
    log: "logs/variants/{sample}_snpeff.log"
    params:
        db = config["snpeff_db"]
    resources:
        mem_mb = 4000
    shell: "snpEff {params.db} {input.vcf} > {output.vcf} 2> {log}"

rule snpsift_filter:
    """Keep only HIGH and MODERATE impact variants for focused analysis."""
    input:
        vcf = "results/variants/{sample}.annotated.vcf"
    output:
        vcf = "results/variants/{sample}.highmod.vcf"
    log: "logs/variants/{sample}_snpsift.log"
    resources:
        mem_mb = 4000
    shell:
        """
        SnpSift filter \
            "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" \
            {input.vcf} > {output.vcf} 2> {log}
        """

# ============================================================
# --- 5. SNIPPY & PHYLOGENY ---
# snippy      : 8 threads, 8 GB  | Max parallel: 7x
# snippy_core : 1 thread,  8 GB  | Runs once (all samples done)
# phylogeny   : 8 threads, 8 GB  | Runs once after snippy_core
# ============================================================
rule snippy:
    """Snippy: bacterial-specific variant calling. Produces clean per-sample VCF
    and tabular SNP summary for easy downstream analysis."""
    input:
        r1  = "results/trimmed/{sample}_R1.fastq.gz",
        r2  = "results/trimmed/{sample}_R2.fastq.gz",
        ref = REF
    output:
        vcf = "results/variants/{sample}_snippy/snps.vcf",
        tab = "results/variants/{sample}_snippy/snps.tab"
    log: "logs/variants/{sample}_snippy.log"
    threads: config["threads"]["align"]
    resources:
        mem_mb = 8000
    shell:
        """
        snippy --cpus {threads} \
               --outdir results/variants/{wildcards.sample}_snippy \
               --ref {input.ref} \
               --R1 {input.r1} \
               --R2 {input.r2} \
               --force >> {log} 2>&1
        """

rule snippy_core:
    """Compute core genome SNP alignment across all samples.
    core.full.aln includes invariant sites; used as input for FastTree."""
    input:
        expand("results/variants/{sample}_snippy/snps.vcf", sample=SAMPLES)
    output:
        aln  = "results/phylogeny/core.aln",
        full = "results/phylogeny/core.full.aln",
        tab  = "results/phylogeny/core.tab"
    log: "logs/phylogeny/snippy_core.log"
    resources:
        mem_mb = 8000
    params:
        ref  = REF,
        dirs = " ".join([f"results/variants/{s}_snippy" for s in SAMPLES])
    shell:
        """
        snippy-core --ref {params.ref} \
            --prefix results/phylogeny/core \
            {params.dirs} >> {log} 2>&1
        """

rule phylogeny:
    """Build maximum-likelihood phylogenetic tree using FastTree (GTR+CAT model).
    Output .tree file is directly uploadable to Microreact / iTOL."""
    input: "results/phylogeny/core.full.aln"
    output: "results/phylogeny/core.tree"
    log: "logs/phylogeny/fasttree.log"
    threads: config["threads"]["align"]
    resources:
        mem_mb = 8000
    shell:
        "FastTree -gtr -nt -threads {threads} {input} > {output} 2> {log}"

# ============================================================
# --- 6. ASSEMBLY & ASSEMBLY QC ---
# SPAdes : 8 threads, 16 GB | Max parallel: 3x (MEMORY BOTTLENECK!)
# quast  : 4 threads,  4 GB | Max parallel: 15x
# ============================================================
rule assemble:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        contigs = "results/amr/{sample}_assembly/contigs.fasta"
    log: "logs/assemble/{sample}.log"
    threads: config["threads"]["assemble"]
    resources:
        mem_mb = config["mem"]["assemble"]
    run:
        # SPAdes flag -m requires memory in GB, convert from MB
        mem_gb = int(resources.mem_mb / 1000)
        shell(
            "spades.py -1 {input.r1} -2 {input.r2} "
            "-o results/amr/{wildcards.sample}_assembly "
            "--threads {threads} -m " + str(mem_gb) + " >> {log} 2>&1"
        )

rule quast:
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: directory("results/qc/{sample}_quast")
    log: "logs/qc/{sample}_quast.log"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb = 4000
    shell: "quast.py {input} -o {output} --threads {threads} 2> {log}"

# ============================================================
# --- 7. AMR PROFILING & TYPING ---
# All tools: 1 thread, 2-4 GB | Max parallel: 15-30x
# ============================================================
rule kleborate:
    """Kleborate: MLST + K/O loci + virulence + resistance scoring for Klebsiella."""
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: "results/typing/{sample}_kleborate.txt"
    log: "logs/typing/{sample}_kleborate.log"
    resources:
        mem_mb = 4000
    shell: "kleborate --all -a {input} -o {output} 2> {log}"

rule mlst_typing:
    """Standalone MLST cross-validation against Klebsiella scheme (Pasteur)."""
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: "results/typing/{sample}_mlst.txt"
    log: "logs/typing/{sample}_mlst.log"
    resources:
        mem_mb = 2000
    shell: "mlst --scheme klebsiella {input} > {output} 2> {log}"

rule abricate:
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: "results/amr/{sample}_abricate.tab"
    log: "logs/amr/{sample}_abricate.log"
    resources:
        mem_mb = 2000
    shell: "abricate --db card {input} > {output} 2> {log}"

rule amrfinderplus:
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: "results/amr/{sample}_amrfinder.txt"
    log: "logs/amr/{sample}_amrfinder.log"
    resources:
        mem_mb = 4000
    shell: "amrfinder -n {input} -O Klebsiella -o {output} 2> {log}"

rule resfinder:
    """ResFinder: detects plasmid-mediated acquired resistance genes.
    Complements ABRicate/AMRFinder for comprehensive resistance profiling."""
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: directory("results/amr/{sample}_resfinder")
    log: "logs/amr/{sample}_resfinder.log"
    params:
        db = config["resfinder_db"]
    resources:
        mem_mb = 4000
    shell:
        """
        python -m resfinder \
            -ifa {input} \
            -db_res {params.db} \
            -o {output} \
            --acquired >> {log} 2>&1
        """

# ============================================================
# --- 8. SUMMARIES ---
# ============================================================
rule abricate_summary:
    input: expand("results/amr/{sample}_abricate.tab", sample=SAMPLES)
    output: "results/amr/summary_abricate.tab"
    log: "logs/amr/abricate_summary.log"
    resources:
        mem_mb = 2000
    shell: "abricate --summary {input} > {output} 2> {log}"

rule multiqc:
    """Aggregate all QC reports: fastp, FastQC, flagstat, Qualimap, QUAST."""
    input:
        expand("results/qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/qc/{sample}_flagstat.txt", sample=SAMPLES),
        expand("results/qc/{sample}_qualimap", sample=SAMPLES),
        expand("results/qc/{sample}_quast", sample=SAMPLES)
    output: "results/qc/multiqc_report.html"
    log: "logs/qc/multiqc.log"
    resources:
        mem_mb = 4000
    shell: "multiqc results/qc/ -o results/qc/ -n multiqc_report.html 2> {log}"
