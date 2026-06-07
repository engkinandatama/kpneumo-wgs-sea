import pandas as pd

# Load config and samples
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample_id", drop=False)
SAMPLES = samples_df.index.tolist()

REF = config["reference"]

rule all:
    input:
        "results/qc/multiqc_report.html",
        "results/amr/summary_abricate.tab",
        expand("results/typing/{sample}_kleborate.txt", sample=SAMPLES),
        expand("results/variants/{sample}.annotated.vcf", sample=SAMPLES),
        expand("results/amr/{sample}_amrfinder.txt", sample=SAMPLES),
        expand("results/qc/{sample}_quast", sample=SAMPLES)

# --- 1. DOWNLOAD DATA ---
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

# --- 2. QC & TRIMMING ---
rule fastp:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"
    output:
        r1 = temp("results/trimmed/{sample}_R1.fastq.gz"),
        r2 = temp("results/trimmed/{sample}_R2.fastq.gz"),
        json = "results/qc/{sample}_fastp.json",
        html = "results/qc/{sample}_fastp.html"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb = 8000
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -j {output.json} -h {output.html} \
              -q 20 -l 36 --thread {threads} --detect_adapter_for_pe
        """

rule fastqc:
    input: "results/trimmed/{sample}_{read}.fastq.gz"
    output: "results/qc/{sample}_{read}_fastqc.html"
    threads: 2
    resources:
        mem_mb = 2000
    shell: "fastqc {input} -o results/qc/ -t {threads}"

# --- 3. ALIGNMENT & VARIANT CALLING ---
rule align:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        ref = REF
    output:
        sorted_bam = "results/aligned/{sample}.sorted.bam"
    log: "logs/align/{sample}.log"
    threads: config["threads"]["align"]
    resources:
        mem_mb = config["mem"]["align"]
    shell:
        """
        bwa mem -t {threads} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
            {input.ref} {input.r1} {input.r2} 2>> {log} \
            | samtools sort -@ {threads} -o {output.sorted_bam} 2>> {log}
        samtools index {output.sorted_bam} 2>> {log}
        """

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

# --- 4. ASSEMBLY & DOWNSTREAM ---
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

rule kleborate:
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: "results/typing/{sample}_kleborate.txt"
    log: "logs/typing/{sample}_kleborate.log"
    resources:
        mem_mb = 4000
    shell: "kleborate --all -a {input} -o {output} 2> {log}"

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

rule quast:
    input: "results/amr/{sample}_assembly/contigs.fasta"
    output: directory("results/qc/{sample}_quast")
    log: "logs/qc/{sample}_quast.log"
    threads: config["threads"]["fastp"]
    resources:
        mem_mb = 4000
    shell: "quast.py {input} -o {output} --threads {threads} 2> {log}"

# --- 5. SUMMARIES ---
rule abricate_summary:
    input: expand("results/amr/{sample}_abricate.tab", sample=SAMPLES)
    output: "results/amr/summary_abricate.tab"
    shell: "abricate --summary {input} > {output}"

rule multiqc:
    input:
        expand("results/qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"])
    output: "results/qc/multiqc_report.html"
    resources:
        mem_mb = 4000
    shell: "multiqc results/qc/ -o results/qc/ -n multiqc_report.html"
