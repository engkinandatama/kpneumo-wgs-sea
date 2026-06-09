#!/usr/bin/env python3
import os
import json
import re
import pandas as pd

def parse_fastp(sample_id):
    json_path = f"results/qc/{sample_id}_fastp.json"
    if not os.path.exists(json_path):
        return None
    
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    before_reads = data["summary"]["before_filtering"]["total_reads"]
    after_reads = data["summary"]["after_filtering"]["total_reads"]
    q30_rate = data["summary"]["after_filtering"]["q30_rate"]
    gc_content = data["summary"]["after_filtering"]["gc_content"]
    
    passed_trim_pct = (after_reads / before_reads) * 100 if before_reads > 0 else 0.0
    q30_pct = q30_rate * 100
    gc_pct = gc_content * 100
    
    return {
        "raw_reads": before_reads,
        "passed_trim_pct": round(passed_trim_pct, 2),
        "q30_pct": round(q30_pct, 2),
        "gc_pct": round(gc_pct, 2)
    }

def parse_flagstat(sample_id):
    flagstat_path = f"results/qc/{sample_id}_flagstat.txt"
    if not os.path.exists(flagstat_path):
        return None
    
    with open(flagstat_path, 'r') as f:
        content = f.read()
        
    # Standard samtools flagstat mapped percentage pattern: "mapped (98.45% : N/A)"
    match = re.search(r"mapped \(([\d\.]+)%\s*:", content)
    if match:
        return round(float(match.group(1)), 2)
    
    # Alternative pattern: "mapped (98.45%)"
    match_alt = re.search(r"mapped \(([\d\.]+)%\)", content)
    if match_alt:
        return round(float(match_alt.group(1)), 2)
        
    return 0.0

def parse_qualimap(sample_id):
    qualimap_path = f"results/qc/{sample_id}_qualimap/genome_results.txt"
    if not os.path.exists(qualimap_path):
        return 0.0
    
    with open(qualimap_path, 'r') as f:
        for line in f:
            if "mean coverageData" in line:
                # E.g., "     mean coverageData = 55.23X" or "     mean coverageData = 55.23"
                match = re.search(r"=\s*([\d\.]+)", line)
                if match:
                    return round(float(match.group(1)), 2)
    return 0.0

def parse_variants(sample_id):
    vcf_path = f"results/variants/{sample_id}.filtered.vcf"
    if not os.path.exists(vcf_path):
        return None
    
    total_variants = 0
    snps = 0
    indels = 0
    transitions = 0
    transversions = 0
    
    transitions_set = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            total_variants += 1
            fields = line.strip().split("\t")
            ref = fields[3].upper()
            alt = fields[4].upper()
            
            # Simple check for SNPs vs Indels (handles single allele calls)
            if len(ref) == 1 and len(alt) == 1:
                snps += 1
                pair = (ref, alt)
                if pair in transitions_set:
                    transitions += 1
                else:
                    transversions += 1
            else:
                indels += 1
                
    ts_tv = transitions / transversions if transversions > 0 else 0.0
    
    # Parse high and moderate impact from the highmod VCF (containing only HIGH/MODERATE)
    highmod_path = f"results/variants/{sample_id}.highmod.vcf"
    high_impact = 0
    moderate_impact = 0
    
    if os.path.exists(highmod_path):
        with open(highmod_path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                # Search for ANN field in the INFO column (field 8)
                # Format: ANN=Allele|Annotation|Impact|Gene|...
                info = line.strip().split("\t")[7]
                if "|HIGH|" in info:
                    high_impact += 1
                elif "|MODERATE|" in info:
                    moderate_impact += 1
                    
    return {
        "total_filtered_variants": total_variants,
        "snps": snps,
        "indels": indels,
        "ts_tv": round(ts_tv, 2),
        "high_impact": high_impact,
        "moderate_impact": moderate_impact
    }

def main():
    samples_file = "config/samples.tsv"
    if not os.path.exists(samples_file):
        print(f"Error: {samples_file} not found.")
        return
        
    df_samples = pd.read_csv(samples_file, sep="\t")
    
    records = []
    print("Gathering statistics for KpSC WGS variant calling pipeline...")
    
    for _, row in df_samples.iterrows():
        sample_id = row["sample_id"]
        country = row["country"]
        
        print(f"Processing sample: {sample_id}...")
        
        # Parse fastp QC stats
        fastp_stats = parse_fastp(sample_id)
        if fastp_stats is None:
            print(f"  Warning: fastp JSON not found for {sample_id}")
            fastp_stats = {"raw_reads": 0, "passed_trim_pct": 0.0, "q30_pct": 0.0, "gc_pct": 0.0}
            
        # Parse mapping stats
        mapped_pct = parse_flagstat(sample_id)
        if mapped_pct is None:
            print(f"  Warning: flagstat file not found for {sample_id}")
            mapped_pct = 0.0
            
        mean_depth = parse_qualimap(sample_id)
        
        # Parse variant statistics
        vcf_stats = parse_variants(sample_id)
        if vcf_stats is None:
            print(f"  Warning: filtered VCF not found for {sample_id}")
            vcf_stats = {
                "total_filtered_variants": 0,
                "snps": 0,
                "indels": 0,
                "ts_tv": 0.0,
                "high_impact": 0,
                "moderate_impact": 0
            }
            
        records.append({
            "Isolate ID": sample_id,
            "Country": country,
            "Raw Reads": fastp_stats["raw_reads"],
            "Passed Trim (%)": fastp_stats["passed_trim_pct"],
            "Q30 Bases (%)": fastp_stats["q30_pct"],
            "GC Content (%)": fastp_stats["gc_pct"],
            "Mapped Reads (%)": mapped_pct,
            "Mean Depth (x)": mean_depth,
            "Total Filtered Variants": vcf_stats["total_filtered_variants"],
            "SNPs": vcf_stats["snps"],
            "Indels": vcf_stats["indels"],
            "Ts/Tv Ratio": vcf_stats["ts_tv"],
            "High Impact": vcf_stats["high_impact"],
            "Moderate Impact": vcf_stats["moderate_impact"]
        })
        
    df_results = pd.DataFrame(records)
    
    # Save output to ready_to_download
    output_dir = "results/ready_to_download"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "wgs_pipeline_statistics.tsv")
    df_results.to_csv(output_file, sep="\t", index=False)
    
    print(f"\nSuccess! All statistics gathered and saved to: {output_file}")
    print("You can download this file to copy the exact figures into your report.")

if __name__ == "__main__":
    main()
