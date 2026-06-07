import urllib.request
import xml.etree.ElementTree as ET
import time
import sys

# List of accessions to check
accessions = [
    "SRR31897984", "SRR31897983", "SRR31897982", "SRR31897981", "SRR21679075",
    "SRR7964123", "SRR7964124", "SRR7964125", "SRR7964126", "SRR7964127",
    "DRR076141", "DRR076142", "DRR076143", "DRR076144", "DRR076145"
]

print("Checking layout for SRA accessions...", flush=True)
print("-" * 50)
print(f"{'Accession':<15} | {'Layout':<10} | {'Status'}")
print("-" * 50)

for acc in accessions:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={acc}&retmode=xml"
    try:
        req = urllib.request.Request(
            url, 
            headers={'User-Agent': 'Mozilla/5.0'}
        )
        with urllib.request.urlopen(req, timeout=10) as response:
            content = response.read().decode('utf-8')
            
        if "<PAIRED" in content or "paired" in content.lower():
            layout = "PAIRED"
        elif "<SINGLE" in content or "single" in content.lower():
            layout = "SINGLE"
        else:
            layout = "UNKNOWN"
            
        print(f"{acc:<15} | {layout:<10} | OK")
    except Exception as e:
        print(f"{acc:<15} | {'ERROR':<10} | {str(e)}")
    
    # Avoid rate-limiting by NCBI
    time.sleep(0.5)

print("-" * 50)
print("Selesai.")
