import urllib.request
import xml.etree.ElementTree as ET
import time
import sys

# List of accessions to check
accessions = [
    "ERR9538944", "ERR9538945", "ERR9538947", "ERR9538952", "ERR9538958",
    "SRR21679075", "SRR31897979", "SRR31897980", "SRR31897981", "SRR31897982",
    "SRR26412217", "SRR26412242", "SRR26412243", "SRR26412254", "SRR26412293",
    "SRR6208298", "SRR6208299", "SRR6208300", "SRR6208301", "SRR6208302"
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
print("Done.")
