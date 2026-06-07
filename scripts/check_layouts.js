const https = require('https');

const accessions = [
    "SRR31897984", "SRR31897983", "SRR31897982", "SRR31897981", "SRR21679075",
    "SRR7964123", "SRR7964124", "SRR7964125", "SRR7964126", "SRR7964127",
    "DRR076141", "DRR076142", "DRR076143", "DRR076144", "DRR076145"
];

console.log("Checking layout for SRA accessions...");
console.log("-".repeat(50));
console.log(`${"Accession".padEnd(15)} | ${"Layout".padEnd(10)} | Status`);
console.log("-".repeat(50));

function checkAccession(index) {
    if (index >= accessions.length) {
        console.log("-".repeat(50));
        console.log("Selesai.");
        return;
    }

    const acc = accessions[index];
    const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${acc}&retmode=xml`;

    const req = https.get(url, { headers: { 'User-Agent': 'Mozilla/5.0' } }, (res) => {
        let data = '';
        res.on('data', (chunk) => { data += chunk; });
        res.on('end', () => {
            let layout = 'UNKNOWN';
            if (data.includes('<PAIRED') || data.toLowerCase().includes('paired')) {
                layout = 'PAIRED';
            } else if (data.includes('<SINGLE') || data.toLowerCase().includes('single')) {
                layout = 'SINGLE';
            }
            console.log(`${acc.padEnd(15)} | ${layout.padEnd(10)} | OK`);
            setTimeout(() => checkAccession(index + 1), 500);
        });
    });

    req.on('error', (err) => {
        console.log(`${acc.padEnd(15)} | ${"ERROR".padEnd(10)} | ${err.message}`);
        setTimeout(() => checkAccession(index + 1), 500);
    });
}

checkAccession(0);
