const https = require('https');

const accessions = [
    "ERR9538944", "ERR9538945", "ERR9538947", "ERR9538952", "ERR9538958",
    "SRR21679075", "SRR31897979", "SRR31897980", "SRR31897981", "SRR31897982",
    "SRR26412217", "SRR26412242", "SRR26412243", "SRR26412254", "SRR26412293",
    "SRR6208298", "SRR6208299", "SRR6208300", "SRR6208301", "SRR6208302"
];

console.log("Checking layout for SRA accessions...");
console.log("-".repeat(50));
console.log(`${"Accession".padEnd(15)} | ${"Layout".padEnd(10)} | Status`);
console.log("-".repeat(50));

function checkAccession(index) {
    if (index >= accessions.length) {
        console.log("-".repeat(50));
        console.log("Done.");
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
