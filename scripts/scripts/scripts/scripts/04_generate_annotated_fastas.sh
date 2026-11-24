U="api.genome.ucsc.edu"

# 1. Fetch GRCh38 reference FASTA slices for all 4 genes
curl -s "https://${U}/getData/sequence?genome=hg38;chrom=chr19;start=44905795;end=44909393" | python3 -c "import sys,json; seq=json.load(sys.stdin)['dna'].upper(); open('apoe_ref_hg38.fa','w').write('>chr19:44905796-44909393\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')"
curl -s "https://${U}/getData/sequence?genome=hg38;chrom=chr19;start=1040100;end=1065571" | python3 -c "import sys,json; seq=json.load(sys.stdin)['dna'].upper(); open('abca7_ref_hg38.fa','w').write('>chr19:1040101-1065571\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')"
curl -s "https://${U}/getData/sequence?genome=hg38;chrom=chr2;start=127048022;end=127107288" | python3 -c "import sys,json; seq=json.load(sys.stdin)['dna'].upper(); open('bin1_ref_hg38.fa','w').write('>chr2:127048023-127107288\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')"
curl -s "https://${U}/getData/sequence?genome=hg38;chrom=chr6;start=41158505;end=41163186" | python3 -c "import sys,json; seq=json.load(sys.stdin)['dna'].upper(); open('trem2_ref_hg38.fa','w').write('>chr6:41158506-41163186\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')"

# 2. Build annotated diploid multi-FASTA files (216 haplotypes per gene)
for GENE in APOE ABCA7 BIN1 TREM2; do
    LOWER=$(echo "$GENE" | tr '[:upper:]' '[:lower:]')
    rm -f "all_ibadan_${LOWER}_tmp.fasta"
    while read -r SAMPLE; do
        bcftools consensus -f "${LOWER}_ref_hg38.fa" -s "$SAMPLE" -H 1 "ibadan_${LOWER}.vcf.gz" 2>/dev/null | sed "s/>.*/>${SAMPLE}_${GENE}_H1/" >> "all_ibadan_${LOWER}_tmp.fasta"
        bcftools consensus -f "${LOWER}_ref_hg38.fa" -s "$SAMPLE" -H 2 "ibadan_${LOWER}.vcf.gz" 2>/dev/null | sed "s/>.*/>${SAMPLE}_${GENE}_H2/" >> "all_ibadan_${LOWER}_tmp.fasta"
    done < yri_samples.txt

    python3 -c "
import csv
with open('ibadan_alzheimers_panel.csv') as f:
    meta = {r['Sample_ID']: r for r in csv.DictReader(f)}
with open('all_ibadan_${LOWER}_tmp.fasta') as fin, open('all_ibadan_${LOWER}_annotated.fasta', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            sid = line[1:].split('_')[0]
            m = meta.get(sid, {})
            fout.write(f'{line.strip()}|APOE={m.get(\"APOE_Diplotype\",\"NA\")}|ABCA7={m.get(\"ABCA7_rs115550680_GT\",\"NA\")}|BIN1={m.get(\"BIN1_rs744373_GT\",\"NA\")}\n')
        else:
            fout.write(line)
" && rm -f "all_ibadan_${LOWER}_tmp.fasta" "${LOWER}_ref_hg38.fa"
done
