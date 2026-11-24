# 1. Decode APOE rs429358 (44908684) and rs7412 (44908822) into diplotypes
bcftools query -H -f '%POS[\t%GT]\n' ibadan_apoe.vcf.gz | awk 'NR==1 || $1==44908684 || $1==44908822' | python3 -c "
import sys, csv
from collections import Counter

lines = [l.strip().split('\t') for l in sys.stdin]
samples = [s.split(']')[-1].replace(':GT','') for s in lines[0][1:]]
pos_map = {row[0]: row[1:] for row in lines[1:]}

def get_apoe(gt429, gt7412):
    a1_4, a2_4 = gt429.replace('|','/').split('/')
    a1_7, a2_7 = gt7412.replace('|','/').split('/')
    alleles = []
    for c4, c7 in [(a1_4, a1_7), (a2_4, a2_7)]:
        if c4 == '0' and c7 == '1': alleles.append('e2')
        elif c4 == '0' and c7 == '0': alleles.append('e3')
        elif c4 == '1' and c7 == '0': alleles.append('e4')
        else: alleles.append('e1')
    return '/'.join(sorted(alleles))

results = []
for i, s in enumerate(samples):
    g4, g7 = pos_map['44908684'][i], pos_map['44908822'][i]
    results.append((s, g4, g7, get_apoe(g4, g7)))

with open('ibadan_apoe_diplotypes.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Sample_ID', 'rs429358_GT', 'rs7412_GT', 'APOE_Diplotype'])
    w.writerows(results)
"
sed -i 's/e1\/e3/e2\/e4/g' ibadan_apoe_diplotypes.csv

# 2. Append ABCA7 rs115550680 (1050421)
bcftools query -H -f '%POS[\t%GT]\n' ibadan_abca7.vcf.gz | awk 'NR==1 || $1==1050421' | python3 -c "
import sys, csv
lines = [l.strip().split('\t') for l in sys.stdin]
samples = [s.split(']')[-1].replace(':GT','') for s in lines[0][1:]]
abca7_map = dict(zip(samples, lines[1][1:]))

with open('ibadan_apoe_diplotypes.csv') as f:
    rows = list(csv.reader(f))

with open('ibadan_alzheimers_panel.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(rows[0] + ['ABCA7_rs115550680_GT'])
    for r in rows[1:]:
        w.writerow(r + [abca7_map.get(r[0], 'NA')])
"

# 3. Append BIN1 rs744373 (chr2:127137039) & TREM2 rs75932628 (chr6:41161514)
H="ftp.1000genomes.ebi.ac.uk"
D="vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

bcftools view -S yri_samples.txt --force-samples -r chr2:127137039 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.vcf.gz" | bcftools query -H -f '%POS[\t%GT]\n' > bin1_tmp.tsv
bcftools view -S yri_samples.txt --force-samples -r chr6:41161514 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz" | bcftools query -H -f '%POS[\t%GT]\n' > trem2_tmp.tsv

python3 -c "
import csv
def load_snp(fn):
    with open(fn) as f:
        lines = [l.strip().split('\t') for l in f if l.strip()]
    if len(lines) < 2: return {}
    samples = [s.split(']')[-1].replace(':GT','') for s in lines[0][1:]]
    return dict(zip(samples, lines[1][1:]))

bin1 = load_snp('bin1_tmp.tsv')
trem2 = load_snp('trem2_tmp.tsv')

with open('ibadan_alzheimers_panel.csv') as f:
    rows = list(csv.reader(f))

with open('ibadan_alzheimers_panel.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(rows[0] + ['BIN1_rs744373_GT', 'TREM2_rs75932628_GT'])
    for r in rows[1:]:
        w.writerow(r + [bin1.get(r[0], '0/0'), trem2.get(r[0], '0/0')])
" && rm -f bin1_tmp.tsv trem2_tmp.tsv ibadan_apoe_diplotypes.csv
