H="ftp.1000genomes.ebi.ac.uk"
P="vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
D="vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

curl -s "http://${H}/${P}" > all_panels.tmp
awk '$2=="CEU" {print $1}' all_panels.tmp > ceu.tmp
awk '$2=="ASW" {print $1}' all_panels.tmp > asw.tmp

for POP in CEU ASW; do
    FILE=$(echo "$POP" | tr '[:upper:]' '[:lower:]').tmp
    bcftools view -S "$FILE" --force-samples -r chr19:1050421,chr19:44908684,chr19:44908822 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.vcf.gz" | bcftools +fill-tags -- -t AF | bcftools query -f "${POP}\t%POS\t%INFO/AF\n" >> pop_comp.tmp
    bcftools view -S "$FILE" --force-samples -r chr2:127137039 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.vcf.gz" | bcftools +fill-tags -- -t AF | bcftools query -f "${POP}\t%POS\t%INFO/AF\n" >> pop_comp.tmp
    bcftools view -S "$FILE" --force-samples -r chr6:41161514 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz" | bcftools +fill-tags -- -t AF | bcftools query -f "${POP}\t%POS\t%INFO/AF\n" >> pop_comp.tmp
done

python3 -c "
import csv
with open('ibadan_alzheimers_panel.csv') as f:
    yri = list(csv.DictReader(f))

def yri_af(col):
    return sum(r[col].replace('|','/').split('/').count('1') for r in yri) / (2 * len(yri)) * 100

other = {}
with open('pop_comp.tmp') as f:
    for line in f:
        pop, pos, af = line.strip().split('\t')
        other[(pop, pos)] = float(af) * 100

snps = [
    ('APOE rs429358 (e4 risk)',  '44908684',  'rs429358_GT'),
    ('APOE rs7412 (e2 prot.)',   '44908822',  'rs7412_GT'),
    ('ABCA7 rs115550680 (Afr)',  '1050421',   'ABCA7_rs115550680_GT'),
    ('BIN1 rs744373 (GWAS)',     '127137039', 'BIN1_rs744373_GT'),
    ('TREM2 rs75932628 (R47H)',  '41161514',  'TREM2_rs75932628_GT'),
]

print(f'{\"Locus / Variant\":26s} | {\"Ibadan (YRI)\":12s} | {\"Afr-Amer (ASW)\":14s} | {\"European (CEU)\":14s}')
print('-'*74)
for label, pos, col in snps:
    print(f'{label:26s} | {yri_af(col):10.1f}%  | {other.get((\"ASW\",pos),0):12.1f}%  | {other.get((\"CEU\",pos),0):12.1f}%')
" && rm -f all_panels.tmp ceu.tmp asw.tmp pop_comp.tmp
