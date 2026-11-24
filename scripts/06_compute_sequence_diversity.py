python3 -c "
import subprocess

genes = [
    ('APOE',  'ibadan_apoe.vcf.gz',  3598),
    ('ABCA7', 'ibadan_abca7.vcf.gz', 25471),
    ('BIN1',  'ibadan_bin1.vcf.gz',  59266),
    ('TREM2', 'ibadan_trem2.vcf.gz', 4681),
]

print(f'{\"Gene\":8s} | {\"Length (bp)\":11s} | {\"Polymorphic Sites\":17s} | {\"Density (/kb)\":13s} | {\"Nucleotide Div (pi)\":19s}')
print('-'*80)

for name, vcf, length in genes:
    cmd = f'bcftools norm -m -any {vcf} 2>/dev/null | bcftools view -c 1 | bcftools +fill-tags -- -t AF | bcftools query -f \"%INFO/AF\\n\"'
    out = subprocess.check_output(cmd, shell=True).decode().strip().splitlines()
    poly = len(out)
    n_chr = 216
    pi_sum = sum(2 * float(af) * (1 - float(af)) * (n_chr / (n_chr - 1)) for af in out if af != '.')
    pi = pi_sum / length
    print(f'{name:8s} | {length:11,d} | {poly:17,d} | {poly/(length/1000):13.2f} | {pi:19.5f}')
"
