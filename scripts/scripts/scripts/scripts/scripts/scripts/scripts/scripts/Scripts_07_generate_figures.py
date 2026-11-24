python3 -c "
import csv
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

variants = ['APOE\nrs429358 (ε4)', 'APOE\nrs7412 (ε2)', 'ABCA7\nrs115550680', 'BIN1\nrs744373', 'TREM2\nrs75932628']
yri_vals = [23.1, 10.2, 14.8, 55.6, 0.0]
asw_vals = [20.5, 13.9,  7.4, 51.6, 0.0]
ceu_vals = [17.7,  6.6,  0.0, 29.8, 1.0]

x = np.arange(len(variants))
w = 0.26

fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=300)
b1 = ax.bar(x - w, yri_vals, w, label='Ibadan, Nigeria (YRI, n=108)', color='#1b4f72')
b2 = ax.bar(x,     asw_vals, w, label='African American (ASW, n=61)', color='#2e86c1')
b3 = ax.bar(x + w, ceu_vals, w, label='European Ancestry (CEU, n=99)', color='#aab7b8')

for bars in [b1, b2, b3]:
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 1.0, f'{h:.1f}%', ha='center', va='bottom', fontsize=7.5)

ax.set_ylabel('Allele Frequency (%)', fontsize=10, fontweight='bold')
ax.set_title('Comparative Allele Frequencies of Key Alzheimer\'s Susceptibility Variants', fontsize=11, fontweight='bold', pad=12)
ax.set_xticks(x)
ax.set_xticklabels(variants, fontsize=9)
ax.set_ylim(0, 65)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', linestyle='--', alpha=0.4)
ax.legend(frameon=False, fontsize=8.5, loc='upper left')
plt.tight_layout()
plt.savefig('Figure1_Allele_Frequencies_YRI_ASW_CEU.png', dpi=300)
plt.close()

with open('ibadan_alzheimers_panel.csv') as f:
    rows = list(csv.DictReader(f))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4), dpi=300)

genes = ['APOE\n(3.6 kb)', 'TREM2\n(4.7 kb)', 'BIN1\n(59.3 kb)', 'ABCA7\n(25.5 kb)']
pi_vals = [0.48, 0.68, 1.55, 3.20]
densities = ['4.7/kb', '6.6/kb', '11.8/kb', '28.0/kb']
colors_a = ['#5d6d7e', '#5d6d7e', '#2874a6', '#117864']

bars_a = ax1.bar(genes, pi_vals, color=colors_a, width=0.55)
for bar, d in zip(bars_a, densities):
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, h + 0.08, f'{h:.2f}\n({d})', ha='center', va='bottom', fontsize=7.8)

ax1.set_ylabel('Nucleotide Diversity (π × 10⁻³)', fontsize=9.5, fontweight='bold')
ax1.set_title('A. Gene-Wide Sequence Diversity in Ibadan (YRI)', fontsize=10, fontweight='bold', pad=10)
ax1.set_ylim(0, 3.9)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.grid(axis='y', linestyle='--', alpha=0.4)

diplotypes = ['e2/e3', 'e2/e4', 'e3/e3', 'e3/e4', 'e4/e4']
abca7_neg = [sum(1 for r in rows if r['APOE_Diplotype']==d and r['ABCA7_rs115550680_GT']=='0/0') for d in diplotypes]
abca7_pos = [sum(1 for r in rows if r['APOE_Diplotype']==d and r['ABCA7_rs115550680_GT']!='0/0') for d in diplotypes]

ax2.bar(diplotypes, abca7_neg, label='ABCA7 rs115550680 Non-Carrier (0/0)', color='#85c1e9', width=0.55)
ax2.bar(diplotypes, abca7_pos, bottom=abca7_neg, label='ABCA7 rs115550680 Carrier (0/1 or 1/1)', color='#1b4f72', width=0.55)

for i, (neg, pos) in enumerate(zip(abca7_neg, abca7_pos)):
    tot = neg + pos
    ax2.text(i, tot + 0.8, f'n={tot}', ha='center', va='bottom', fontsize=8.5, fontweight='bold')

ax2.set_ylabel('Number of Ibadan Individuals (N = 108)', fontsize=9.5, fontweight='bold')
ax2.set_xlabel('APOE Diplotype', fontsize=9.5, fontweight='bold')
ax2.set_title('B. APOE Diplotypes & ABCA7 Co-Occurrence', fontsize=10, fontweight='bold', pad=10)
ax2.set_ylim(0, 55)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.grid(axis='y', linestyle='--', alpha=0.4)
ax2.legend(frameon=False, fontsize=7.8, loc='upper left')

plt.tight_layout()
plt.savefig('Figure2_Ibadan_Diversity_and_CoOccurrence.png', dpi=300)
plt.close()
"
