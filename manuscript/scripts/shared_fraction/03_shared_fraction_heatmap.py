"""Heatmap: top 30 shared species (by prevalence in both compartments) per sample, grouped by ward / ward E room."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]          # manuscript/
DATA = str(ROOT / 'data') + '/'
OUTD = str(ROOT / 'notes' / 'shared_fraction') + '/'
Path(OUTD).mkdir(parents=True, exist_ok=True)

import pandas as pd, numpy as np, matplotlib, re
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = DATA
D = OUTD
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8})
sh = pd.read_csv(D + 'shared_fraction_all_species.csv', index_col=0)
sh['minprev'] = np.minimum(sh.pA, sh.pS)
top = sh.sort_values(['minprev', 'aA'], ascending=False).head(30)
m = pd.read_csv(P + 'metadata/hospital_metadata.csv').set_index('sample')
d = pd.read_csv(P + 'processed/taxa_counts_floor10.csv.gz')
keep = set(pd.read_csv(P + 'processed/alpha_diversity.csv')['sample'])
d = d[(d.tax_rank == 'S') & d['sample'].isin(keep)].copy()
d['rel'] = d['count'] / d.groupby('sample')['count'].transform('sum') * 100
mat = d[d.name.isin(top.index)].pivot_table(index='name', columns='sample', values='rel', aggfunc='sum').reindex(index=top.index, columns=sorted(keep)).fillna(0)
def room(s):
    x = re.split(r'[,;]', m.description[s])[0].strip()
    return {'Chemo Freezer': 'Nurse bay', 'Helmer Blood bank fridge': 'Blood bank', 'Stem cell fridge': 'Stem cell', 'Staff room freezer': 'Staff room', 'Corridor Vent': 'Corridor'}.get(x, x)
grp = {s: ('Air', m.ward[s]) if s.startswith('A') else ('Surface', room(s)) for s in mat.columns}
order = sorted(mat.columns, key=lambda s: (grp[s][0], grp[s][1], int(re.sub(r'\D', '', s))))
mat = mat[order]
tier = top.tier.str.split(' ').str[0]
fig, ax = plt.subplots(figsize=(15, 7.2))
v = np.log10(mat.values + 1e-3); v[mat.values == 0] = np.nan
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('paper', ['#FFF5F0', '#FCBBA1', '#FB6A4A', '#CB181D', '#67000D']); cmap.set_bad('#f2f2f2')
im = ax.imshow(v, aspect='auto', cmap=cmap, vmin=-3, vmax=2, interpolation='nearest')
ax.set_yticks(range(len(top))); ax.set_yticklabels([n.split()[0][0] + '. ' + ' '.join(n.split()[1:]) for n in top.index], style='italic')
for t, lab in zip(ax.get_yticklabels(), tier): t.set_color('#CB181D' if lab == 'Widespread' else 'k')
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6.5)
# group separators and labels
prev = None; start = 0
groups = [grp[s] for s in order] + [None]
for i, g in enumerate(groups):
    if g != prev and prev is not None:
        ax.axvline(i - .5, c='white', lw=1.5 if g is None or g[0] == prev[0] else 0)
        if g is not None and g[0] != prev[0]: ax.axvline(i - .5, c='k', lw=1.5)
        lab = ('Ward ' + prev[1]) if prev[0] == 'Air' else prev[1]
        ax.text((start + i - 1) / 2, -0.9, lab, rotation=60, ha='left', va='bottom', fontsize=7)
        start = i
    prev = g
n_air = sum(1 for s in order if s.startswith('A'))
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label('log$_{10}$ relative abundance (%)')
ax.set_title('Top 30 shared species by prevalence in both compartments (red: widespread tier).  Left: air samples by ward  |  right: ward E surface samples by room', loc='left', fontweight='bold', pad=62)
fig.savefig(D + 'shared_fraction_heatmap_top30.png', dpi=300, bbox_inches='tight')
top[['nA', 'nS', 'aA', 'aS', 'tier']].to_csv(D + 'shared_fraction_top30.csv')
print(tier.value_counts().to_dict()); print(top[['nA', 'nS', 'phylum']].to_string())
