"""Schematic ward map: linkage of the widespread shared species between air (per ward) and ward E surface rooms."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]          # manuscript/
DATA = str(ROOT / 'data') + '/'
OUTD = str(ROOT / 'notes' / 'shared_fraction') + '/'
Path(OUTD).mkdir(parents=True, exist_ok=True)

import pandas as pd, numpy as np, matplotlib, itertools, re
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
P = DATA
OUT = OUTD
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8.5})
ws = pd.read_csv(OUT + 'shared_fraction_widespread_species.csv', index_col=0).index.tolist()
m = pd.read_csv(P + 'metadata/hospital_metadata.csv')
d = pd.read_csv(P + 'processed/taxa_counts_floor10.csv.gz')
keep = set(pd.read_csv(P + 'processed/alpha_diversity.csv')['sample'])
d = d[(d.tax_rank == 'S') & d['sample'].isin(keep) & d['name'].isin(ws)]
room = {}
for _, r in m[m.type == 'Surface'].iterrows():
    x = re.split(r'[,;]', r.description)[0].strip()
    x = {'Chemo Freezer': 'Nurse bay', 'Helmer Blood bank fridge': 'Blood bank', 'Stem cell fridge': 'Stem cell unit',
         'Staff room freezer': 'Staff room', 'Corridor Vent': 'Corridor'}.get(x, x)
    room[r['sample']] = 'E: ' + x
site = {s: (room[s] if s.startswith('S') else 'Air ' + m.set_index('sample').ward[s]) for s in d['sample'].unique()}
d['site'] = d['sample'].map(site)
pa = d.groupby(['name', 'site']).size().unstack(fill_value=0).gt(0).astype(int)
nsamp = pd.Series(site).value_counts()
mm = m.set_index('sample')
site.update({x: 'Air ' + mm.ward[x] for x in keep if x.startswith('A')}); site.update({x: room[x] for x in keep if x.startswith('S')})
for w in sorted({v for v in site.values() if v.startswith('Air')}):
    if w not in pa.columns: pa[w] = 0
nsamp = pd.Series(site).value_counts()
air = sorted([c for c in pa.columns if c.startswith('Air')]); surf = sorted([c for c in pa.columns if c.startswith('E:')])

fig = plt.figure(figsize=(12, 6.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.25)
ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off()
# layout: air wards on a left arc, ward E (box) on the right with rooms on a grid
pos = {}
for i, w in enumerate(air):
    pos[w] = (0.20, 0.93 - i * 0.088)
cols = 3
for i, r_ in enumerate(surf):
    pos[r_] = (0.58 + (i % cols) * 0.16, 0.80 - (i // cols) * 0.19)
ax.add_patch(FancyBboxPatch((0.50, 0.05), 0.48, 0.90, boxstyle='round,pad=0.01', fc='#fdf0f6', ec='#ed8cb6', lw=1.2))
ax.text(0.74, 0.97, 'Ward E (haematology/oncology): surfaces', ha='center', fontweight='bold', color='#b8457a')
ax.text(0.13, 0.99, 'Air, by ward', ha='center', fontweight='bold', color='#5a5b5d')
# edges: air ward -- surface room, weight = number of widespread species shared
for a in air:
    for s in surf:
        k = int((pa[a] & pa[s]).sum())
        if k >= 5:
            (x1, y1), (x2, y2) = pos[a], pos[s]
            ax.plot([x1, x2], [y1, y2], c='#CB181D', lw=0.3 + 0.5 * (k - 4), alpha=0.15 + 0.06 * (k - 4), zorder=1)
for n, (x, y) in pos.items():
    k = int(pa[n].sum()); c = '#828385' if n.startswith('Air') else '#ed8cb6'
    ax.scatter(x, y, s=25 + 30 * k, c=c, edgecolors='k', lw=.6, zorder=3)
    if n.startswith('Air'):
        ax.text(x - 0.05, y, f"{n}  ({k} spp; n = {nsamp[n]})", ha='right', va='center', fontsize=7, zorder=4)
    else:
        ax.text(x, y - 0.05, f"{n.replace('E: ', '')}\n({k} spp; n = {nsamp[n]})", ha='center', va='top', fontsize=7, zorder=4)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.set_title('A  Air wards linked to ward E rooms by ≥5 widespread shared species', loc='left', fontweight='bold')
for k in (5, 8, 11):
    ax.plot([], [], c='#CB181D', lw=0.3 + 0.5 * (k - 4), alpha=0.25 + 0.06 * (k - 4), label=f'{k} species shared')
ax.legend(loc='lower left', frameon=False, fontsize=7, bbox_to_anchor=(0.0, -0.02))

# B: presence matrix species x site
ax = fig.add_subplot(gs[0, 1])
order = pa.loc[ws].sum(axis=1).sort_values(ascending=False).index
mat = pa.loc[order, air + surf]
ax.imshow(np.where(mat.values, 1, 0), cmap=matplotlib.colors.ListedColormap(['#f2f2f2', '#CB181D']), aspect='auto')
ax.set_yticks(range(len(order))); ax.set_yticklabels([n.split()[0][0] + '. ' + n.split()[1] for n in order], style='italic')
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels([c.replace('E: ', '') for c in mat.columns], rotation=90, fontsize=7.5)
ax.axvline(len(air) - .5, c='k', lw=1)
ax.text((len(air) - 1) / 2, -0.75, 'Air (ward)', ha='center', fontweight='bold', color='#5a5b5d')
ax.text(len(air) + (len(surf) - 1) / 2, -0.75, 'Ward E surfaces (room)', ha='center', fontweight='bold', color='#b8457a')
ax.set_title('B  Where each widespread shared species was detected', loc='left', fontweight='bold', pad=24)
fig.savefig(OUT + 'shared_fraction_ward_map.png', dpi=300, bbox_inches='tight')
pa.loc[order, air + surf].to_csv(OUT + 'shared_fraction_site_presence.csv')
print('air sites', air); print('surface rooms', surf)
print('edges', sum(int((pa[a] & pa[s]).sum()) > 0 for a in air for s in surf), 'of', len(air) * len(surf))
print((pa.loc[order, air].sum(axis=1)).to_dict()); print((pa.loc[order, surf].sum(axis=1)).to_dict())
