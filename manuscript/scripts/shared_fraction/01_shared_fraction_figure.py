"""Shared air-surface species: prevalence-abundance figure + widespread-species table (for collaborator note)."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]          # manuscript/
DATA = str(ROOT / 'data') + '/'
OUTD = str(ROOT / 'notes' / 'shared_fraction') + '/'
Path(OUTD).mkdir(parents=True, exist_ok=True)

import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.lines import Line2D
P = DATA + 'processed/'
OUT = OUTD
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9})

d = pd.read_csv(P + 'taxa_counts_floor10.csv.gz')
keep = set(pd.read_csv(P + 'alpha_diversity.csv')['sample'])
d = d[(d.tax_rank == 'S') & d['sample'].isin(keep)].copy()
d['grp'] = np.where(d['sample'].str.startswith('A'), 'Air', 'Surface')
d['rel'] = d['count'] / d.groupby('sample')['count'].transform('sum') * 100
ns = d.groupby('grp')['sample'].nunique()
pres = d.groupby(['name', 'grp'])['sample'].nunique().unstack(fill_value=0)
ab = d.groupby(['name', 'grp'])['rel'].sum().unstack(fill_value=0) / ns
phy = d.drop_duplicates('name').set_index('name')['phylum']
t = pd.DataFrame({'nA': pres.Air, 'nS': pres.Surface, 'pA': pres.Air / ns['Air'] * 100, 'pS': pres.Surface / ns['Surface'] * 100,
                  'aA': ab.Air, 'aS': ab.Surface, 'phylum': phy})
t['cat'] = np.select([(t.nA > 0) & (t.nS == 0), (t.nA == 0) & (t.nS > 0)], ['Air only', 'Surface only'], 'Shared')
sh = t[t.cat == 'Shared'].copy()
sh['tier'] = np.select([(sh.nA >= 5) & (sh.nS >= 5), (sh.nA >= 2) & (sh.nS >= 2), (sh.nA == 1) | (sh.nS == 1)],
                       ['Widespread (\u22655 samples in both)', 'Recurrent (2\u20134 in both)', 'Sporadic (1 sample in one compartment)'], 'other')
assert (sh.tier != 'other').all()
core = sh[sh.tier.str.startswith('Widespread')].copy()
print(t.cat.value_counts().to_dict(), sh.tier.value_counts().to_dict())

# genetic elements carried by core species
e = pd.read_csv(P + 'abri_kraken2_cleaned.csv', usecols=['sample', 'GENE', 'DATABASE', 'name'])
e = e[e['sample'].isin(keep) & e['name'].isin(core.index)]
db = {'card': 'ARGs', 'vfdb': 'VFs', 'plasmidfinder': 'MGEs'}
e['cls'] = e.DATABASE.map(db)
g = e.groupby(['name', 'cls'])['GENE'].nunique().unstack(fill_value=0).reindex(columns=['ARGs', 'VFs', 'MGEs'], fill_value=0)
core = core.join(g).fillna({'ARGs': 0, 'VFs': 0, 'MGEs': 0})
ex = e.groupby('name')['GENE'].apply(lambda s: ', '.join(sorted(s.unique())[:6]) + (' \u2026' if s.nunique() > 6 else ''))
core['examples'] = ex
core['eA'] = e[e['sample'].str.startswith('A')].groupby('name')['sample'].nunique()
core['eS'] = e[e['sample'].str.startswith('S')].groupby('name')['sample'].nunique()
core = core.fillna({'eA': 0, 'eS': 0}).sort_values('pS', ascending=False)
core.to_csv(OUT + 'shared_fraction_widespread_species.csv')
sh.to_csv(OUT + 'shared_fraction_all_species.csv')

# ---------- figure ----------
C = {'Air only': '#828385', 'Surface only': '#ed8cb6', 'Shared': '#b88a9e'}
TC = {'Widespread (\u22655 samples in both)': '#CB181D', 'Recurrent (2\u20134 in both)': '#FCBBA1', 'Sporadic (1 sample in one compartment)': '#D9D9D9'}
fig = plt.figure(figsize=(11, 8.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], hspace=0.38, wspace=0.28)

# A: prevalence-abundance, all species by Venn category
ax = fig.add_subplot(gs[0, 0])
t['prev'] = (t.nA + t.nS) / ns.sum() * 100
t['abun'] = (t.aA * ns['Air'] + t.aS * ns['Surface']) / ns.sum()
for c in ['Air only', 'Surface only', 'Shared']:
    s = t[t.cat == c]; ax.scatter(s.prev, np.log10(s.abun), s=10, c=C[c], alpha=.6, lw=0, label=f'{c} (n = {len(s)})')
ax.scatter(t.loc[core.index, 'prev'], np.log10(t.loc[core.index, 'abun']), s=28, facecolors='none', edgecolors='#CB181D', lw=1.2, label='Widespread shared species (n = 13)')
tx = [ax.text(t.prev[n], np.log10(t.abun[n]), n.split()[0][0] + '. ' + n.split()[1], fontsize=7.5, style='italic')
      for n in ['Rahnella bonaserana', 'Acinetobacter lwoffii', 'Stutzerimonas stutzeri', 'Pseudomonas aeruginosa', 'Paracoccus marcusii', 'Escherichia coli']]
adjust_text(tx, ax=ax, arrowprops=dict(arrowstyle='-', color='grey', lw=.5))
ax.set_xlabel('Prevalence, all samples (%)'); ax.set_ylabel('log$_{10}$ mean relative abundance (%)')
ax.set_title('A  Prevalence and abundance of all species', loc='left', fontweight='bold'); ax.legend(fontsize=7, frameon=False, loc='lower right')

# B: shared species, air vs surface prevalence
ax = fig.add_subplot(gs[0, 1])
for k, col in TC.items():
    s = sh[sh.tier == k]; ax.scatter(s.pA, s.pS, s=8 + 6 * np.sqrt(s.aA + s.aS) * 4, c=col, alpha=.75, edgecolors='grey', lw=.3, label=f'{k.split(" (")[0]} (n = {len(s)})')
tx = [ax.text(r.pA, r.pS, n.split()[0][0] + '. ' + n.split()[1], fontsize=7, style='italic') for n, r in core.iterrows()]
adjust_text(tx, ax=ax, arrowprops=dict(arrowstyle='-', color='grey', lw=.5), expand=(1.3, 1.6))
ax.axvline(5 / ns['Air'] * 100, ls=':', c='grey', lw=.8); ax.axhline(5 / ns['Surface'] * 100, ls=':', c='grey', lw=.8)
ax.set_xlabel(f'Prevalence in air (% of {ns["Air"]} samples)'); ax.set_ylabel(f'Prevalence on surfaces (% of {ns["Surface"]} samples)')
ax.set_title('B  The shared fraction (590 species)', loc='left', fontweight='bold'); ax.legend(fontsize=7, frameon=False, loc='center left', bbox_to_anchor=(0.0, 0.55)); ax.set_ylim(0, 75)

# C: breakdown of shared species: tier x phylum
ax = fig.add_subplot(gs[1, 0])
top = sh.phylum.value_counts().head(5).index
sh['phy'] = np.where(sh.phylum.isin(top), sh.phylum, 'Other')
order = list(TC)[::-1]
ct = sh.groupby(['tier', 'phy']).size().unstack(fill_value=0).reindex(order)
cols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#DDCC77', '#CC6677', '#882255'][:ct.shape[1]]
left = np.zeros(len(ct))
for i, p in enumerate(ct.columns):
    ax.barh([o.split(' (')[0] for o in order], ct[p], left=left, color=cols[i], label=p); left += ct[p].values
for y, v in enumerate(left): ax.text(v + 5, y, str(int(v)), va='center', fontsize=8)
ax.set_xlabel('Number of shared species'); ax.legend(fontsize=7, frameon=False, ncol=2, loc='upper right')
ax.set_title('C  Shared fraction by prevalence tier and phylum', loc='left', fontweight='bold')

# D: core species: prevalence + elements
ax = fig.add_subplot(gs[1, 1])
cs = core.sort_values('pS')
y = np.arange(len(cs))
ax.barh(y + .2, cs.pA, .4, color=C['Air only'], label='Air'); ax.barh(y - .2, cs.pS, .4, color=C['Surface only'], label='Surface')
ax.set_yticks(y); ax.set_yticklabels([n.split()[0][0] + '. ' + n.split()[1] for n in cs.index], style='italic', fontsize=8)
for i, (n, r) in enumerate(cs.iterrows()):
    ax.text(102, i, f'{int(r.ARGs)}/{int(r.VFs)}/{int(r.MGEs)}', va='center', fontsize=7.5)
ax.text(102, len(cs) - .2, 'ARG/VF/MGE', fontsize=7, fontweight='bold')
ax.set_xlim(0, 125); ax.set_xlabel('Prevalence (% of samples)'); ax.legend(fontsize=7, frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.17), ncol=2)
ax.set_title('D  Widespread shared species: prevalence and elements', loc='left', fontweight='bold')
for a in fig.axes: a.spines[['top', 'right']].set_visible(False)
fig.savefig(OUT + 'shared_fraction_prevalence_abundance.png', dpi=300, bbox_inches='tight')
print(core[['nA', 'nS', 'aA', 'aS', 'ARGs', 'VFs', 'MGEs', 'eA', 'eS']].round(2).to_string())
print(core['examples'].to_string())
print('shared tiers x phylum'); print(ct.to_string())
print('abundance share of widespread species:', (d[d.name.isin(core.index)].groupby('grp')['rel'].sum() / ns).round(1).to_dict())
