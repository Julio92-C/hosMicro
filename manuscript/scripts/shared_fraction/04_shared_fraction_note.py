"""Collaborator note: the air-surface shared fraction (docx)."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]          # manuscript/
DATA = str(ROOT / 'data') + '/'
OUTD = str(ROOT / 'notes' / 'shared_fraction') + '/'
Path(OUTD).mkdir(parents=True, exist_ok=True)

import re, pandas as pd
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
D = OUTD
OUT = D + 'Shared_fraction_note_2026-10-29.docx'
ws = pd.read_csv(D + 'shared_fraction_widespread_species.csv', index_col=0)

doc = Document(); sec = doc.sections[0]
sec.top_margin = sec.bottom_margin = Cm(1.8); sec.left_margin = sec.right_margin = Cm(2)
st = doc.styles['Normal']; st.font.name = 'Calibri'; st.font.size = Pt(10.5)
ITAL = re.compile(r'\b((?:Rahnella|Paracoccus|Pseudomonas|Stutzerimonas|Acinetobacter|Stenotrophomonas|Escherichia|Serratia|Pantoea|Chlamydomonas|'
                  r'Microbacterium|Gordonia|Mycolicibacterium|Rhodococcus|Nocardia|Psychrobacter)(?: [a-z]+)?|[A-Z]\. [a-z]{3,})\b')
def add(p, text, bold=False):
    pos = 0
    for m in ITAL.finditer(text):
        if m.start() > pos: r = p.add_run(text[pos:m.start()]); r.bold = bold
        r = p.add_run(m.group(0)); r.italic = True; r.bold = bold; pos = m.end()
    if pos < len(text): r = p.add_run(text[pos:]); r.bold = bold
def para(text, bold=False, size=None, after=4, style=None):
    p = doc.add_paragraph(style=style); p.paragraph_format.space_after = Pt(after); add(p, text, bold)
    if size:
        for r in p.runs: r.font.size = Pt(size)
    return p
def head(text): p = para(text, bold=True, size=12, after=3); p.paragraph_format.space_before = Pt(8)

para('The air–surface shared fraction: a closer look', bold=True, size=15, after=2)
para('Hospital microbiome manuscript (v2) — note for co-authors, 29 October 2026. For discussion; not yet in the manuscript.', size=9.5, after=8)

head('Why this note')
para('The manuscript reports that 590 of 2,255 species (26.2%) were detected in both air and surface samples (Figure 1C). This shared fraction '
     'is larger than the "ecologically distinct" framing suggests, so we broke it down by how often each species occurs in each compartment, '
     'what it is, and where in the hospital it was found. The same species counts and filters as the manuscript are used (23 air and 35 surface '
     'samples; species with more than 10 reads; 696 air-only / 590 shared / 969 surface-only).')

head('Main observations')
for b in [
    'Most of the shared fraction is sporadic. 496 of the 590 shared species occur in only one sample of at least one compartment, and no '
    'species is found in half or more of both air and surface samples. About half of these sporadic species are Actinomycetota '
    '(Microbacterium, Gordonia, Mycolicibacterium, Rhodococcus, Nocardia), typical dust and built-environment genera (Figure 1C).',
    'Its large share of reads reflects dominant species detected at trace levels in the other compartment. Shared species account for 83% of '
    'air reads and 55% of surface reads, but mostly because, for example, R. bonaserana (28.6% of air, 0.01% of surface reads) and A. lwoffii '
    '(8.9% of surface, 0.02% of air) were each detected once or a few times in the other compartment (Figure 1A–B).',
    'A small widespread shared fraction stands out. Thirteen species occur in at least five air and five surface samples (Figure 1B, D; Table 1). '
    'Twelve are Gram-negative bacteria, mostly opportunistic or hospital-associated: Pseudomonas aeruginosa, P. putida, P. fluorescens, '
    'Stutzerimonas stutzeri, Stenotrophomonas maltophilia, Acinetobacter baumannii, Escherichia coli, Serratia marcescens, S. proteamaculans, '
    'S. quinivorans, Pantoea agglomerans and Rahnella aceris; the thirteenth is the alga Chlamydomonas reinhardtii. Together they make up '
    '6.4% of air and 18.2% of surface reads.',
    'Two patterns within it. The Pseudomonas group (P. aeruginosa, P. putida, P. fluorescens, S. stutzeri) is surface-led (57–66% of surface '
    'samples, 22–39% of air) and carries most of the virulence genes (S. stutzeri: 42 VF genes). The Enterobacterales (R. aceris, '
    'S. marcescens, P. agglomerans, E. coli) are air-led (57–65% of air samples, 14–29% of surfaces).',
    'Hospital linkage (Figure 2). These species link air from most wards to almost every ward E room: 82 of the 110 possible air-ward / '
    'surface-room pairs share at least one of them. The strongest links (up to 12 species) connect air from wards D, G, H and I with the ward E '
    'patient rooms (SR 16, SR 26), the kitchen and the beverage bay. Air from ward F carried none of them, and ward E air (one sample) only one.',
    'Sample-level view (Figure 3). The Pseudomonas group is spread across almost all ward E rooms and appears in air mainly from wards D, G, H '
    'and I, whereas the Enterobacterales (Serratia, Pantoea, Rahnella) occur in air from nearly every ward but on surfaces at only a handful of '
    'sites: the three Serratia species in the same five samples (the SR 16 and beverage-bay fridges, an SR 26 room vent and a corridor vent) and '
    'P. agglomerans mainly on kitchen vents. A. lwoffii and Paracoccus marcusii, although shared, are clearly anchored in one compartment '
    '(surfaces and air, respectively).',
]:
    p = para('', after=3, style='List Bullet'); add(p, b)

head('Figures')
doc.add_picture(D + 'shared_fraction_prevalence_abundance.png', width=Cm(17))
para('Figure 1. The air–surface shared fraction. (A) Prevalence across all 58 samples against mean relative abundance for every species, coloured '
     'by Venn category (air only, surface only, shared); the 13 widespread shared species are circled. (B) Prevalence in air against prevalence '
     'on surfaces for the 590 shared species; point size reflects mean relative abundance; dotted lines mark five samples in each compartment. '
     'Tiers: widespread (≥5 samples in both), recurrent (2–4 samples in both), sporadic (one sample in at least one compartment). (C) Number of '
     'shared species per tier, by phylum. (D) Prevalence of the 13 widespread shared species in air and on surfaces, with the number of distinct '
     'ARGs, VFs and MGEs detected on their reads.', size=9, after=8)
doc.add_picture(D + 'shared_fraction_ward_map.png', width=Cm(17))
para('Figure 2. Where the widespread shared species occur in the hospital (schematic; positions are not to scale). (A) Air samples grouped by '
     'ward (left) and ward E surface samples grouped by room (right); node size is the number of widespread shared species detected; lines join '
     'an air ward and a ward E room that share at least five of them, with line width proportional to the number shared; n is the number of '
     'samples. (B) Presence of each widespread shared species per air ward and per ward E room.', size=9, after=8)

doc.add_picture(D + 'shared_fraction_heatmap_top30.png', width=Cm(17))
para('Figure 3. Relative abundance of the 30 shared species with the highest prevalence in both compartments, per sample. Air samples are '
     'grouped by ward (left) and ward E surface samples by room (right); grey cells mark no detection (≤10 reads). Species names in red are '
     'the 13 widespread shared species; the other 17 belong to the recurrent tier.', size=9, after=8)

head('Table 1. The 13 widespread shared species')
t = doc.add_table(rows=1, cols=6); t.style = 'Light List Accent 1'
for c, h in zip(t.rows[0].cells, ['Species', 'Air samples (of 23)', 'Surface samples (of 35)', 'Mean abundance air / surface (%)', 'ARG / VF / MGE', 'Example genes']):
    c.text = h
for n, r in ws.iterrows():
    cells = t.add_row().cells
    cells[0].text = ''; cells[0].paragraphs[0].add_run(n).italic = True
    cells[1].text = str(int(r.nA)); cells[2].text = str(int(r.nS))
    cells[3].text = f'{r.aA:.2f} / {r.aS:.2f}'; cells[4].text = f'{int(r.ARGs)} / {int(r.VFs)} / {int(r.MGEs)}'
    cells[5].text = str(r.examples) if isinstance(r.examples, str) else '–'
for row in t.rows:
    for c in row.cells:
        for p in c.paragraphs:
            for run in p.runs: run.font.size = Pt(8)

head('Caveats')
for b in [
    'Gene–taxon links come from Kraken2 classification of element-carrying reads. Some assignments are implausible and probably reflect '
    'misclassification between related taxa or of short element-carrying reads: for example, the Rahnella RAHN-1 β-lactamase on Pantoea and Serratia '
    'reads, and the bacterial regulator CRP on a Chlamydomonas read. The element counts in Table 1 should be read with this in mind.',
    'Presence is based on more than 10 reads per sample; the sporadic tier is sensitive to this floor and to sequencing depth.',
    'Air and surfaces were not sampled symmetrically (surfaces only on ward E; one ward E air sample), so Figure 2 shows co-occurrence, not transfer.',
]:
    p = para('', after=3, style='List Bullet'); add(p, b)

head('Questions for co-authors')
for b in [
    'Is this worth adding to the manuscript? A minimal option is one Results sentence and one Discussion sentence on the widespread shared '
    'fraction, with Figure 1 (or panels B and D) as a new supplementary figure.',
    'Should Figure 2 (ward linkage) and Figure 3 (heatmap) go in as well, or are they better kept for a follow-up on transmission routes?',
    'Is the five-sample threshold for "widespread" reasonable, or would you prefer a prevalence-based cut-off (e.g. ≥20% of samples in both)?',
]:
    p = para('', after=3, style='List Number'); add(p, b)
doc.save(OUT); print('saved', OUT)
