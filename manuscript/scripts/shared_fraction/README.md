# Shared fraction of air and surface communities

Exploratory analysis behind the collaborator note *The air–surface shared fraction* (29 October 2026). It is not in the manuscript yet.

Run in order from any directory (Python 3.12, with pandas, matplotlib, adjustText and python-docx installed):

| Script | Output (in `manuscript/notes/shared_fraction/`) |
|---|---|
| `01_shared_fraction_figure.py` | `shared_fraction_prevalence_abundance.png`; tables of all shared species and of the 13 widespread species with the genes they carry |
| `02_shared_fraction_ward_map.py` | `shared_fraction_ward_map.png` (air by ward linked to ward E surfaces by room); per-site presence table |
| `03_shared_fraction_heatmap.py` | `shared_fraction_heatmap_top30.png`; top-30 table |
| `04_shared_fraction_note.py` | `Shared_fraction_note_2026-10-29.docx` |

Inputs are the processed tables in `manuscript/data/`: the >10-read taxon table, the analysed sample list (`alpha_diversity.csv`), the element-carrying reads (`abri_kraken2_cleaned.csv`) and the metadata. Tiers: widespread = at least 5 samples in each compartment; recurrent = 2–4 in each; sporadic = one sample in at least one compartment. Colours follow the manuscript (air `#828385`, surface `#ed8cb6`).
