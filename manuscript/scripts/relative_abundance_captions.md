# Relative abundance figures — captions

Figures produced by `polish_relative_abundance.R` (run via `render_relative_abundance.bat`).
Outputs: `metaOmics/test_run/Figures/relative_abundance/relative_abundance_{species_top30,genus_top10}.{png,tiff}`;
tables: `metaOmics/test_run/Datasets/relative_abundance/`.

## Species level (`relative_abundance_species_top30`)

**Figure X. Species-level relative abundance of the hospital air and surface microbiome.**
Stacked bars show the relative abundance of the 30 most abundant microbial species in each sample,
grouped by sampling type (Air, n = 23; Surface, n = 35). Taxonomic profiles were obtained with Kraken2;
only species-rank assignments with more than 17 reads per sample were retained, and non-microbial taxa
(land plants and animals) and likely false-positive assignments (*Plasmodium* spp.) were removed based
on their Kraken2 lineage. Species were ranked by mean relative abundance across all samples; remaining
species are pooled as "Other". Species of the same genus share a colour hue. Air samples were dominated
by *Paracoccus marcusii* (mean relative abundance 14.3%), *Rahnella aquatilis* (7.8%) and the fungi
*Cryptococcus neoformans* (6.6%) and *Aspergillus chevalieri* (6.0%), whereas surface samples were
dominated by *Acinetobacter lwoffii* (9.4%), *Stutzerimonas stutzeri* (8.2%) and *Psychrobacter* sp.
van23A (6.9%).

Notes: replace "Figure X". Short-caption option: stop after "share a colour hue" and move the
dominant-species sentences to Results. "Other" ≈ 43% (Air) / 41% (Surface) of reads.
Values from `species_top30_long.csv` (run 2026-09-24).
