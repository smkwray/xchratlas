# chrXatlas

**[Browse the interactive atlas](https://smkwray.github.io/xchratlas/)**

chrXatlas ranks human traits by the strength and quality of their genetic associations on chromosome X, which many genome-wide association study (GWAS) resources omit or handle inconsistently. It extracts chromosome X results from public Pan-UK Biobank (Pan-UKB) GWAS summary statistics, groups significant variants into loci, maps candidate genes with Ensembl GRCh37 annotations, and looks up expression quantitative trait locus (eQTL) evidence in the eQTL Catalogue. The repository is `xchratlas` and the Python package is `xatlas`.

It is a discovery and ranking catalog built on Pan-UKB's sex-combined results. It does not fit its own association models, analyze the sexes separately, or estimate how much of a trait chromosome X explains, and its score is a ranking aid, not a biological measurement.

## Current build

The build of 19 March 2026 has three curated panels and a discovery pool.

| Panel | Traits | Loci | eQTL lookup-hit traits | eQTL-supported traits |
|---|---:|---:|---:|---:|
| Broad Atlas | 48 | 751 | 47 | 40 |
| Mind & Behavior | 10 | 26 | 9 | 4 |
| Biochemistry Deep Dive | 21 | 332 | 21 | 18 |
| Curated panels | 79 | 1,109 | 77 | 62 |
| Discovery Pool | 150 | 643 | 72 | 34 |

The curated panels span 26 biological domains, 21 of them in the Broad Atlas. Mind & Behavior covers risk-taking, smoking, alcohol, mood, irritability, reaction time, neuroticism, insomnia, depression and fluid intelligence; the Biochemistry Deep Dive covers lipoproteins, kidney markers, endocrine biomarkers, minerals and liver enzymes. The Discovery Pool, a source for curation, holds the Pan-UKB max independent set of traits that pass quality control (QC) in at least two populations. In the Broad Atlas, standing height has the most chromosome X loci (51), followed by body mass index (37), insulin-like growth factor 1 (36) and whole body fat mass (35).

## Method

The pipeline reads only the chromosome X rows of each trait's Pan-UKB summary statistics (sex-combined meta-analysis, GRCh37) through remote tabix queries. Significant variants (p < 5 × 10⁻⁸) no more than 500 kb apart are chained into one locus, never across the boundaries between the two pseudoautosomal regions (PAR1, PAR2) and the rest of the chromosome. Candidate genes come from the Ensembl GRCh37 catalog, preferring a gene that overlaps the lead variant, then the locus, then the nearest gene. Lead variants are looked up by rsID (recovered through Ensembl's variant recoder when missing) in a prioritized list of eQTL Catalogue datasets. A locus is an eQTL lookup-hit when that lookup returns any association, and eQTL-supported when one of its candidate genes has an eQTL p-value of 1 × 10⁻⁵ or lower. Each trait also gets a coverage grade (A, B, C or U) that rates its upstream GWAS metadata, not its eQTL follow-up.

The X-evidence score adds the points below (weights in `config/scoring_weights.yml`) and is clipped to 0–100.

| Component | Points |
|---|---:|
| Trait is in the Pan-UKB max independent set | 12 |
| Trait passes QC in two or more populations | 8 |
| High-quality cohort of 50,000 or more | 10 |
| At least one significant chromosome X locus | 12 |
| Each additional locus | 4, up to 16 |
| Lead p-value beyond genome-wide significance | 1 per order of magnitude, up to 20 |
| A candidate gene overlaps the lead variant or locus (nearest gene only: 8) | 12 |
| At least one eQTL-supported locus | 15 |
| A locus outside the pseudoautosomal regions (pseudoautosomal only: 3) | 10 |
| No significant locus | −10 |

## Data sources

| Source | Used for |
|---|---|
| [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org) phenotype manifest and per-phenotype summary statistics | Trait metadata, trait selection and chromosome X associations |
| [Ensembl GRCh37 REST API](https://grch37.rest.ensembl.org) | Chromosome X gene catalog and rsID recovery |
| [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) | Dataset index and eQTL lookups |

The pipeline never stores full summary-statistic files; it keeps a reduced set of columns for chromosome X rows.

## Reproduce

Requires Python 3.10 or later and network access to Pan-UKB, Ensembl and the eQTL Catalogue.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .

python scripts/00_fetch_panukb_manifests.py
python scripts/01_select_seed_traits.py --panel config/panel_expanded.csv
python scripts/02_extract_panukb_x.py
python scripts/03_build_x_gene_catalog.py
python scripts/04_prepare_eqtl_index.py
python scripts/05_call_x_loci.py
python scripts/06_fetch_eqtl_region_hits.py
python scripts/07_map_loci_to_genes.py
python scripts/08_export_release_tables.py
pytest
```

Pipeline data go to `data/`. These commands build one panel from its trait list in `config/`; script `10` prepares the Discovery Pool trait list, and scripts `11` and `12` shard a large panel and merge the results. `scripts/14_export_frontend_bundle.py` reads one release directory per panel (by default `data/release/<panel_id>/`) and writes the site's JSON bundle to `data/release/frontend/`. To update and preview the site:

```bash
cp -R data/release/frontend/. site/data/
python3 -m http.server -d site 8891
```

## Project structure

```
src/xatlas/     Python package: loci, regions, gene mapping, eQTL lookups, scoring, bundle export
scripts/        Numbered pipeline scripts (00–14)
config/         Panel trait lists, scoring weights, region boundaries, eQTL dataset priorities
site/           Interactive atlas (HTML, CSS, JavaScript) and its JSON data
docs/           Frontend data contract
tests/          pytest suite
data/           Pipeline data (not committed)
```

Each release directory holds `trait_scores.tsv`, `x_loci.tsv` and `x_gene_candidates.tsv`; [`docs/FRONTEND_DATA_CONTRACT.md`](docs/FRONTEND_DATA_CONTRACT.md) describes the JSON bundle.

## Limitations

- Loci come from distance-based merging, not linkage-disequilibrium-aware fine-mapping, and there is no colocalization analysis.
- Dosage compensation and X inactivation are not modeled; pseudoautosomal-only signals should be read with care.
- eQTL datasets come from a fixed priority list, not a model of which tissues matter for each trait, and loci whose lead variant has no recoverable rsID are marked as not safely assessed.
- 22 of 48 Broad Atlas traits and 19 of 21 Biochemistry Deep Dive traits reach the 100-point score cap, so the score does not order the strongest traits.

## Disclaimer and license

This is independent research, not affiliated with or endorsed by UK Biobank, the Pan-UKB team, the eQTL Catalogue or any other data provider. No software license file is included yet, so reuse terms are not specified.
