# chrXatlas

**A systematic, evidence-scored atlas of chromosome X genetic associations in human traits.**

Chromosome X remains under-analyzed or omitted in many GWAS summary-stat resources. The analytical complications — different copy numbers in males and females, dosage compensation uncertainty, non-standard LD structure — mean that an entire chromosome encoding ~800 protein-coding genes is still handled inconsistently in the era of biobank-scale genomics.

chrXatlas is a chromosome X evidence atlas built from publicly available Pan-UK Biobank summary statistics. It ranks human traits by the strength and quality of their chrX associations, groups significant variants into loci, maps candidate genes, and adds targeted eQTL follow-up. It is a transparent discovery resource, not a sex-stratified, LD-aware, or causal model of chromosome X biology.

**[Browse the atlas](https://smkwray.github.io/xchratlas/)**

---

## Current build (2026-03-19)

| Metric | Value |
|--------|-------|
| Curated panels | 3 (Broad Atlas, Mind & Risk, Biochemistry Deep Dive) |
| Curated traits | 79 across panels |
| Chromosome X loci mapped | 1,109 |
| eQTL lookup-hit coverage | 77 of 79 curated traits |
| Strict eQTL-supported traits | 62 of 79 curated traits |
| Biological domains | 26 unique |
| Discovery pool | 150 traits, 643 loci |

### Curated panels

- **Broad Atlas** — 48 traits, 751 loci, 47 lookup-hit and 40 strictly eQTL-supported. The featured all-purpose panel with broad trait coverage across 21 biological domains.
- **Mind & Risk** — 10 traits, 26 loci, 9 lookup-hit and 4 strictly eQTL-supported. Focused behavior and cognition panel: risk-taking, smoking, alcohol, mood, irritability, reaction time, neuroticism, insomnia, depression, fluid intelligence.
- **Biochemistry Deep Dive** — 21 traits, 332 loci, 21 lookup-hit and 18 strictly eQTL-supported. High-yield blood biochemistry panel covering lipoproteins, kidney markers, endocrine biomarkers, minerals, and liver enzymes.
- **Discovery Pool** — 150 traits, 643 loci, 72 lookup-hit and 34 strictly eQTL-supported. The discovery pool is built from the full Pan-UKB max independent set filtered to `num_pops_pass_qc >= 2`, so the broad build starts from non-redundant traits with at least minimal multi-population QC rather than the noisiest single-population results. It is kept as a curation reservoir, not a featured panel.

---

## Why this project exists

Chromosome X is a persistent blind spot of modern human genetics. It carries genes that influence immune function, metabolism, brain development, blood chemistry, and dozens of other processes. But because males carry one copy (XY) and females carry two (XX), many GWAS pipelines and downstream summary-stat resources still omit chrX or handle it with untested assumptions.

The result is a systematic gap: thousands of GWAS have been published, but chromosome X associations remain poorly cataloged, rarely compared across traits, and difficult to find in existing resources.

chrXatlas addresses this by building a structured, ranked catalog from public data. It does not claim to measure how much of a trait is caused by chromosome X. It identifies **where** on chrX we see statistically significant associations and attempts to connect those signals to candidate genes through gene mapping and eQTL follow-up.

---

## Methods

### Model provenance

chrXatlas does not fit its own GWAS association model. It consumes existing summary statistics from Pan-UKB and adds downstream extraction, locus calling, gene mapping, and evidence scoring. The division of labor:

- **Association model**: inherited from Pan-UKB (linear/logistic mixed models across multiple ancestries)
- **chrXatlas model**: distance-based locus grouping + rule-based gene mapping + heuristic evidence score + targeted eQTL lookup
- **Not modeled here**: sex-stratified effects, dosage compensation, X-inactivation, LD-aware fine-mapping, colocalization, tissue-relevance weighting

### Overview

The pipeline takes publicly available GWAS summary statistics, extracts chromosome X associations, calls loci, maps candidate genes, and scores each trait by the strength and quality of its chrX evidence. Every step is designed to be conservative and transparent.

### Data source

All GWAS data comes from the **Pan-UK Biobank** (Pan-UKB) project, which provides multi-ancestry GWAS summary statistics for thousands of traits. Pan-UKB per-phenotype files use **GRCh37/hg19** coordinates. chrXatlas uses sex-combined summary statistics (`pheno_sex = both_sexes`).

### Trait selection

Traits are selected from the Pan-UKB phenotype manifest using structured queries against `description`, `category`, `trait_type`, `phenocode`, and `modifier`. Selection strongly prefers:

- Membership in the Pan-UKB max independent set (`in_max_independent_set = True`)
- Larger high-quality cohort sample sizes
- Multi-population QC passing (`num_pops_pass_qc >= 2`)
- Clear, common trait definitions

Traits are organized into curated panels by biological domain. Each panel is defined by a CSV config file specifying trait IDs and domain labels.

### Chromosome X extraction

For each trait, the pipeline extracts only chromosome X rows from the Pan-UKB per-phenotype files using remote tabix queries (or HTTP streaming fallback). Only a reduced column subset is retained: chromosome, position, alleles, rsID/varid, and the best available meta-analysis effect size, standard error, and p-value columns. Full source files are never stored locally.

### P-value handling

Pan-UKB changed its p-value schema over time. The pipeline detects the actual column format and scale:

1. `neglog10_pval_meta_hq` (preferred)
2. `neglog10_pval_meta` (fallback)
3. Raw `pval_*` columns (if present)

Both the original column name and inferred scale (`raw`, `neglog10`, or `ln`) are recorded in outputs.

### Chromosome X region model

Every variant, locus, and gene is labeled by its X-chromosome region using GRCh37 boundaries:

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 60,001 | 2,699,520 |
| nonPAR | 2,699,521 | 154,931,043 |
| PAR2 | 154,931,044 | 155,260,560 |

Loci are never merged across region boundaries. PAR regions behave like autosomes; nonPAR regions have the sex-linked dosage differences that make chrX analysis complicated.

### Locus calling

Loci are called using a simple, transparent distance-based strategy:

1. Retain only variants passing genome-wide significance (p < 5 &times; 10<sup>-8</sup>) on chromosome X
2. Sort by region then position
3. Merge significant variants into a locus when within a **500 kb** window
4. Never merge across PAR1/nonPAR/PAR2 boundaries

This is deliberately conservative and easy to explain. Each locus records its genomic interval, lead variant (position, alleles, rsID, p-value), and region label.

### Candidate gene mapping

For each locus, candidate genes are mapped using the GRCh37 Ensembl chrX gene catalog with the following priority hierarchy:

1. **Lead overlap** — gene overlaps the lead variant position
2. **Locus overlap** — gene overlaps the locus interval
3. **Nearest gene** — nearest gene to the lead variant (within ±100 kb, then beyond)

The mapping relation is recorded explicitly (e.g., `lead_overlap`, `locus_overlap`, `nearest_42317bp`). Protein-coding genes are preferred in ranking, but all gene biotypes are retained with their mapping relation visible.

### eQTL follow-up

eQTL evidence is incorporated through targeted queries against the **eQTL Catalogue** using a conservative lookup order:

1. **rsID-first** — if a lead variant has an rsID, query the eQTL Catalogue dataset-scoped REST API by rsID
2. **Variant recoder rescue** — if no rsID is available but GRCh37 alleles are, attempt rsID recovery through the GRCh37 Ensembl `variant_recoder`
3. **Stop if no rsID** — if no rsID can be recovered, mark the locus as not safely assessed (coordinate-based region queries are provisional and explicitly flagged)

For each locus, the pipeline records:
- Whether the lookup returned any eQTL associations (`eQTL lookup-hit`)
- Whether any positional candidate gene passes the strict support rule (`eQTL-supported`)
- Which studies and datasets contributed evidence
- Best observed eQTL p-value
- Lookup mode used (`rsid`, `variant_recoder+rsid`, or `region_provisional`)

A locus is **eQTL lookup-hit** when its rsID-based follow-up against prioritized eQTL Catalogue datasets returns one or more associations. A locus is **eQTL-supported** only when at least one pre-defined candidate gene at that locus has aggregated eQTL evidence with `best_eqtl_pvalue <= 1e-5` in the current build. A trait is **eQTL-supported** when at least one of its chrX loci has eQTL-supported candidate-gene evidence.

Note: the current pipeline does not perform chrX LD modeling, colocalization analysis, or trait-specific tissue relevance matching. Study selection uses a prioritized dataset list (see `config/eqtl_priority_studies.csv`), not a tissue-relevance model. LD-aware colocalization and tissue-specific weighting are future work.

### Evidence scoring

Each trait receives an **X-evidence score** (0–100) combining multiple weighted components:

| Component | Weight | What it measures |
|-----------|--------|-----------------|
| Independent set membership | 12 | Trait is in the Pan-UKB max independent set |
| Multi-population QC | 8 | Passes QC in 2+ populations |
| Large high-quality cohort | 10 | Sample size in the high-quality cohort |
| Any significant locus | 12 | At least one genome-wide significant chrX locus |
| Per-extra locus | 4 (capped at 16) | Additional loci beyond the first |
| Lead signal strength | up to 20 | Strength of the lead variant p-value |
| Nearest gene mapping | 8 | Lead variant maps to a nearby gene |
| Overlapping gene | 12 | Lead variant falls within a gene |
| eQTL support | 15 | At least one locus has eQTL evidence |
| nonPAR bonus | 10 | Signal observed outside pseudo-autosomal regions |
| PAR-only bonus | 3 | Signal observed only in pseudo-autosomal regions |
| No-locus penalty | -10 | No genome-wide significant loci found |

The score is a **ranking aid**, not a biological measurement. It is normalized to 0–100 and accompanied by confidence notes.

### Coverage grades

Each trait receives a **coverage grade** reflecting the quality and completeness of its upstream GWAS metadata, computed from Pan-UKB trait-level fields (`scoring.py:coverage_grade`):

- **Grade A** — trait is in the Pan-UKB max independent set, passes QC in 2+ populations, and has a high-quality cohort of at least 50,000
- **Grade B** — trait is in the independent set with a cohort of at least 10,000, or passes QC in 2+ populations
- **Grade C** — high-quality cohort of at least 1,000
- **Grade U** — does not meet the above thresholds; coverage is uncertain

The coverage grade reflects **GWAS data quality**, not eQTL follow-up completeness. A trait can have Grade A coverage (strong upstream GWAS) but no strict eQTL support.

### Gene annotation

The chromosome X gene catalog is built from the **GRCh37 Ensembl REST API** (`grch37.rest.ensembl.org`) to ensure coordinate alignment with Pan-UKB positions. Gene records include Ensembl ID, gene name, biotype, and genomic coordinates.

---

## What the atlas is not

- It does **not** estimate the fraction of trait variance explained by chromosome X.
- It does **not** perform sex-stratified analysis or model dosage compensation explicitly.
- It uses **sex-combined** GWAS summary statistics from Pan-UKB.
- It does **not** claim autosomal tools like LDSC apply to chrX.
- Distance-based locus merging is a conservative v1 strategy, not LD-aware fine-mapping.

The atlas identifies where on chromosome X we see associations and connects those signals to candidate genes. Interpreting the biological importance of those associations requires additional context.

---

## Frontend

The interactive frontend lives in `site/` and reads the JSON data bundle from `site/data/`.

### Pages

- **Homepage** (`index.html`) — Hero, metrics, curated panel cards, domain composition, top traits rail, panel comparison, discovery pool, methods explainer, search
- **Panel detail** (`panel.html?id=...`) — Full trait table with column sorting, domain filter chips, text search, eQTL-supported-only toggle
- **Trait detail** (`trait.html?panel=...&trait=...`) — Evidence score, eQTL lookup-hit versus support explanation card, trait metadata grid, expandable locus cards with candidate gene tables, notes

### Features

- Light/dark theme with system preference detection, manual override, localStorage persistence
- Cross-panel trait search from the homepage
- SVG sun/moon toggle icon in the header
- Colored panel quick-nav buttons (Broad Atlas, Mind & Risk, Biochemistry) in the header
- Responsive design (desktop + mobile)
- No build step, no framework dependencies — plain HTML/CSS/JS reading static JSON

### Running locally

```bash
cd site && python3 -m http.server 8891
```

### Data refresh

The frontend reads from `site/data/`, which is a copy of the backend's `data/release/frontend/` output. To refresh after a backend rebuild:

```bash
cp -R "$PROJ_SHARED_DATA_ROOT/release/frontend/." site/data/
```

---

## Backend pipeline

### Scripts (in order)

| Script | Purpose |
|--------|---------|
| `00_fetch_panukb_manifests.py` | Download Pan-UKB phenotype manifest and eQTL metadata |
| `01_select_seed_traits.py` | Select traits from manifest using panel config |
| `02_extract_panukb_x.py` | Extract chrX-only GWAS slices via remote tabix |
| `03_build_x_gene_catalog.py` | Build GRCh37 chrX gene catalog from Ensembl |
| `04_prepare_eqtl_index.py` | Prepare eQTL Catalogue study index |
| `05_call_x_loci.py` | Call chrX loci with distance-based merging |
| `06_fetch_eqtl_region_hits.py` | Targeted eQTL follow-up around lead loci |
| `07_map_loci_to_genes.py` | Map loci to candidate genes |
| `08_export_release_tables.py` | Export release TSV tables |
| `09_compare_release_panels.py` | Compare release panels |
| `10_prepare_independent_set_panel.py` | Build the broader discovery pool |
| `11_split_selected_traits.py` | Split trait selection by panel shards |
| `12_merge_panel_shards.py` | Merge panel shards |
| `13_generate_static_site.py` | Generate static HTML site |
| `14_export_frontend_bundle.py` | Export JSON bundle for the interactive frontend |

### Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .

python scripts/00_fetch_panukb_manifests.py
python scripts/01_select_seed_traits.py --panel config/panel_expanded.csv
python scripts/02_extract_panukb_x.py
python scripts/03_build_x_gene_catalog.py
python scripts/04_prepare_eqtl_index.py          # required before eQTL follow-up
python scripts/05_call_x_loci.py
python scripts/06_fetch_eqtl_region_hits.py
python scripts/07_map_loci_to_genes.py
python scripts/08_export_release_tables.py
python scripts/14_export_frontend_bundle.py
```

### Release outputs

```
data/release/
├── trait_scores.tsv
├── x_loci.tsv
├── x_gene_candidates.tsv
├── manifest_snapshot.tsv
├── panel_expanded/
├── panel_mind_risk_core/
├── panel_blood_biochemistry_core/
├── panel_independent_set/
├── panel_release_comparison.tsv
└── frontend/
    ├── manifest.json
    └── panels/{panel_id}/summary.json, traits.json, traits/*.json
```

---

## Data sources

| Source | Artifact / Endpoint | Role | Build | Used by |
|--------|---------------------|------|-------|---------|
| [Pan-UKB](https://pan.ukbb.broadinstitute.org) | `phenotype_manifest.tsv.bgz` | Trait metadata, selection, manifest fields | GRCh37 | `00`, `01` |
| Pan-UKB | Per-phenotype summary stat files (remote tabix) | chrX GWAS extraction | GRCh37 | `02` |
| Pan-UKB | `full_variant_qc_metrics.txt.bgz` (optional) | Variant-level QC reference | GRCh37 | — |
| [Ensembl GRCh37 REST](https://grch37.rest.ensembl.org) | `/overlap/region/human/X:{start}-{end}?feature=gene` | chrX gene catalog | GRCh37 | `03` |
| Ensembl GRCh37 REST | `/variant_recoder/human/{varid}` | rsID recovery for lead variants lacking rsIDs | GRCh37 | `06` |
| [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) | chrX genotypes metadata (`chrX_genotypes.tsv`) | Identify chrX-capable eQTL studies | — | `04` |
| eQTL Catalogue | Tabix FTP paths metadata (`tabix_ftp_paths.tsv`) | Study/dataset index for lookups | — | `04` |
| eQTL Catalogue | Dataset-scoped REST API (rsID queries) | Targeted eQTL follow-up per locus | rsID-based | `06` |

See `docs/DATA_SOURCES.md` for full endpoint URLs and access notes.

---

## Storage

The sub-50 GB constraint is met by:

1. Never storing full Pan-UKB summary-stat files locally
2. Using remote tabix on per-phenotype files
3. Keeping only reduced chrX column subsets
4. Doing targeted eQTL queries, not bulk downloads

---

## Project structure

```
chrXatlas/
├── src/xatlas/          # Python package (scoring, loci, eQTL, site generation, bundle export)
├── scripts/             # Numbered pipeline scripts (00–14)
├── config/              # Panel definitions, scoring weights, region boundaries
├── tests/               # Test suite
├── site/                # Interactive frontend (HTML/CSS/JS + JSON data)
│   ├── index.html       # Landing page
│   ├── panel.html       # Panel detail page
│   ├── trait.html       # Trait detail page
│   ├── css/style.css    # Design system
│   ├── js/              # theme.js, app.js, panel.js, trait.js
│   └── data/            # Copied frontend JSON bundle
├── docs/                # Project documentation
└── data/                # Pipeline data (raw, interim, processed, release)
```

---

## Disclaimer and attribution

This project is independent research. It is not affiliated with or endorsed by UK Biobank, the Pan-UKB team, the eQTL Catalogue, or any of the data providers.

No software license file is included in this repository yet, so reuse terms are not currently specified.

Data sources: Pan-UK Biobank, eQTL Catalogue, Ensembl GRCh37.
