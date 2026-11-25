# chrXatlas

**[Browse the interactive atlas](https://smkwray.github.io/xchratlas/)**

chrXatlas summarizes chromosome X associations across selected UK Biobank traits using public Pan-UK Biobank (Pan-UKB) summary statistics. It groups significant variants by distance, maps Ensembl GRCh37 positional candidate genes, attaches eQTL Catalogue lookups, and computes a 0–100 ranking score. The repository is `xchratlas` and the Python package is `xatlas`.

## Snapshot

The saved snapshot (exported 19 March 2026) contains three overlapping curated panels and a 150-trait discovery pool.

| Panel | Traits | Trait–locus entries | Traits with eQTL records | eQTL-supported traits |
|---|---:|---:|---:|---:|
| Broad Atlas | 48 | 751 | 47 | 40 |
| Mind and Behavior | 10 | 26 | 9 | 4 |
| Biochemistry Deep Dive | 21 | 332 | 21 | 18 |
| Curated panel entries | 79 | 1,109 | 77 | 62 |
| Discovery Pool | 150 | 643 | 72 | 34 |

The curated panels hold 57 distinct traits across 79 panel entries and 26 trait categories (21 in Broad Atlas); traits and their loci can appear in more than one panel. Mind and Behavior covers risk-taking, smoking, alcohol, mood, irritability, reaction time, neuroticism, insomnia, frequency of depressed mood in the past two weeks, and fluid intelligence. Biochemistry Deep Dive covers lipoproteins, kidney markers, endocrine biomarkers, minerals, and liver enzymes. The Discovery Pool holds 150 Pan-UKB max independent-set traits passing QC in at least two ancestry groups, including 66 with no significant chromosome X loci. In Broad Atlas, standing height has 51 loci, body mass index 37, insulin-like growth factor 1 has 36, and whole body fat mass 35.

## Method

- **Loci & threshold:** Pan-UKB `both_sexes` GRCh37 summary statistics are filtered at p < 5 × 10⁻⁸ per trait, with no correction across traits. Consecutive significant variants within 500 kb are chained into a locus within PAR1, non-PAR, and PAR2 separately; these are distance groups, not LD-independent signals.
- **Candidate genes:** Up to five positional candidate genes per locus are selected from Ensembl GRCh37: genes overlapping the locus, otherwise within 100 kb of the lead variant, otherwise the nearest genes. Genes returned only by eQTL lookups are labeled `Follow-up` and excluded from trait-level support.
- **eQTL lookups & support:** Lead variants are queried by rsID (with Ensembl GRCh37 rsID recovery when needed, without allele-level matching) in a fixed priority list of eQTL Catalogue datasets. **Returned records** means at least one association was returned. **eQTL support** means at least one positional candidate gene has an eQTL p ≤ 10⁻⁵ (nominal, uncorrected, no colocalization or tissue-relevance model).
- **Coverage grades:** Grades A, B, C, and U rate upstream Pan-UKB metadata (not eQTL completeness) by first match, where N is HQ cases for binary traits or HQ measured participants for quantitative traits: **A** = max independent set, QC in ≥2 populations, N ≥ 50,000; **B** = max independent set with N ≥ 10,000, or QC in ≥2 populations; **C** = N ≥ 1,000; **U** = below A–C.

The evidence score sums the weights below (`config/scoring_weights.yml`) and caps the total at 0–100.

| Component | Points |
|---|---:|
| Trait is in the Pan-UKB max independent set | 12 |
| Trait passes QC in two or more populations | 8 |
| At least 50,000 HQ cases (binary) or measured participants (quantitative) | 10 |
| At least one significant chromosome X locus | 12 |
| Each additional locus | 4, up to 16 |
| Lead p-value beyond genome-wide significance | 1 per order of magnitude, up to 20 |
| Any positional candidate overlaps the lead variant or locus (nearby or nearest only: 8) | 12 |
| At least one eQTL-supported locus | 15 |
| A locus outside the pseudoautosomal regions (pseudoautosomal only: 3) | 10 |
| No significant locus | −10 |

## Data sources

| Source | Used for |
|---|---|
| [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org) | Trait metadata, trait selection, and GRCh37 chromosome X summary statistics |
| [Ensembl GRCh37 REST API](https://grch37.rest.ensembl.org) | Chromosome X gene catalog and rsID recovery |
| [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) | Dataset index and historical eQTL lookups |

## Preview the saved site

```bash
python3 -m http.server --directory site 8891
```

Open `http://localhost:8891/`.

## Reproduction and provenance

The [eQTL Catalogue REST API](https://www.ebi.ac.uk/eqtl/Data_access/) used for the original lookups (scripts `04` and `06`) has been retired, so the pipeline cannot be rerun end-to-end against the live API. The saved site data in `site/data/` is the reference snapshot.

From compatible cached processed inputs for each panel (Python 3.10+):

```bash
panel_id=panel_expanded  # set to the panel represented by the cached inputs
python scripts/08_export_release_tables.py --outdir "data/release/${panel_id}"
```

With release tables for all four panels (`panel_expanded`, `panel_mind_risk_core`, `panel_blood_biochemistry_core`, `panel_independent_set`):

```bash
python scripts/14_export_frontend_bundle.py
cp -R data/release/frontend/. site/data/
python scripts/15_build_gene_index.py
```

## Project structure

```
src/xatlas/     Python package: loci, regions, gene mapping, eQTL lookups, scoring, bundle export
scripts/        Pipeline, release export, and gene-index scripts (00–15)
config/         Panel trait lists, scoring weights, region boundaries, eQTL dataset priorities
site/           Static site (HTML, CSS, JavaScript) and saved JSON data
docs/           Frontend data contract
tests/          pytest suite
data/           Local pipeline data (not committed)
```

[`docs/FRONTEND_DATA_CONTRACT.md`](docs/FRONTEND_DATA_CONTRACT.md) describes the JSON bundle schema.

## Limitations

- Source files are Pan-UKB `both_sexes` summary statistics (some phenotypes are sex-limited); sex differences, X inactivation, and dosage compensation are not tested or modeled.
- Pan-UKB samples are mostly of European ancestry, and contributing ancestries vary by trait.
- Three saved loci reached the p-value conversion floor and display `−log10(p) ≥ 300` with a recorded representative variant; gene mapping and eQTL results refer to that representative.
- 26 of 48 Broad Atlas traits and 19 of 21 Biochemistry Deep Dive traits reach the 100-point score cap; ties are ordered by locus count and returned eQTL records, which does not represent genetic effect size.

## Disclaimer

Independent research, unaffiliated with or endorsed by UK Biobank, the Pan-UKB team, the eQTL Catalogue, or Ensembl.

## License

MIT. See [LICENSE](LICENSE). Source data remain under their providers' terms.
