# Frontend Data Contract

This document defines the JSON bundle contract produced by:

`scripts/14_export_frontend_bundle.py`

The bundle is additive. Existing release TSV files remain canonical for scientific auditing.

## Output root

Default output root:

`../data/release/frontend/`

Layout:

- `manifest.json`
- `panels/<panel_id>/summary.json`
- `panels/<panel_id>/traits.json`
- `panels/<panel_id>/traits/<trait_id>.json`

## Schema versioning

- `schema_version` is currently `2`.
- `schema_version` is included in:
  - `manifest.json`
  - each panel `summary.json`
  - each trait-detail JSON
- If shape changes are not backward compatible, increment `schema_version`.

## Panels in v1

V1 exports exactly these panel IDs:

- `panel_expanded` (`featured: true`)
- `panel_mind_risk_core` (`featured: true`)
- `panel_blood_biochemistry_core` (`featured: true`)
- `panel_independent_set` (`featured: false`)

`panel_independent_set` is intentionally marked non-featured so frontend can demote it by default.

## Field semantics

## `manifest.json`

- `generated_at`: UTC ISO timestamp
- `panels`: one entry per exported panel
- each panel entry includes:
  - `panel_id`, `slug`, `title`, `featured`, `description`
  - `trait_count`, `lookup_hit_trait_count`, `supported_trait_count`, `zero_locus_trait_count`, `locus_count`
  - `files` with relative paths to:
    - `summary`
    - `traits`
    - `trait_details_dir`

## `summary.json`

- `panel`: panel metadata copied from manifest entry
- `metrics`: aggregate counts
- `domain_summary`: grouped by `domain`
- `top_traits`: top compact trait rows
- `unsupported_traits`: traits without strict eQTL support

## `traits.json`

Compact list for grid/table UIs. Each row includes:

- `trait_id`
- `query_id`
- `slug`
- `description`
- `domain`
- `x_evidence_score`
- `coverage_grade`
- `n_loci`
- `eqtl_lookup_hit`
- `eqtl_lookup_hit_locus_count`
- `eqtl_lookup_hit_count`
- `eqtl_supported`
- `eqtl_supported_locus_count`
- `eqtl_total_hit_count`
- `top_candidate_genes` (array)
- `confidence_notes`
- `eqtl_lookup_note`

## `traits/<trait_id>.json`

- `trait`: compact trait row plus phenotype metadata (`trait_type`, `phenocode`, `pheno_sex`, `coding`, `modifier`, `category`)
- `loci`: sorted by strongest signal first
- each locus includes:
  - `locus_id`, `x_region`, `locus_start`, `locus_end`
  - `lead_pos`, `lead_rsid`, `lead_varid`, `lead_neglog10_pvalue`
  - `top_gene_id`, `top_gene_name`, `best_mapping_relation`
  - `eqtl_lookup_hit`, `eqtl_lookup_hit_gene_count`
  - `eqtl_supported`, `eqtl_supported_gene_count`
  - `eqtl_lookup_status`, `eqtl_lookup_mode`, `eqtl_lookup_n_hits`
  - `candidate_genes` array
- each candidate gene includes:
  - `candidate_rank`, `gene_id`, `gene_name`, `gene_biotype`
  - `mapping_relation`, `distance_to_lead_bp`, `eqtl_gene_role`
  - `eqtl_lookup_hit`, `eqtl_lookup_hit_count`
  - `eqtl_supported`, `eqtl_study_count`, `eqtl_dataset_count`, `best_eqtl_pvalue`

## Evidence semantics

- `eqtl_lookup_hit` means the prioritized rsID-based follow-up returned one or more associations.
- `eqtl_supported` is stricter and means at least one candidate gene passed the build's current support rule after aggregation across queried datasets.
- Trait-level `eqtl_supported` does not imply the positional top gene is supported.
- `eqtl_gene_role = candidate` means the gene came from the positional mapping hierarchy.
- `eqtl_gene_role = followup` means the gene was added only because it appeared in eQTL lookup results.

## Normalization rules

- Booleans are JSON booleans, not string literals.
- Missing numeric/string values are `null`.
- `top_candidate_genes` is always an array (split from comma-delimited TSV field).
- Sorting:
  - traits: `x_evidence_score desc`, `n_loci desc`, `eqtl_total_hit_count desc`, then stable string tie-breakers
  - loci in trait detail: `lead_neglog10_pvalue desc`, `n_sig_variants desc`, then genomic start
  - candidate genes per locus: `candidate_rank asc`, then lexical tie-breakers
