# chrXatlas frontend

Interactive static site for the Chromosome X Evidence Atlas.

## Pages

- `index.html` — Landing page with hero, metrics, curated panels, domain composition, top traits, search
- `panel.html?id=...` — Panel detail with sortable/filterable trait table
- `trait.html?panel=...&trait=...` — Trait detail with loci, candidate genes, eQTL support explanation

- `gene.html` — Gene search by symbol or Ensembl ID

## Local development

```bash
cd site && python3 -m http.server 8891
```

## Data

The frontend reads JSON from `data/`, a copy of the pipeline's `data/release/frontend/` output. To refresh it after a rebuild, run from the repository root:

```bash
python scripts/14_export_frontend_bundle.py
cp -R data/release/frontend/. site/data/
python scripts/15_build_gene_index.py
```

## Stack

Plain HTML, CSS, and vanilla JS. No build step, no framework dependencies.
