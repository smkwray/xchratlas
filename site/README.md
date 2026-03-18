# chrXatlas frontend

Interactive static site for the Chromosome X Evidence Atlas.

## Pages

- `index.html` — Landing page with hero, metrics, curated panels, domain composition, top traits, search
- `panel.html?id=...` — Panel detail with sortable/filterable trait table
- `trait.html?panel=...&trait=...` — Trait detail with loci, candidate genes, eQTL support explanation

## Local development

```bash
cd site && python3 -m http.server 8891
```

## Data

The frontend reads JSON from `data/`, which mirrors the backend's `release/frontend/` output. To refresh after a pipeline rebuild:

```bash
cp -R "$PROJ_SHARED_DATA_ROOT/release/frontend/." site/data/
```

## Stack

Plain HTML, CSS, and vanilla JS. No build step, no framework dependencies.
