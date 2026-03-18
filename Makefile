PYTHON ?= python

setup:
	$(PYTHON) -m venv .venv
	. .venv/bin/activate && pip install -r requirements.txt && pip install -e .

fetch-manifests:
	$(PYTHON) scripts/00_fetch_panukb_manifests.py
	$(PYTHON) scripts/04_prepare_eqtl_index.py

select-traits:
	$(PYTHON) scripts/01_select_seed_traits.py --panel config/panel_small.csv

extract-x:
	$(PYTHON) scripts/02_extract_panukb_x.py

genes:
	$(PYTHON) scripts/03_build_x_gene_catalog.py

loci:
	$(PYTHON) scripts/05_call_x_loci.py

eqtl:
	$(PYTHON) scripts/06_fetch_eqtl_region_hits.py

map-genes:
	$(PYTHON) scripts/07_map_loci_to_genes.py

release:
	$(PYTHON) scripts/08_export_release_tables.py

test:
	$(PYTHON) -m pytest
