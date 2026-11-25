#!/usr/bin/env python3
"""Build gene_index.json from existing trait detail JSONs in site/data/.

Scans every trait detail file across all panels and inverts the data:
  gene_name -> list of trait/locus appearances with eQTL evidence.

Output: site/data/gene_index.json
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "site" / "data"


def display_trait_name(trait: dict) -> str:
    if trait.get("query_id") == "operation_code" and trait.get("coding_description"):
        return str(trait["coding_description"])
    return str(trait.get("description") or trait.get("query_id") or trait.get("trait_id") or "Unlabeled trait")


def main() -> None:
    manifest_path = DATA_DIR / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    # gene_name -> { gene_id, biotype, hits: [...] }
    genes: dict[str, dict] = defaultdict(lambda: {"gene_id": None, "biotype": None, "hits": []})

    for panel in manifest["panels"]:
        panel_id = panel["panel_id"]
        panel_title = panel["title"]
        traits_dir = DATA_DIR / panel["files"]["trait_details_dir"]

        if not traits_dir.is_dir():
            print(f"[skip] {traits_dir} not found")
            continue

        for trait_path in sorted(traits_dir.glob("*.json")):
            detail = json.loads(trait_path.read_text(encoding="utf-8"))
            trait = detail["trait"]

            for locus in detail.get("loci", []):
                for gene in locus.get("candidate_genes", []):
                    gene_name = gene.get("gene_name")
                    if not gene_name:
                        continue

                    entry = genes[gene_name]
                    if entry["gene_id"] is None and gene.get("gene_id"):
                        entry["gene_id"] = gene["gene_id"]
                    if entry["biotype"] is None and gene.get("gene_biotype"):
                        entry["biotype"] = gene["gene_biotype"]

                    entry["hits"].append({
                        "trait_id": trait["trait_id"],
                        "trait_desc": display_trait_name(trait),
                        "panel_id": panel_id,
                        "panel_title": panel_title,
                        "domain": trait.get("domain", ""),
                        "locus_region": locus.get("x_region", ""),
                        "lead_rsid": locus.get("lead_rsid", ""),
                        "lead_pos": locus.get("lead_pos"),
                        "mapping_relation": gene.get("mapping_relation", ""),
                        "eqtl_gene_role": gene.get("eqtl_gene_role", ""),
                        "eqtl_lookup_hit": gene.get("eqtl_lookup_hit", False),
                        "eqtl_supported": gene.get("eqtl_supported", False),
                        "best_eqtl_pvalue": gene.get("best_eqtl_pvalue"),
                        "eqtl_study_count": gene.get("eqtl_study_count", 0),
                        "eqtl_dataset_count": gene.get("eqtl_dataset_count", 0),
                        "eqtl_lookup_hit_count": gene.get("eqtl_lookup_hit_count", 0),
                    })

    # Deduplicate hits (same gene can appear in overlapping panels for same trait/locus)
    for gene_name, entry in genes.items():
        seen = set()
        unique_hits = []
        for hit in entry["hits"]:
            key = (hit["trait_id"], hit["panel_id"], hit["lead_rsid"])
            if key not in seen:
                seen.add(key)
                unique_hits.append(hit)
        # Sort: eQTL-supported first, then by best p-value
        unique_hits.sort(key=lambda h: (
            not h["eqtl_supported"],
            not h["eqtl_lookup_hit"],
            h["best_eqtl_pvalue"] if h["best_eqtl_pvalue"] is not None else 1.0,
        ))
        entry["hits"] = unique_hits

    # Sort genes: most eQTL-supported hits first, then by total hits
    sorted_genes = dict(sorted(
        genes.items(),
        key=lambda kv: (
            -sum(1 for h in kv[1]["hits"] if h["eqtl_supported"]),
            -sum(1 for h in kv[1]["hits"] if h["eqtl_lookup_hit"]),
            -len(kv[1]["hits"]),
            kv[0],
        ),
    ))

    outpath = DATA_DIR / "gene_index.json"
    outpath.write_text(
        json.dumps(sorted_genes, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )

    total_genes = len(sorted_genes)
    with_eqtl = sum(1 for g in sorted_genes.values() if any(h["eqtl_lookup_hit"] for h in g["hits"]))
    with_support = sum(1 for g in sorted_genes.values() if any(h["eqtl_supported"] for h in g["hits"]))
    print(f"[done] gene_index.json: {total_genes} genes, {with_eqtl} with eQTL hits, {with_support} with strict support")


if __name__ == "__main__":
    main()
