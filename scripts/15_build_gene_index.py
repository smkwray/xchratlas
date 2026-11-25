#!/usr/bin/env python3
"""Build a versionless Ensembl gene index from the listed panel traits."""
from __future__ import annotations

import json
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

    genes = {}
    seen = set()
    for panel in manifest["panels"]:
        panel_id = panel["panel_id"]
        panel_title = panel["title"]
        traits_dir = DATA_DIR / panel["files"]["trait_details_dir"]

        rows = json.loads((DATA_DIR / panel["files"]["traits"]).read_text(encoding="utf-8"))
        for row in rows:
            trait_path = traits_dir / (row["trait_id"] + ".json")
            detail = json.loads(trait_path.read_text(encoding="utf-8"))
            trait = detail["trait"]
            if trait["trait_id"] != row["trait_id"]:
                raise ValueError(f"Unexpected trait in {trait_path}")
            for locus in detail["loci"]:
                for gene in locus["candidate_genes"]:
                    gene_id = gene["gene_id"].split(".")[0]
                    key = (gene_id, trait["trait_id"], panel_id, locus["locus_id"])
                    if key in seen:
                        continue
                    seen.add(key)
                    entry = genes.setdefault(gene_id, {
                        "gene_id": gene_id, "gene_name": gene["gene_name"],
                        "biotype": gene.get("gene_biotype"), "hits": [],
                    })
                    entry["hits"].append({
                        "trait_id": trait["trait_id"],
                        "trait_desc": display_trait_name(trait),
                        "panel_id": panel_id,
                        "panel_title": panel_title,
                        "domain": trait.get("domain", ""),
                        "locus_id": locus["locus_id"],
                        "eqtl_lookup_status": locus.get("eqtl_lookup_status"),
                        **{k: locus[k] for k in ("lead_selection", "lead_neglog10_is_lower_bound",
                                                  "lead_neglog10_pvalue", "lead_rsid", "lead_pos") if locus.get("lead_neglog10_is_lower_bound")},
                        "mapping_relation": gene.get("mapping_relation", ""),
                        "eqtl_gene_role": gene.get("eqtl_gene_role", ""),
                        "eqtl_lookup_hit": gene.get("eqtl_lookup_hit", False),
                        "eqtl_supported": gene.get("eqtl_supported", False),
                        "best_eqtl_pvalue": gene.get("best_eqtl_pvalue"),
                    })

    for entry in genes.values():
        hits = entry["hits"]
        entry["returned_trait_count"] = len({h["trait_id"] for h in hits if h["eqtl_lookup_hit"]})
        entry["threshold_trait_count"] = len({h["trait_id"] for h in hits if h["eqtl_supported"]})
        entry["candidate_threshold_trait_count"] = len({
            h["trait_id"] for h in hits if h["eqtl_supported"] and h["eqtl_gene_role"] != "followup"
        })
    outpath = DATA_DIR / "gene_index.json"
    outpath.write_text(json.dumps(dict(sorted(genes.items())), separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"[done] {len(genes)} genes")


if __name__ == "__main__":
    main()
