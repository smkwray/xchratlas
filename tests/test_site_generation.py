from pathlib import Path

import pandas as pd

from xatlas.site import PanelSpec, build_site, load_panel_release


def _write_release_dir(base: Path, rows: list[dict]) -> Path:
    base.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(base / "trait_scores.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"trait_id": "t1", "locus_id": "l1", "eqtl_supported": True},
            {"trait_id": "t1", "locus_id": "l2", "eqtl_supported": False},
            {"trait_id": "t2", "locus_id": "l3", "eqtl_supported": True},
        ]
    ).to_csv(base / "x_loci.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"trait_id": "t1", "locus_id": "l1", "gene_id": "g1"},
            {"trait_id": "t1", "locus_id": "l2", "gene_id": "g2"},
            {"trait_id": "t2", "locus_id": "l3", "gene_id": "g3"},
            {"trait_id": "t2", "locus_id": "l3", "gene_id": "g4"},
        ]
    ).to_csv(base / "x_gene_candidates.tsv", sep="\t", index=False)
    return base


def test_load_panel_release_summarizes_counts(tmp_path: Path) -> None:
    release_dir = _write_release_dir(
        tmp_path / "panel_alpha",
        [
            {
                "query_id": "trait_a",
                "description": "Trait A",
                "domain": "alpha",
                "n_loci": 3,
                "eqtl_supported": "True",
                "eqtl_total_hit_count": 10,
                "x_evidence_score": 100.0,
            },
            {
                "query_id": "trait_b",
                "description": "Trait B",
                "domain": "alpha",
                "n_loci": 0,
                "eqtl_supported": "False",
                "eqtl_total_hit_count": 0,
                "x_evidence_score": 20.0,
            },
        ],
    )
    panel = load_panel_release(PanelSpec("alpha", "Alpha", "Alpha panel", release_dir))

    assert panel["metrics"]["n_traits"] == 2
    assert panel["metrics"]["n_supported_traits"] == 1
    assert panel["metrics"]["n_zero_loci_traits"] == 1
    assert panel["metrics"]["n_loci"] == 3
    assert panel["metrics"]["n_gene_rows"] == 4
    assert panel["top_traits"].iloc[0]["query_id"] == "trait_a"


def test_build_site_writes_expected_pages(tmp_path: Path) -> None:
    release_a = _write_release_dir(
        tmp_path / "panel_a",
        [
            {
                "query_id": "trait_a",
                "description": "Trait A",
                "domain": "alpha",
                "n_loci": 2,
                "eqtl_supported": "True",
                "eqtl_total_hit_count": 7,
                "x_evidence_score": 88.0,
            }
        ],
    )
    release_b = _write_release_dir(
        tmp_path / "panel_b",
        [
            {
                "query_id": "trait_b",
                "description": "Trait B",
                "domain": "beta",
                "n_loci": 1,
                "eqtl_supported": "False",
                "eqtl_total_hit_count": 0,
                "x_evidence_score": 42.0,
            }
        ],
    )
    outdir = tmp_path / "site"
    build_site(
        [
            PanelSpec("panel-a", "Panel A", "First panel", release_a),
            PanelSpec("panel-b", "Panel B", "Second panel", release_b),
        ],
        outdir,
    )

    assert (outdir / "index.html").exists()
    assert (outdir / "methods.html").exists()
    assert (outdir / "panel-a.html").exists()
    assert (outdir / "panel-b.html").exists()
    assert (outdir / "style.css").exists()

    index_html = (outdir / "index.html").read_text(encoding="utf-8")
    assert "Panel A" in index_html
    assert "Panel B" in index_html

    panel_html = (outdir / "panel-a.html").read_text(encoding="utf-8")
    assert "Trait A" in panel_html
