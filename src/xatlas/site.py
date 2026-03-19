from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Iterable

import pandas as pd

from xatlas.io import ensure_dir

STYLE_CSS = """
:root {
  --bg: #f6f1e8;
  --card: #fffaf2;
  --ink: #1d1b18;
  --muted: #6d655b;
  --line: #d8cbb6;
  --accent: #9c4f2f;
  --accent-soft: #efe0d2;
  --good: #236a48;
  --warn: #8f5b19;
  --mono: "SFMono-Regular", "Menlo", "Monaco", monospace;
  --serif: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", Georgia, serif;
}

* { box-sizing: border-box; }
body {
  margin: 0;
  color: var(--ink);
  background:
    radial-gradient(circle at top left, rgba(156, 79, 47, 0.10), transparent 38%),
    radial-gradient(circle at top right, rgba(35, 106, 72, 0.09), transparent 30%),
    linear-gradient(180deg, #f7f1e7 0%, #f2eadf 100%);
  font-family: var(--serif);
  line-height: 1.55;
}

a { color: var(--accent); }
a:hover { color: #6f381f; }

.shell {
  max-width: 1180px;
  margin: 0 auto;
  padding: 28px 20px 56px;
}

.hero {
  padding: 28px;
  border: 1px solid var(--line);
  border-radius: 22px;
  background: linear-gradient(145deg, rgba(255,250,242,0.98), rgba(244,233,220,0.92));
  box-shadow: 0 12px 40px rgba(86, 63, 30, 0.08);
}

.eyebrow {
  font-family: var(--mono);
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.14em;
  color: var(--muted);
}

h1, h2, h3 {
  margin: 0 0 12px;
  line-height: 1.1;
  font-weight: 700;
}

h1 { font-size: clamp(2.3rem, 4vw, 4.4rem); max-width: 10ch; }
h2 { font-size: 1.7rem; margin-top: 30px; }
h3 { font-size: 1.15rem; }

p.lead {
  max-width: 68ch;
  font-size: 1.08rem;
  color: var(--muted);
}

.nav {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin: 18px 0 0;
}

.nav a,
.pill {
  display: inline-block;
  padding: 9px 13px;
  border-radius: 999px;
  border: 1px solid var(--line);
  background: rgba(255,255,255,0.62);
  text-decoration: none;
  color: var(--ink);
  font-family: var(--mono);
  font-size: 0.88rem;
}

.grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 14px;
  margin-top: 20px;
}

.card {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 18px;
  padding: 18px;
  box-shadow: 0 8px 28px rgba(70, 52, 28, 0.05);
}

.metric {
  font-size: 2rem;
  font-weight: 700;
  margin: 6px 0 4px;
}

.metric-label,
.muted {
  color: var(--muted);
}

.table-wrap {
  margin-top: 18px;
  overflow-x: auto;
  background: rgba(255,255,255,0.72);
  border: 1px solid var(--line);
  border-radius: 18px;
}

table {
  width: 100%;
  border-collapse: collapse;
  min-width: 760px;
}

th, td {
  padding: 12px 14px;
  border-bottom: 1px solid rgba(216, 203, 182, 0.85);
  vertical-align: top;
}

th {
  text-align: left;
  font-family: var(--mono);
  font-size: 0.82rem;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--muted);
}

tbody tr:hover { background: rgba(239, 224, 210, 0.28); }

.ok { color: var(--good); font-weight: 700; }
.warn { color: var(--warn); font-weight: 700; }

.two-up {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 18px;
}

.footer {
  margin-top: 32px;
  padding-top: 18px;
  border-top: 1px solid var(--line);
  color: var(--muted);
  font-size: 0.95rem;
}

code {
  font-family: var(--mono);
  font-size: 0.95em;
  background: rgba(239, 224, 210, 0.45);
  padding: 0.1em 0.35em;
  border-radius: 0.3em;
}
""".strip()


@dataclass(frozen=True)
class PanelSpec:
    slug: str
    title: str
    blurb: str
    release_dir: Path
    featured: bool = True


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _fmt_int(value: object) -> str:
    try:
        number = int(float(value))
    except (TypeError, ValueError):
        return ""
    return f"{number:,}"


def _relative_release_link(site_outdir: Path, release_dir: Path, filename: str) -> str:
    return Path("../") / release_dir.name / filename


def _render_nav(panel_specs: Iterable[PanelSpec], *, current: str) -> str:
    links = ['<a href="index.html">Overview</a>']
    specs = list(panel_specs)
    featured_specs = [spec for spec in specs if spec.featured]
    current_spec = next((spec for spec in specs if spec.slug == current), None)
    if current_spec is not None and not current_spec.featured:
        featured_specs.append(current_spec)

    for spec in featured_specs:
        href = f"{spec.slug}.html"
        if current == spec.slug:
            links.append(f'<span class="pill">{escape(spec.title)}</span>')
        else:
            links.append(f'<a href="{href}">{escape(spec.title)}</a>')
    if current == "methods":
        links.append('<span class="pill">Methods</span>')
    else:
        links.append('<a href="methods.html">Methods</a>')
    return f'<div class="nav">{"".join(links)}</div>'


def _render_table(df: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    if df.empty:
        return '<p class="muted">No rows.</p>'
    header = "".join(f"<th>{escape(label)}</th>" for _, label in columns)
    body_rows = []
    for _, row in df.iterrows():
        cells = []
        for key, _ in columns:
            value = row.get(key, "")
            cells.append(f"<td>{escape(str(value if pd.notna(value) else ''))}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")
    return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{"".join(body_rows)}</tbody></table></div>'


def load_panel_release(spec: PanelSpec) -> dict[str, object]:
    trait_scores = pd.read_csv(spec.release_dir / "trait_scores.tsv", sep="\t")
    loci = pd.read_csv(spec.release_dir / "x_loci.tsv", sep="\t")
    gene_candidates = pd.read_csv(spec.release_dir / "x_gene_candidates.tsv", sep="\t")

    trait_scores = trait_scores.copy()
    trait_scores["eqtl_lookup_hit_bool"] = trait_scores.get("eqtl_lookup_hit", False).map(_truthy)
    trait_scores["eqtl_supported_bool"] = trait_scores["eqtl_supported"].map(_truthy)
    trait_scores["n_loci"] = pd.to_numeric(trait_scores.get("n_loci"), errors="coerce").fillna(0).astype(int)
    trait_scores["eqtl_total_hit_count"] = pd.to_numeric(trait_scores.get("eqtl_total_hit_count"), errors="coerce").fillna(0).astype(int)
    trait_scores["x_evidence_score"] = pd.to_numeric(trait_scores.get("x_evidence_score"), errors="coerce").fillna(0.0)

    top_traits = (
        trait_scores[
            ["query_id", "description", "domain", "n_loci", "eqtl_lookup_hit_bool", "eqtl_supported_bool", "eqtl_total_hit_count", "x_evidence_score"]
        ]
        .sort_values(["x_evidence_score", "n_loci", "eqtl_total_hit_count"], ascending=[False, False, False])
        .head(20)
        .copy()
    )
    top_traits["eqtl_lookup_hit_bool"] = top_traits["eqtl_lookup_hit_bool"].map(lambda value: "yes" if value else "no")
    top_traits["eqtl_supported_bool"] = top_traits["eqtl_supported_bool"].map(lambda value: "yes" if value else "no")
    top_traits["x_evidence_score"] = top_traits["x_evidence_score"].map(lambda value: f"{value:.2f}")

    unsupported = (
        trait_scores.loc[~trait_scores["eqtl_supported_bool"], ["query_id", "description", "domain", "n_loci", "x_evidence_score"]]
        .sort_values(["x_evidence_score", "n_loci"], ascending=[False, False])
        .head(20)
        .copy()
    )
    unsupported["x_evidence_score"] = unsupported["x_evidence_score"].map(lambda value: f"{value:.2f}")

    domains = (
        trait_scores.groupby("domain", dropna=False)
        .agg(
            n_traits=("query_id", "count"),
            lookup_hit_traits=("eqtl_lookup_hit_bool", "sum"),
            supported_traits=("eqtl_supported_bool", "sum"),
            total_loci=("n_loci", "sum"),
        )
        .reset_index()
        .sort_values(["supported_traits", "total_loci", "n_traits"], ascending=[False, False, False])
        .head(16)
    )

    return {
        "spec": spec,
        "trait_scores": trait_scores,
        "top_traits": top_traits,
        "unsupported": unsupported,
        "domains": domains,
        "metrics": {
            "n_traits": len(trait_scores),
            "n_lookup_hit_traits": int(trait_scores["eqtl_lookup_hit_bool"].sum()),
            "n_supported_traits": int(trait_scores["eqtl_supported_bool"].sum()),
            "n_zero_loci_traits": int((trait_scores["n_loci"] == 0).sum()),
            "n_loci": len(loci),
            "n_gene_rows": len(gene_candidates),
        },
    }


def _panel_summary_cards(panels: list[dict[str, object]], site_outdir: Path) -> str:
    cards = []
    for panel in panels:
        spec = panel["spec"]
        if not spec.featured:
            continue
        metrics = panel["metrics"]
        cards.append(
            "\n".join(
                [
                    '<div class="card">',
                    f"<h3>{escape(spec.title)}</h3>",
                    f"<p class=\"muted\">{escape(spec.blurb)}</p>",
                    f"<div class=\"metric\">{_fmt_int(metrics['n_lookup_hit_traits'])}/{_fmt_int(metrics['n_traits'])}</div>",
                    '<div class="metric-label">eQTL lookup-hit traits</div>',
                    f"<div class=\"metric\">{_fmt_int(metrics['n_supported_traits'])}/{_fmt_int(metrics['n_traits'])}</div>",
                    '<div class="metric-label">eQTL-supported traits</div>',
                    '<div class="grid">',
                    f'<div><div class="metric">{_fmt_int(metrics["n_loci"])}</div><div class="metric-label">chrX loci</div></div>',
                    f'<div><div class="metric">{_fmt_int(metrics["n_zero_loci_traits"])}</div><div class="metric-label">Zero-locus traits</div></div>',
                    f'<div><div class="metric">{_fmt_int(metrics["n_gene_rows"])}</div><div class="metric-label">Locus-gene rows</div></div>',
                    f'<div><div class="metric"><a href="{escape(spec.slug)}.html">open</a></div><div class="metric-label">Panel page</div></div>',
                    "</div>",
                    "</div>",
                ]
            )
        )
    return '<div class="grid">' + "".join(cards) + "</div>"


def _panel_comparison_table(panels: list[dict[str, object]]) -> str:
    rows = []
    for panel in panels:
        spec = panel["spec"]
        if not spec.featured:
            continue
        metrics = panel["metrics"]
        rows.append(
            {
                "panel": spec.title,
                "focus": spec.blurb,
                "traits": _fmt_int(metrics["n_traits"]),
                "lookup_hit": f'{_fmt_int(metrics["n_lookup_hit_traits"])}/{_fmt_int(metrics["n_traits"])}',
                "supported": f'{_fmt_int(metrics["n_supported_traits"])}/{_fmt_int(metrics["n_traits"])}',
                "loci": _fmt_int(metrics["n_loci"]),
                "zero_loci": _fmt_int(metrics["n_zero_loci_traits"]),
            }
        )
    return _render_table(pd.DataFrame(rows), [("panel", "Panel"), ("focus", "Focus"), ("traits", "Traits"), ("lookup_hit", "eQTL lookup-hit"), ("supported", "eQTL-supported"), ("loci", "Loci"), ("zero_loci", "Zero-locus")])


def _render_discovery_block(panels: list[dict[str, object]]) -> str:
    secondary = [panel for panel in panels if not panel["spec"].featured]
    if not secondary:
        return ""

    cards = []
    for panel in secondary:
        spec = panel["spec"]
        metrics = panel["metrics"]
        cards.append(
            "\n".join(
                [
                    '<div class="card">',
                    f"<h3>{escape(spec.title)}</h3>",
                    f"<p class=\"muted\">{escape(spec.blurb)}</p>",
                    f"<p><span class=\"warn\">{_fmt_int(metrics['n_zero_loci_traits'])}</span> zero-locus traits and <span class=\"warn\">{_fmt_int(metrics['n_traits'] - metrics['n_supported_traits'])}</span> traits without strict eQTL support. Use this page as a scouting surface, not as the main story.</p>",
                    f'<p><a href="{escape(spec.slug)}.html">Open the discovery page</a></p>',
                    "</div>",
                ]
            )
        )
    return "<h2>Discovery Pool</h2><p class=\"lead\">These broader pages are useful for curation, but they are intentionally kept secondary so the public-facing story stays anchored on the cleaner curated panels.</p><div class=\"grid\">" + "".join(cards) + "</div>"


def _render_panel_page(panel: dict[str, object], panel_specs: list[PanelSpec], site_outdir: Path) -> str:
    spec = panel["spec"]
    metrics = panel["metrics"]
    trait_link = _relative_release_link(site_outdir, spec.release_dir, "trait_scores.tsv")
    loci_link = _relative_release_link(site_outdir, spec.release_dir, "x_loci.tsv")
    genes_link = _relative_release_link(site_outdir, spec.release_dir, "x_gene_candidates.tsv")
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
            f"<title>{escape(spec.title)} | chrXatlas</title>",
            '<link rel="stylesheet" href="style.css"></head><body><div class="shell">',
            '<section class="hero">',
            '<div class="eyebrow">chrXatlas panel</div>',
            f"<h1>{escape(spec.title)}</h1>",
            f'<p class="lead">{escape(spec.blurb)}</p>',
            _render_nav(panel_specs, current=spec.slug),
            "</section>",
            '<section class="grid">',
            f'<div class="card"><div class="metric">{_fmt_int(metrics["n_lookup_hit_traits"])}/{_fmt_int(metrics["n_traits"])}</div><div class="metric-label">eQTL lookup-hit traits</div></div>',
            f'<div class="card"><div class="metric">{_fmt_int(metrics["n_supported_traits"])}/{_fmt_int(metrics["n_traits"])}</div><div class="metric-label">eQTL-supported traits</div></div>',
            f'<div class="card"><div class="metric">{_fmt_int(metrics["n_loci"])}</div><div class="metric-label">chrX loci</div></div>',
            f'<div class="card"><div class="metric">{_fmt_int(metrics["n_zero_loci_traits"])}</div><div class="metric-label">Zero-locus traits</div></div>',
            "</section>",
            '<section class="two-up">',
            '<div class="card">',
            "<h2>Underlying tables</h2>",
            f'<p><a href="{escape(str(trait_link))}">trait_scores.tsv</a><br><a href="{escape(str(loci_link))}">x_loci.tsv</a><br><a href="{escape(str(genes_link))}">x_gene_candidates.tsv</a></p>',
            "</div>",
            '<div class="card">',
            "<h2>Plain English</h2>",
            "<p>This panel is a curated slice of the chrX atlas. Each trait was scanned for strong chromosome X associations, nearby hits were merged into loci, and loci were then linked to nearby genes plus outside gene-activity evidence when safely available.</p>",
            "</div>",
            "</section>",
            "<h2>Top Traits</h2>",
            _render_table(
                panel["top_traits"],
                [
                    ("query_id", "Trait"),
                    ("description", "Description"),
                    ("domain", "Domain"),
                    ("n_loci", "Loci"),
                    ("eqtl_lookup_hit_bool", "eQTL lookup-hit"),
                    ("eqtl_supported_bool", "eQTL-supported"),
                    ("eqtl_total_hit_count", "eQTL Hits"),
                    ("x_evidence_score", "Score"),
                ],
            ),
            "<h2>Domain Summary</h2>",
            _render_table(panel["domains"], [("domain", "Domain"), ("n_traits", "Traits"), ("lookup_hit_traits", "eQTL lookup-hit"), ("supported_traits", "eQTL-supported"), ("total_loci", "Loci")]),
            "<h2>Unsupported Or Weak Traits</h2>",
            _render_table(panel["unsupported"], [("query_id", "Trait"), ("description", "Description"), ("domain", "Domain"), ("n_loci", "Loci"), ("x_evidence_score", "Score")]),
            '<div class="footer">Generated from release tables. eQTL lookup-hit means returned associations were observed in the prioritized lookup path; eQTL-supported means a candidate gene met the strict support rule in the current build.</div>',
            "</div></body></html>",
        ]
    )


def _render_methods_page(panel_specs: list[PanelSpec]) -> str:
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
            "<title>Methods | chrXatlas</title>",
            '<link rel="stylesheet" href="style.css"></head><body><div class="shell">',
            '<section class="hero">',
            '<div class="eyebrow">chrXatlas methods</div>',
            "<h1>How chrXatlas works</h1>",
            '<p class="lead">In plain English: chrXatlas looks for strong statistical signals on chromosome X, groups nearby hits into loci, then tries to connect those loci to likely genes and outside gene-activity evidence.</p>',
            _render_nav(panel_specs, current="methods"),
            "</section>",
            '<section class="two-up">',
            '<div class="card"><h2>Signal finding</h2><p>For each public trait, the pipeline uses Pan-UKB summary statistics and scans only chromosome X. Strongly associated variants are kept, nearby hits are merged into one locus, and the lead variant from each locus is recorded.</p></div>',
            '<div class="card"><h2>Gene follow-up</h2><p>Each locus is mapped to nearby genes in GRCh37 coordinates. eQTL follow-up is added through rsID-first lookups into the eQTL Catalogue, with GRCh37 Ensembl variant recoder rescue for some missing rsIDs.</p></div>',
            "</section>",
            '<section class="two-up">',
            '<div class="card"><h2>Robustness checks</h2><p>The build uses a strict genome-wide significance threshold, keeps PAR and nonPAR regions separate, infers the p-value schema from the actual Pan-UKB file, and avoids unsafe direct GRCh37-to-GRCh38 region lookups by default.</p></div>',
            '<div class="card"><h2>How to read eQTL labels</h2><p>An eQTL lookup-hit means returned associations were observed in the prioritized datasets. eQTL-supported is narrower: at least one candidate gene passed the strict support rule in the current build. Neither label is a causal claim.</p></div>',
            "</section>",
            '<div class="footer">The curated panels are intended for presentation. The broad independent-set panel is best treated as a discovery pool for further curation.</div>',
            "</div></body></html>",
        ]
    )


def build_site(panel_specs: list[PanelSpec], outdir: str | Path) -> Path:
    site_outdir = ensure_dir(outdir)
    panels = [load_panel_release(spec) for spec in panel_specs]

    (site_outdir / "style.css").write_text(STYLE_CSS + "\n", encoding="utf-8")

    index_html = "\n".join(
        [
            "<!doctype html>",
            '<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
            "<title>chrXatlas</title>",
            '<link rel="stylesheet" href="style.css"></head><body><div class="shell">',
            '<section class="hero">',
            '<div class="eyebrow">chrXatlas</div>',
            "<h1>Chromosome X signals, made legible.</h1>",
            '<p class="lead">This presentation layer turns the release tables into a simple, browsable summary. The curated panels are the main story; the broad independent-set build is included as a discovery pool, not as the cleanest public-facing tier.</p>',
            _render_nav(panel_specs, current="index"),
            "</section>",
            _panel_summary_cards(panels, site_outdir),
            "<h2>Panel Comparison</h2>",
            _panel_comparison_table(panels),
            _render_discovery_block(panels),
            '<section class="two-up">',
            '<div class="card"><h2>What counts as a chrX signal?</h2><p>A chrX signal is a spot on chromosome X where the public genetics data shows a strong statistical link to a trait. Nearby hits are grouped into one locus so a single region is not counted over and over.</p></div>',
            '<div class="card"><h2>How to read eQTL evidence</h2><p>eQTL lookup-hit means the rsID-based follow-up returned associations in the prioritized datasets. eQTL-supported is stricter and only counts candidate-gene evidence that passes the current support rule. Both are discovery signals, not proof of direct causation.</p></div>',
            "</section>",
            '<div class="footer">Static HTML generated from release tables. Open the panel pages above for top traits, domain summaries, and direct links to the TSV outputs.</div>',
            "</div></body></html>",
        ]
    )
    (site_outdir / "index.html").write_text(index_html + "\n", encoding="utf-8")
    (site_outdir / "methods.html").write_text(_render_methods_page(panel_specs) + "\n", encoding="utf-8")

    for panel in panels:
        spec = panel["spec"]
        (site_outdir / f"{spec.slug}.html").write_text(_render_panel_page(panel, panel_specs, site_outdir) + "\n", encoding="utf-8")

    return site_outdir
