#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

from xatlas.settings import RELEASE_DIR
from xatlas.site import PanelSpec, build_site


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a small static chrXatlas site from release tables.")
    parser.add_argument("--expanded-release-dir", default=str(RELEASE_DIR / "panel_expanded"))
    parser.add_argument("--mind-release-dir", default=str(RELEASE_DIR / "panel_mind_risk_core"))
    parser.add_argument("--biochem-release-dir", default=str(RELEASE_DIR / "panel_blood_biochemistry_core"))
    parser.add_argument("--independent-release-dir", default=str(RELEASE_DIR / "panel_independent_set"))
    parser.add_argument("--outdir", default=str(RELEASE_DIR / "site"))
    args = parser.parse_args()

    candidates = [
        PanelSpec(
            slug="broad-atlas",
            title="Broad Atlas",
            blurb="Curated all-purpose release with broad trait coverage and low zero-locus drift.",
            release_dir=Path(args.expanded_release_dir),
        ),
        PanelSpec(
            slug="mind-and-risk",
            title="Mind And Risk",
            blurb="Focused behavior and cognition companion panel built around the strongest mind-adjacent chrX traits.",
            release_dir=Path(args.mind_release_dir),
        ),
        PanelSpec(
            slug="biochemistry-deep-dive",
            title="Biochemistry Deep Dive",
            blurb="High-yield blood biochemistry panel carved out of the independent-set discovery pool.",
            release_dir=Path(args.biochem_release_dir),
        ),
        PanelSpec(
            slug="discovery-pool",
            title="Discovery Pool",
            blurb="Full Pan-UKB max independent-set build kept as a broad curation source rather than a polished release tier.",
            release_dir=Path(args.independent_release_dir),
            featured=False,
        ),
    ]

    panel_specs = [spec for spec in candidates if (spec.release_dir / "trait_scores.tsv").exists()]
    if not panel_specs:
        raise SystemExit("No release tables found for site generation.")

    outdir = build_site(panel_specs, args.outdir)
    print(f"[done] static site written to {outdir}/")


if __name__ == "__main__":
    main()
