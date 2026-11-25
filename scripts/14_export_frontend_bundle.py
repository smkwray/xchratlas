#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

from xatlas.frontend_bundle import default_panel_specs, export_frontend_bundle
from xatlas.settings import RELEASE_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="Export frontend-ready JSON bundles from chrXatlas release tables.")
    parser.add_argument("--release-root", default=str(RELEASE_DIR))
    parser.add_argument("--outdir", default=str(RELEASE_DIR / "frontend"))
    parser.add_argument("--expanded-release-dir", default=None)
    parser.add_argument("--mind-release-dir", default=None)
    parser.add_argument("--biochem-release-dir", default=None)
    parser.add_argument("--independent-release-dir", default=None)
    args = parser.parse_args()

    specs = default_panel_specs(args.release_root)
    override_map = {
        "panel_expanded": args.expanded_release_dir,
        "panel_mind_risk_core": args.mind_release_dir,
        "panel_blood_biochemistry_core": args.biochem_release_dir,
        "panel_independent_set": args.independent_release_dir,
    }
    specs = [
        spec if override_map.get(spec.panel_id) is None else type(spec)(
            panel_id=spec.panel_id,
            slug=spec.slug,
            title=spec.title,
            description=spec.description,
            featured=spec.featured,
            release_dir=Path(override_map[spec.panel_id]),
        )
        for spec in specs
    ]

    manifest = export_frontend_bundle(specs, args.outdir)
    print(f"[done] frontend bundle exported to {args.outdir}/")
    print(f"[done] panels: {len(manifest['panels'])}")


if __name__ == "__main__":
    main()
