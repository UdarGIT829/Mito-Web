#!/usr/bin/env python3
"""Register the project's EV and lupusN sources in the analysis catalog."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from mito_viewer.catalog import (  # noqa: E402
    CatalogRepository,
    register_current_source_cohorts,
)


DEFAULT_CATALOG_PATH = PROJECT_ROOT / "analysis_catalog.sqlite"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Register EV_Study.sqlite and lupusN.sqlite as read-only "
            "source cohorts."
        )
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DEFAULT_CATALOG_PATH,
        help=f"Writable catalog path. Default: {DEFAULT_CATALOG_PATH}",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help=f"Directory containing the source databases. Default: {PROJECT_ROOT}",
    )
    args = parser.parse_args()

    with CatalogRepository(args.catalog) as catalog:
        results = register_current_source_cohorts(
            catalog,
            args.project_root,
        )

    for result in results:
        status = "registered" if result.created else "already registered"
        print(
            f"{result.cohort.name}: {status}; "
            f"{len(result.samples)} observed samples; "
            f"{result.cohort.source_database_fingerprint}"
        )
    print(f"Catalog: {args.catalog.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
