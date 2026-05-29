#!/usr/bin/env python3
"""Run the perturbagen pipeline end to end for the current targets.yml."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "targets.yml"
SCRIPTS = ROOT / "scripts"
RUNS_DIR = ROOT / "runs"

ARCHIVE_FILES = [
    ROOT / "data" / "processed" / "perturbagen_source_records.csv",
    ROOT / "data" / "processed" / "perturbagens_merged.csv",
    ROOT / "data" / "processed" / "perturbagen_master_table.csv",
    ROOT / "output" / "perturbagen_master_table.csv",
    ROOT / "output" / "stem_cell_screen_ready_list.csv",
]

ARCHIVE_RAW_DIRS = [
    ROOT / "data" / "raw" / "chembl",
    ROOT / "data" / "raw" / "iuphar",
]


def run_step(label: str, command: list[str]) -> None:
    print(f"\n== {label} ==", flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def slugify(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "unnamed_run"


def load_config() -> dict:
    if not CONFIG.exists():
        raise FileNotFoundError(f"Missing config file: {CONFIG}")
    with open(CONFIG, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def infer_run_name(config: dict) -> str:
    explicit = config.get("run_name") or config.get("cell_type") or config.get("name")
    if explicit:
        return slugify(str(explicit))

    pathways = sorted((config.get("pathways") or {}).keys())
    if pathways:
        return slugify("_".join(pathways[:3]))
    return "targets"


def copy_file(src: Path, dest: Path) -> None:
    if not src.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def copy_tree(src: Path, dest: Path) -> None:
    if not src.exists():
        return
    if dest.exists():
        shutil.rmtree(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)


def archive_run(run_name: str, archive_raw: bool = True) -> Path:
    run_dir = RUNS_DIR / run_name
    copy_file(CONFIG, run_dir / "config" / "targets.yml")

    for src in ARCHIVE_FILES:
        relative = src.relative_to(ROOT)
        copy_file(src, run_dir / relative)

    if archive_raw:
        for src in ARCHIVE_RAW_DIRS:
            relative = src.relative_to(ROOT)
            copy_tree(src, run_dir / relative)

    summary = run_dir / "README.txt"
    summary.write_text(
        "\n".join(
            [
                f"run_name: {run_name}",
                f"config: {run_dir / 'config' / 'targets.yml'}",
                f"master_table: {run_dir / 'output' / 'perturbagen_master_table.csv'}",
                f"screen_ready_list: {run_dir / 'output' / 'stem_cell_screen_ready_list.csv'}",
                "",
                "This folder is a snapshot from scripts/run_pipeline.py.",
                "The top-level output/ directory always contains the latest run.",
            ]
        ),
        encoding="utf-8",
    )
    return run_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=["all", "lincs", "chembl", "iuphar", "selleck_l2100", "vendors"],
        default="all",
        help="Raw source to refresh before merge. Use all after changing targets.yml.",
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Reuse existing raw/normalized source files and only rebuild merged/scored outputs.",
    )
    parser.add_argument(
        "--skip-score",
        action="store_true",
        help="Stop after rebuilding the merged master table.",
    )
    parser.add_argument(
        "--run-name",
        help="Archive results under runs/<run-name>. Defaults to run_name/cell_type in config/targets.yml.",
    )
    parser.add_argument(
        "--no-archive-raw",
        action="store_true",
        help="Do not copy target-dependent raw ChEMBL/IUPHAR files into runs/<run-name>/data/raw.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    python = sys.executable
    config = load_config()
    run_name = slugify(args.run_name) if args.run_name else infer_run_name(config)

    if not args.skip_fetch:
        run_step(
            "Fetch and normalize source data",
            [python, str(SCRIPTS / "01_fetch_sources.py"), "--source", args.source],
        )

    run_step("Merge and deduplicate into master schema", [python, str(SCRIPTS / "02_normalize_merge.py")])

    if not args.skip_score:
        run_step("Score and sort screen-ready list", [python, str(SCRIPTS / "04_score_for_ipsc.py")])

    run_dir = archive_run(run_name, archive_raw=not args.no_archive_raw)

    print("\nPipeline complete.", flush=True)
    print(f"Run archive: {run_dir}", flush=True)
    print(f"Master table: {ROOT / 'output' / 'perturbagen_master_table.csv'}", flush=True)
    print(f"Screen-ready list: {ROOT / 'output' / 'stem_cell_screen_ready_list.csv'}", flush=True)


if __name__ == "__main__":
    main()
