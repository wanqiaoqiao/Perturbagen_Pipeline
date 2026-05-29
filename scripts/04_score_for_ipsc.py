#!/usr/bin/env python3
"""Score merged perturbagens for stem-cell and iPSC screening priority."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "targets.yml"
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "output"
PROCESSED_MASTER_TABLE = PROCESSED_DIR / "perturbagen_master_table.csv"
OUTPUT_MASTER_TABLE = OUTPUT_DIR / "perturbagen_master_table.csv"
SCREEN_READY_LIST = OUTPUT_DIR / "stem_cell_screen_ready_list.csv"

DEV_KEYWORDS = [
    "Wnt",
    "GSK3",
    "TGF",
    "Activin",
    "BMP",
    "FGF",
    "Hedgehog",
    "SMO",
    "Notch",
    "retinoic",
    "JAK",
    "STAT",
    "HDAC",
    "DNMT",
    "ROCK",
]

def text_value(row: pd.Series, *names: str) -> str:
    values = []
    for name in names:
        value = row.get(name)
        if value is None or pd.isna(value):
            continue
        text = str(value).strip()
        if text and text.lower() not in {"nan", "none", "null"}:
            values.append(text)
    return " ".join(values)


def load_scoring_keywords() -> list[str]:
    keywords = set()
    if not CONFIG.exists():
        return sorted(DEV_KEYWORDS, key=str.lower)

    with open(CONFIG, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    for pathway, payload in (config.get("pathways") or {}).items():
        keywords.add(str(pathway).replace("_", " "))
        keywords.add(str(pathway))
        for gene in payload.get("genes", []) or []:
            keywords.add(str(gene))
        for keyword in payload.get("keywords", []) or []:
            keywords.add(str(keyword))

    if not keywords:
        keywords.update(DEV_KEYWORDS)
    return sorted({keyword for keyword in keywords if keyword}, key=str.lower)


def score_row(row: pd.Series) -> int:
    return score_row_with_keywords(row, load_scoring_keywords())


def score_row_with_keywords(row: pd.Series, keywords: list[str]) -> int:
    score = 0

    source = text_value(row, "source", "sources")
    library = text_value(row, "library_name", "library_names")
    pathway = text_value(row, "pathway", "pathways")
    target = text_value(row, "target", "targets", "target_gene", "target_genes")
    moa = text_value(row, "mode_of_action", "moa", "mechanism", "mechanisms", "actions")

    if "Selleck" in source or "Tocris" in source:
        score += 3

    if "LINCS" in source:
        score += 2

    if "ChEMBL" in source:
        score += 2

    if "IUPHAR" in source:
        score += 2

    text = " ".join([pathway, target, moa, library]).lower()
    for keyword in keywords:
        if keyword.lower() in text:
            score += 2

    return score


def main() -> None:
    master_path = PROCESSED_MASTER_TABLE if PROCESSED_MASTER_TABLE.exists() else OUTPUT_MASTER_TABLE
    if not master_path.exists():
        raise FileNotFoundError(f"Run scripts/02_normalize_merge.py first: {master_path}")

    master = pd.read_csv(master_path, low_memory=False)
    scoring_keywords = load_scoring_keywords()
    master["priority_score"] = master.apply(lambda row: score_row_with_keywords(row, scoring_keywords), axis=1)

    sort_columns = ["priority_score"]
    ascending = [False]
    if "raw_record_count" in master.columns:
        sort_columns.append("raw_record_count")
        ascending.append(False)
    if "canonical_name" in master.columns:
        sort_columns.append("canonical_name")
        ascending.append(True)

    master = master.sort_values(sort_columns, ascending=ascending, na_position="last")
    master.to_csv(OUTPUT_MASTER_TABLE, index=False)

    screen_ready = master.sort_values("priority_score", ascending=False).copy()
    screen_ready.to_csv(SCREEN_READY_LIST, index=False)

    print(f"Read {master_path}, rows={len(master)}")
    print(f"Wrote {OUTPUT_MASTER_TABLE}, rows={len(master)}")
    print(f"Wrote {SCREEN_READY_LIST}, rows={len(screen_ready)}")
    print(f"Top priority score={master['priority_score'].max()}")
    print(f"Scoring keywords={len(scoring_keywords)} from {CONFIG}")


if __name__ == "__main__":
    main()
