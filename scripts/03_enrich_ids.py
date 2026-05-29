#!/usr/bin/env python3
"""Enrich merged perturbagens with cross-database identifiers."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ROOT / "data" / "processed"
ENRICHED_TABLE = PROCESSED_DIR / "perturbagens_enriched.csv"

HEADER = [
    "compound_name",
    "source",
    "source_id",
    "inchikey",
    "smiles",
    "target",
    "mechanism",
    "chembl_id",
    "pubchem_cid",
    "drugbank_id",
    "vendor_ids",
]


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    if not ENRICHED_TABLE.exists():
        ENRICHED_TABLE.write_text(",".join(HEADER) + "\n", encoding="utf-8")
    print(f"Wrote/verified: {ENRICHED_TABLE}")


if __name__ == "__main__":
    main()
