#!/usr/bin/env python3
"""Normalize source-specific perturbagen tables into one merged master schema."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog("rdApp.*")


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "output"

SOURCE_RECORDS = PROCESSED_DIR / "perturbagen_source_records.csv"
MERGED_TABLE = PROCESSED_DIR / "perturbagens_merged.csv"
PROCESSED_MASTER_TABLE = PROCESSED_DIR / "perturbagen_master_table.csv"
OUTPUT_MASTER_TABLE = OUTPUT_DIR / "perturbagen_master_table.csv"

MASTER_COLUMNS = [
    "compound_key",
    "canonical_name",
    "synonyms",
    "sources",
    "source_ids",
    "catalog_ids",
    "inchikey",
    "inchi",
    "smiles",
    "cas_number",
    "perturbagen_type",
    "targets",
    "target_genes",
    "pathways",
    "mechanisms",
    "actions",
    "activity_types",
    "activity_values",
    "activity_units",
    "assay_descriptions",
    "molecular_weight",
    "alogp",
    "hba",
    "hbd",
    "vendor",
    "library_ids",
    "library_names",
    "approved",
    "max_phase",
    "raw_record_count",
]

SOURCE_COLUMNS = [
    "source",
    "source_file",
    "source_id",
    "catalog_id",
    "compound_name",
    "synonyms",
    "inchikey",
    "inchi",
    "smiles",
    "cas_number",
    "perturbagen_type",
    "target",
    "target_gene",
    "pathway",
    "mechanism",
    "action",
    "activity_type",
    "activity_value",
    "activity_units",
    "assay_description",
    "molecular_weight",
    "alogp",
    "hba",
    "hbd",
    "vendor",
    "library_id",
    "library_name",
    "approved",
    "max_phase",
]


def clean_value(value):
    if value is None or pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"} or text == "-666":
        return pd.NA
    return text


def clean_series(series: pd.Series) -> pd.Series:
    return series.map(clean_value)


def normalize_name(value) -> str | None:
    value = clean_value(value)
    if pd.isna(value):
        return None
    text = str(value).lower()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text or None


def inchikey_from_smiles(smiles) -> str | None:
    smiles = clean_value(smiles)
    if pd.isna(smiles):
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


def first_existing(df: pd.DataFrame, names: Iterable[str]) -> pd.Series:
    for name in names:
        if name in df.columns:
            return clean_series(df[name])
    return pd.Series(pd.NA, index=df.index)


def scalar(value, index) -> pd.Series:
    return pd.Series(value, index=index)


def ensure_columns(df: pd.DataFrame) -> pd.DataFrame:
    for column in SOURCE_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA
    return df[SOURCE_COLUMNS]


def load_lincs() -> pd.DataFrame:
    path = RAW / "lincs" / "lincs_perturbagens_normalized.csv"
    if not path.exists():
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    src = pd.read_csv(path)
    out = pd.DataFrame(index=src.index)
    out["source"] = first_existing(src, ["source"]).fillna("LINCS/CMap")
    out["source_file"] = first_existing(src, ["source_file"])
    out["source_id"] = first_existing(src, ["source_id"])
    out["compound_name"] = first_existing(src, ["canonical_name"])
    out["inchikey"] = first_existing(src, ["inchi_key"])
    out["smiles"] = first_existing(src, ["canonical_smiles"])
    out["perturbagen_type"] = first_existing(src, ["perturbagen_type"])
    out["target"] = first_existing(src, ["target"])
    out["mechanism"] = first_existing(src, ["moa"])
    out["vendor"] = first_existing(src, ["vendor"])
    out["library_name"] = first_existing(src, ["library_name"])
    return ensure_columns(out)


def load_chembl() -> pd.DataFrame:
    path = RAW / "chembl" / "chembl_bioactive_molecules_normalized.csv"
    if not path.exists():
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    src = pd.read_csv(path)
    out = pd.DataFrame(index=src.index)
    out["source"] = "ChEMBL"
    out["source_file"] = path.name
    out["source_id"] = first_existing(src, ["molecule_chembl_id"])
    out["compound_name"] = first_existing(src, ["pref_name", "molecule_chembl_id"])
    out["inchikey"] = first_existing(src, ["standard_inchi_key"])
    out["inchi"] = first_existing(src, ["standard_inchi"])
    out["smiles"] = first_existing(src, ["canonical_smiles"])
    out["perturbagen_type"] = first_existing(src, ["molecule_type"])
    out["target"] = first_existing(src, ["target_pref_name", "pref_name_target"])
    out["target_gene"] = first_existing(src, ["query_gene"])
    out["pathway"] = first_existing(src, ["pathway"])
    out["mechanism"] = first_existing(src, ["mechanism_of_action"])
    out["action"] = first_existing(src, ["action_type"])
    out["activity_type"] = first_existing(src, ["standard_type"])
    out["activity_value"] = first_existing(src, ["standard_value", "pchembl_value"])
    out["activity_units"] = first_existing(src, ["standard_units"])
    out["assay_description"] = first_existing(src, ["assay_description"])
    out["molecular_weight"] = first_existing(src, ["full_mwt"])
    out["alogp"] = first_existing(src, ["alogp"])
    out["hba"] = first_existing(src, ["hba"])
    out["hbd"] = first_existing(src, ["hbd"])
    out["approved"] = first_existing(src, ["therapeutic_flag"])
    out["max_phase"] = first_existing(src, ["max_phase"])
    return ensure_columns(out)


def load_iuphar() -> pd.DataFrame:
    path = RAW / "iuphar" / "iuphar_ligand_target_interactions_normalized.csv"
    if not path.exists():
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    src = pd.read_csv(path)
    out = pd.DataFrame(index=src.index)
    out["source"] = "IUPHAR/GtoPdb"
    out["source_file"] = path.name
    out["source_id"] = first_existing(src, ["ligand_id"])
    out["compound_name"] = first_existing(src, ["ligand_name"])
    out["inchikey"] = first_existing(src, ["inchi_key"])
    out["inchi"] = first_existing(src, ["inchi"])
    out["smiles"] = first_existing(src, ["smiles"])
    out["perturbagen_type"] = first_existing(src, ["ligand_type"])
    out["target"] = first_existing(src, ["target_name"])
    out["target_gene"] = first_existing(src, ["query_gene"])
    out["pathway"] = first_existing(src, ["pathway"])
    out["mechanism"] = first_existing(src, ["interaction_type"])
    out["action"] = first_existing(src, ["action"])
    out["activity_type"] = first_existing(src, ["affinity_parameter", "original_affinity_type"])
    out["activity_value"] = first_existing(src, ["affinity", "original_affinity"])
    out["assay_description"] = first_existing(src, ["assay_description"])
    out["approved"] = first_existing(src, ["approved"])
    return ensure_columns(out)


def load_selleck() -> pd.DataFrame:
    path = RAW / "vendors" / "selleck" / "selleck_l2100_normalized.csv"
    if not path.exists():
        return pd.DataFrame(columns=SOURCE_COLUMNS)
    src = pd.read_csv(path)
    out = pd.DataFrame(index=src.index)
    out["source"] = "Selleck"
    out["source_file"] = first_existing(src, ["source_file"]).fillna(path.name)
    out["catalog_id"] = first_existing(src, ["catalog_id"])
    out["source_id"] = out["catalog_id"]
    out["compound_name"] = first_existing(src, ["compound_name"])
    out["synonyms"] = first_existing(src, ["Synonyms"])
    out["smiles"] = first_existing(src, ["smiles"])
    out["cas_number"] = first_existing(src, ["cas_number"])
    out["perturbagen_type"] = first_existing(src, ["Form"])
    out["target"] = first_existing(src, ["target"])
    out["pathway"] = first_existing(src, ["pathway"])
    out["assay_description"] = first_existing(src, ["description"])
    out["molecular_weight"] = first_existing(src, ["molecular_weight"])
    out["alogp"] = first_existing(src, ["ALogP"])
    out["hba"] = first_existing(src, ["HBA_Count"])
    out["hbd"] = first_existing(src, ["HBD_Count"])
    out["vendor"] = "Selleck"
    out["library_id"] = first_existing(src, ["library_id"])
    out["library_name"] = first_existing(src, ["library_name"])
    out["inchikey"] = out["smiles"].map(inchikey_from_smiles)
    return ensure_columns(out)


def compound_key(row: pd.Series) -> str:
    for prefix, column in (
        ("inchikey", "inchikey"),
        ("cas", "cas_number"),
        ("name", "compound_name"),
        ("source", "source_id"),
    ):
        value = clean_value(row.get(column))
        if not pd.isna(value):
            if column == "compound_name":
                value = normalize_name(value)
            if value:
                return f"{prefix}:{value}"
    return "unknown"


def join_unique(values: pd.Series, limit: int = 80) -> str | pd.NA:
    cleaned = []
    for value in values:
        value = clean_value(value)
        if pd.isna(value):
            continue
        for part in str(value).split("|"):
            part = clean_value(part)
            if not pd.isna(part):
                cleaned.append(str(part))
    unique = sorted(set(cleaned), key=lambda item: item.lower())
    if not unique:
        return pd.NA
    if len(unique) > limit:
        unique = unique[:limit] + [f"...(+{len(unique) - limit} more)"]
    return "|".join(unique)


def first_non_null(values: pd.Series):
    for value in values:
        value = clean_value(value)
        if not pd.isna(value):
            return value
    return pd.NA


def merge_sources(records: pd.DataFrame) -> pd.DataFrame:
    records = records.copy()
    records["compound_key"] = records.apply(compound_key, axis=1)
    records = records[records["compound_key"] != "unknown"].copy()

    grouped = records.groupby("compound_key", dropna=False)
    merged = grouped.agg(
        canonical_name=("compound_name", first_non_null),
        synonyms=("synonyms", join_unique),
        sources=("source", join_unique),
        source_ids=("source_id", join_unique),
        catalog_ids=("catalog_id", join_unique),
        inchikey=("inchikey", first_non_null),
        inchi=("inchi", first_non_null),
        smiles=("smiles", first_non_null),
        cas_number=("cas_number", first_non_null),
        perturbagen_type=("perturbagen_type", join_unique),
        targets=("target", join_unique),
        target_genes=("target_gene", join_unique),
        pathways=("pathway", join_unique),
        mechanisms=("mechanism", join_unique),
        actions=("action", join_unique),
        activity_types=("activity_type", join_unique),
        activity_values=("activity_value", join_unique),
        activity_units=("activity_units", join_unique),
        assay_descriptions=("assay_description", join_unique),
        molecular_weight=("molecular_weight", first_non_null),
        alogp=("alogp", first_non_null),
        hba=("hba", first_non_null),
        hbd=("hbd", first_non_null),
        vendor=("vendor", join_unique),
        library_ids=("library_id", join_unique),
        library_names=("library_name", join_unique),
        approved=("approved", join_unique),
        max_phase=("max_phase", first_non_null),
        raw_record_count=("source", "size"),
    ).reset_index()

    return merged[MASTER_COLUMNS].sort_values(["canonical_name", "compound_key"], na_position="last")


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_tables = [load_lincs(), load_chembl(), load_iuphar(), load_selleck()]
    records = pd.concat(source_tables, ignore_index=True)
    records = ensure_columns(records)
    records.to_csv(SOURCE_RECORDS, index=False)

    merged = merge_sources(records)
    merged.to_csv(MERGED_TABLE, index=False)
    merged.to_csv(PROCESSED_MASTER_TABLE, index=False)
    merged.to_csv(OUTPUT_MASTER_TABLE, index=False)

    print(f"Wrote {SOURCE_RECORDS}, rows={len(records)}")
    print(f"Wrote {MERGED_TABLE}, rows={len(merged)}")
    print(f"Wrote {PROCESSED_MASTER_TABLE}, rows={len(merged)}")
    print(f"Wrote {OUTPUT_MASTER_TABLE}, rows={len(merged)}")


if __name__ == "__main__":
    main()
