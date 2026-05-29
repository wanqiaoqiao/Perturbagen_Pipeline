#!/usr/bin/env python3
"""Fetch and normalize raw perturbagen source files."""

import argparse
from pathlib import Path
from typing import Iterable, Optional
from urllib.parse import urljoin

import pandas as pd
import requests
import yaml
from bs4 import BeautifulSoup


def log(message: str) -> None:
    print(message, flush=True)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CONFIG = ROOT / "config" / "targets.yml"
RAW.mkdir(parents=True, exist_ok=True)

HEADERS = {
    "User-Agent": "Mozilla/5.0 perturbagen-research-pipeline/0.1",
}
BIOACTIVITY_TYPES = {
    "AC50",
    "EC50",
    "IC50",
    "Ki",
    "Kd",
    "Potency",
    "XC50",
}
MAX_ACTIVITIES_PER_TARGET = 1000
CHUNK_SIZE = 100
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/"
CHEMBL_TIMEOUT = 60
INCLUDED_CHEMBL_TARGET_TYPES = {"SINGLE PROTEIN"}
GTOP_BASE_URL = "https://www.guidetopharmacology.org/services/"
GTOP_TIMEOUT = 60
SELLECK_L2100_PAGE = "https://www.selleckchem.com/screening/stem-cell-compound-library.html"
SELLECK_L2100_EXPECTED_COMPOUNDS = 1243
CHEMBL_COLLECTION_KEYS = {
    "activity": "activities",
    "drug_mechanism": "drug_mechanisms",
    "mechanism": "mechanisms",
    "molecule": "molecules",
    "target": "targets",
}


def download(url: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        log(f"Exists: {out}")
        return out

    log(f"Downloading: {url}")
    with requests.get(url, headers=HEADERS, stream=True, timeout=120) as response:
        response.raise_for_status()
        with open(out, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return out


def fetch_lincs_geo() -> dict[str, Path]:
    urls = {
        "GSE92742_pert_info": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz",
        "GSE92742_sig_info": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz",
        "GSE70138_pert_info": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_pert_info.txt.gz",
        "GSE70138_sig_info": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
    }

    files = {}
    for name, url in urls.items():
        files[name] = download(url, RAW / "lincs" / f"{name}.txt.gz")
    return files


def normalize_lincs_pert_info(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    df.columns = [column.strip() for column in df.columns]

    # Different LINCS releases use slightly different field names.
    out = pd.DataFrame(index=df.index)
    out["source"] = "LINCS/CMap"
    out["source_file"] = path.name
    out["source_id"] = df.get("pert_id")
    out["canonical_name"] = df.get("pert_iname", df.get("cmap_name"))
    out["perturbagen_type"] = df.get("pert_type")
    out["inchi_key"] = df.get("inchi_key")
    out["canonical_smiles"] = df.get("canonical_smiles")
    out["moa"] = df.get("moa")
    out["target"] = df.get("target")
    out["vendor"] = None
    out["library_name"] = None
    return out.dropna(how="all")


def fetch_and_normalize_lincs() -> None:
    files = fetch_lincs_geo()

    lincs_tables = []
    for key, file in files.items():
        if "pert_info" in key:
            lincs_tables.append(normalize_lincs_pert_info(file))

    lincs = pd.concat(lincs_tables, ignore_index=True).drop_duplicates()
    out = RAW / "lincs" / "lincs_perturbagens_normalized.csv"
    lincs.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(lincs)}")


def load_target_genes() -> pd.DataFrame:
    with open(CONFIG, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    rows = []
    for pathway, payload in config.get("pathways", {}).items():
        for gene in payload.get("genes", []):
            rows.append({"pathway": pathway, "gene": str(gene).upper()})
    return pd.DataFrame(rows).drop_duplicates()


def chunks(values: list[str], size: int = CHUNK_SIZE) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def component_genes(target_record: dict) -> set[str]:
    genes = set()
    for component in target_record.get("target_components") or []:
        for synonym in component.get("target_component_synonyms") or []:
            value = synonym.get("component_synonym")
            if value:
                genes.add(str(value).upper())
    return genes


def serialize_record(record: dict) -> dict:
    row = {}
    for key, value in record.items():
        if isinstance(value, (dict, list)):
            row[key] = yaml.safe_dump(value, sort_keys=True, default_flow_style=True).strip()
        else:
            row[key] = value
    return row


def chembl_get_collection(
    collection: str,
    params: dict,
    max_records: Optional[int] = None,
    session: Optional[requests.Session] = None,
) -> list[dict]:
    session = session or requests.Session()
    url = urljoin(CHEMBL_BASE_URL, f"{collection}.json")
    rows = []
    limit = min(max_records or 1000, 1000)
    request_params = {"limit": limit, **params}

    while url:
        response = session.get(url, params=request_params, headers=HEADERS, timeout=CHEMBL_TIMEOUT)
        response.raise_for_status()
        payload = response.json()
        rows.extend(payload.get(CHEMBL_COLLECTION_KEYS[collection], []))

        if max_records is not None and len(rows) >= max_records:
            return rows[:max_records]

        next_page = payload.get("page_meta", {}).get("next")
        url = urljoin(CHEMBL_BASE_URL, next_page) if next_page else None
        request_params = None
    return rows


def gtop_get_json(
    path: str,
    params: Optional[dict] = None,
    session: Optional[requests.Session] = None,
    empty_on_404: bool = False,
):
    session = session or requests.Session()
    url = urljoin(GTOP_BASE_URL, path)
    response = session.get(url, params=params, headers=HEADERS, timeout=GTOP_TIMEOUT)
    if empty_on_404 and response.status_code == 404:
        return []
    response.raise_for_status()
    return response.json()


def slim_target_record(record: dict, pathway: str, query_gene: str) -> dict:
    genes = component_genes(record)
    accessions = sorted(
        {
            component.get("accession")
            for component in record.get("target_components") or []
            if component.get("accession")
        }
    )
    return {
        "pathway": pathway,
        "query_gene": query_gene,
        "target_chembl_id": record.get("target_chembl_id"),
        "pref_name": record.get("pref_name"),
        "organism": record.get("organism"),
        "tax_id": record.get("tax_id"),
        "target_type": record.get("target_type"),
        "component_genes": "|".join(sorted(genes)),
        "component_accessions": "|".join(accessions),
    }


def fetch_chembl_targets(target_genes: pd.DataFrame) -> pd.DataFrame:
    session = requests.Session()
    rows = []
    for item in target_genes.itertuples(index=False):
        gene = item.gene
        log(f"Searching ChEMBL targets for {gene}")
        seen_ids = set()
        candidates = []
        candidates.extend(
            chembl_get_collection(
                "target",
                {"target_synonym__icontains": gene},
                max_records=50,
                session=session,
            )
        )

        for record in candidates:
            target_id = record.get("target_chembl_id")
            if not target_id or target_id in seen_ids:
                continue
            seen_ids.add(target_id)

            organism = record.get("organism")
            target_type = record.get("target_type")
            genes = component_genes(record)
            if organism != "Homo sapiens" or gene not in genes:
                continue
            if target_type not in INCLUDED_CHEMBL_TARGET_TYPES:
                continue

            rows.append(slim_target_record(record, item.pathway, gene))

    targets = pd.DataFrame(rows)
    if not targets.empty:
        targets = targets.drop_duplicates(subset=["pathway", "query_gene", "target_chembl_id"])
    out = RAW / "chembl" / "chembl_targets.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(targets)}")
    return targets


def fetch_chembl_activities(targets: pd.DataFrame) -> pd.DataFrame:
    if targets.empty:
        activities = pd.DataFrame()
        out = RAW / "chembl" / "chembl_activities.csv"
        activities.to_csv(out, index=False)
        log(f"Wrote {out}, rows=0")
        return activities

    session = requests.Session()
    target_lookup = targets[["pathway", "query_gene", "target_chembl_id"]].drop_duplicates()
    rows = []

    for item in target_lookup.itertuples(index=False):
        log(f"Fetching ChEMBL activities for {item.query_gene} / {item.target_chembl_id}")
        query = chembl_get_collection(
            "activity",
            {
                "target_chembl_id": item.target_chembl_id,
                "pchembl_value__isnull": "false",
                "standard_type__in": ",".join(sorted(BIOACTIVITY_TYPES)),
            },
            max_records=MAX_ACTIVITIES_PER_TARGET,
            session=session,
        )

        for record in query:
            if record.get("standard_type") not in BIOACTIVITY_TYPES:
                continue
            row = {
                "activity_id": record.get("activity_id"),
                "assay_chembl_id": record.get("assay_chembl_id"),
                "assay_description": record.get("assay_description"),
                "assay_type": record.get("assay_type"),
                "document_chembl_id": record.get("document_chembl_id"),
                "molecule_chembl_id": record.get("molecule_chembl_id"),
                "pchembl_value": record.get("pchembl_value"),
                "standard_relation": record.get("standard_relation"),
                "standard_type": record.get("standard_type"),
                "standard_units": record.get("standard_units"),
                "standard_value": record.get("standard_value"),
                "target_chembl_id": record.get("target_chembl_id"),
                "target_organism": record.get("target_organism"),
                "target_pref_name": record.get("target_pref_name"),
            }
            row["pathway"] = item.pathway
            row["query_gene"] = item.query_gene
            rows.append(row)

    activities = pd.DataFrame(rows)
    if not activities.empty:
        activities = activities.drop_duplicates()
    out = RAW / "chembl" / "chembl_activities.csv"
    activities.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(activities)}")
    return activities


def fetch_chembl_molecules(molecule_ids: list[str]) -> pd.DataFrame:
    session = requests.Session()
    rows = []
    for batch in chunks(molecule_ids):
        log(f"Fetching ChEMBL molecule metadata for {len(batch)} molecules")
        query = chembl_get_collection(
            "molecule",
            {"molecule_chembl_id__in": ",".join(batch)},
            session=session,
        )
        for record in query:
            structures = record.get("molecule_structures") or {}
            properties = record.get("molecule_properties") or {}
            rows.append(
                {
                    "molecule_chembl_id": record.get("molecule_chembl_id"),
                    "pref_name": record.get("pref_name"),
                    "molecule_type": record.get("molecule_type"),
                    "max_phase": record.get("max_phase"),
                    "therapeutic_flag": record.get("therapeutic_flag"),
                    "black_box_warning": record.get("black_box_warning"),
                    "first_approval": record.get("first_approval"),
                    "canonical_smiles": structures.get("canonical_smiles"),
                    "standard_inchi": structures.get("standard_inchi"),
                    "standard_inchi_key": structures.get("standard_inchi_key"),
                    "full_mwt": properties.get("full_mwt"),
                    "alogp": properties.get("alogp"),
                    "psa": properties.get("psa"),
                    "hba": properties.get("hba"),
                    "hbd": properties.get("hbd"),
                    "ro5_violations": properties.get("num_ro5_violations"),
                }
            )

    molecules = pd.DataFrame(rows)
    if not molecules.empty:
        molecules = molecules.drop_duplicates(subset=["molecule_chembl_id"])
    out = RAW / "chembl" / "chembl_molecules.csv"
    molecules.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(molecules)}")
    return molecules


def fetch_chembl_mechanisms(molecule_ids: list[str]) -> pd.DataFrame:
    session = requests.Session()
    rows = []
    for batch in chunks(molecule_ids):
        log(f"Fetching ChEMBL mechanisms for {len(batch)} molecules")
        try:
            query = chembl_get_collection(
                "mechanism",
                {"molecule_chembl_id__in": ",".join(batch)},
                session=session,
            )
        except requests.HTTPError as error:
            if error.response.status_code != 404:
                raise
            query = chembl_get_collection(
                "drug_mechanism",
                {"molecule_chembl_id__in": ",".join(batch)},
                session=session,
            )
        rows.extend(serialize_record(record) for record in query)
    mechanisms = pd.DataFrame(rows)
    if not mechanisms.empty:
        mechanisms = mechanisms.drop_duplicates()

    out = RAW / "chembl" / "chembl_mechanisms.csv"
    mechanisms.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(mechanisms)}")
    return mechanisms


def normalize_chembl_bioactive_molecules(
    targets: pd.DataFrame,
    activities: pd.DataFrame,
    molecules: pd.DataFrame,
    mechanisms: pd.DataFrame,
) -> pd.DataFrame:
    if activities.empty:
        normalized = pd.DataFrame()
    else:
        normalized = activities.merge(molecules, on="molecule_chembl_id", how="left")
        if not mechanisms.empty and "molecule_chembl_id" in mechanisms.columns:
            mechanism_cols = [
                col
                for col in [
                    "molecule_chembl_id",
                    "mechanism_of_action",
                    "action_type",
                    "mechanism_refs",
                ]
                if col in mechanisms.columns
            ]
            mechanism_summary = mechanisms[mechanism_cols].drop_duplicates()
            value_cols = [col for col in mechanism_cols if col != "molecule_chembl_id"]
            mechanism_summary = (
                mechanism_summary.groupby("molecule_chembl_id", as_index=False)[value_cols]
                .agg(lambda values: "|".join(sorted({str(value) for value in values.dropna()})))
            )
            normalized = normalized.merge(mechanism_summary, on="molecule_chembl_id", how="left")

        target_names = targets[["target_chembl_id", "pref_name", "target_type"]].drop_duplicates()
        normalized = normalized.merge(target_names, on="target_chembl_id", how="left", suffixes=("", "_target"))

    out = RAW / "chembl" / "chembl_bioactive_molecules_normalized.csv"
    normalized.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(normalized)}")
    return normalized


def fetch_and_normalize_chembl() -> None:
    target_genes = load_target_genes()
    targets = fetch_chembl_targets(target_genes)
    activities = fetch_chembl_activities(targets)
    molecule_ids = sorted(set(activities.get("molecule_chembl_id", pd.Series(dtype=str)).dropna()))
    molecules = fetch_chembl_molecules(molecule_ids) if molecule_ids else pd.DataFrame()
    mechanisms = fetch_chembl_mechanisms(molecule_ids) if molecule_ids else pd.DataFrame()
    normalize_chembl_bioactive_molecules(targets, activities, molecules, mechanisms)


def fetch_iuphar_targets(target_genes: pd.DataFrame) -> pd.DataFrame:
    session = requests.Session()
    rows = []
    for item in target_genes.itertuples(index=False):
        gene = item.gene
        log(f"Searching IUPHAR targets for {gene}")
        records = gtop_get_json("targets", {"geneSymbol": gene}, session=session, empty_on_404=True)
        for record in records:
            row = {
                "pathway": item.pathway,
                "query_gene": gene,
                "target_id": record.get("targetId"),
                "target_name": record.get("name"),
                "target_abbreviation": record.get("abbreviation"),
                "target_type": record.get("type"),
                "family_ids": "|".join(str(value) for value in record.get("familyIds", [])),
                "subunit_ids": "|".join(str(value) for value in record.get("subunitIds", [])),
                "complex_ids": "|".join(str(value) for value in record.get("complexIds", [])),
            }
            rows.append(row)

    targets = pd.DataFrame(rows)
    if not targets.empty:
        targets = targets.drop_duplicates(subset=["pathway", "query_gene", "target_id"])

    out = RAW / "iuphar" / "iuphar_targets.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(targets)}")
    return targets


def slim_iuphar_interaction(record: dict, pathway: str, query_gene: str) -> dict:
    refs = record.get("refs") or []
    pmids = sorted({str(ref.get("pmid")) for ref in refs if ref.get("pmid")})
    reference_ids = sorted({str(ref.get("referenceId")) for ref in refs if ref.get("referenceId")})
    article_titles = sorted({str(ref.get("articleTitle")) for ref in refs if ref.get("articleTitle")})
    return {
        "pathway": pathway,
        "query_gene": query_gene,
        "interaction_id": record.get("interactionId"),
        "target_id": record.get("targetId"),
        "target_name": record.get("targetName"),
        "target_species": record.get("targetSpecies"),
        "primary_target": record.get("primaryTarget"),
        "ligand_id": record.get("ligandId"),
        "ligand_name": record.get("ligandName"),
        "ligand_context": record.get("ligandContext"),
        "endogenous": record.get("endogenous"),
        "interaction_type": record.get("type"),
        "action": record.get("action"),
        "selectivity": record.get("selectivity"),
        "affinity": record.get("affinity"),
        "affinity_parameter": record.get("affinityParameter"),
        "original_affinity": record.get("originalAffinity"),
        "original_affinity_type": record.get("originalAffinityType"),
        "original_affinity_relation": record.get("originalAffinityRelation"),
        "assay_description": record.get("assayDescription"),
        "assay_conditions": record.get("assayConditions"),
        "reference_ids": "|".join(reference_ids),
        "pmids": "|".join(pmids),
        "article_titles": "|".join(article_titles),
    }


def fetch_iuphar_interactions(targets: pd.DataFrame) -> pd.DataFrame:
    session = requests.Session()
    rows = []
    if targets.empty:
        interactions = pd.DataFrame()
    else:
        for item in targets[["pathway", "query_gene", "target_id"]].drop_duplicates().itertuples(index=False):
            log(f"Fetching IUPHAR interactions for {item.query_gene} / target {item.target_id}")
            records = gtop_get_json(
                f"targets/{int(item.target_id)}/interactions",
                {"species": "Human"},
                session=session,
                empty_on_404=True,
            )
            for record in records:
                rows.append(slim_iuphar_interaction(record, item.pathway, item.query_gene))
        interactions = pd.DataFrame(rows)
        if not interactions.empty:
            interactions = interactions.drop_duplicates()

    out = RAW / "iuphar" / "iuphar_interactions.csv"
    interactions.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(interactions)}")
    return interactions


def fetch_iuphar_ligands(ligand_ids: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    session = requests.Session()
    ligand_rows = []
    structure_rows = []
    for ligand_id in ligand_ids:
        log(f"Fetching IUPHAR ligand metadata for {ligand_id}")
        ligand = gtop_get_json(f"ligands/{ligand_id}", session=session)
        ligand_rows.append(
            {
                "ligand_id": ligand.get("ligandId"),
                "ligand_name": ligand.get("name"),
                "ligand_type": ligand.get("type"),
                "abbreviation": ligand.get("abbreviation"),
                "inn": ligand.get("inn"),
                "approval_source": ligand.get("approvalSource"),
                "approved": ligand.get("approved"),
                "who_essential": ligand.get("whoEssential"),
                "withdrawn": ligand.get("withdrawn"),
                "antibacterial": ligand.get("antibacterial"),
                "immuno": ligand.get("immuno"),
                "malaria": ligand.get("malaria"),
                "labelled": ligand.get("labelled"),
                "radioactive": ligand.get("radioactive"),
            }
        )

        try:
            structure = gtop_get_json(f"ligands/{ligand_id}/structure", session=session)
        except (requests.HTTPError, ValueError) as error:
            if isinstance(error, requests.HTTPError) and error.response.status_code != 404:
                raise
            structure = {"ligandId": ligand_id}
        structure_rows.append(
            {
                "ligand_id": structure.get("ligandId", ligand_id),
                "structure_ligand_name": structure.get("ligandName"),
                "iupac_name": structure.get("iupacName"),
                "smiles": structure.get("smiles"),
                "inchi": structure.get("inchi"),
                "inchi_key": structure.get("inchiKey"),
            }
        )

    ligands = pd.DataFrame(ligand_rows).drop_duplicates(subset=["ligand_id"])
    structures = pd.DataFrame(structure_rows).drop_duplicates(subset=["ligand_id"])

    ligand_out = RAW / "iuphar" / "iuphar_ligands.csv"
    structure_out = RAW / "iuphar" / "iuphar_ligand_structures.csv"
    ligands.to_csv(ligand_out, index=False)
    structures.to_csv(structure_out, index=False)
    log(f"Wrote {ligand_out}, rows={len(ligands)}")
    log(f"Wrote {structure_out}, rows={len(structures)}")
    return ligands, structures


def normalize_iuphar_ligand_target_interactions(
    targets: pd.DataFrame,
    interactions: pd.DataFrame,
    ligands: pd.DataFrame,
    structures: pd.DataFrame,
) -> pd.DataFrame:
    if interactions.empty:
        normalized = pd.DataFrame()
    else:
        normalized = interactions.merge(ligands, on=["ligand_id", "ligand_name"], how="left")
        normalized = normalized.merge(structures, on="ligand_id", how="left")
        target_cols = [
            "target_id",
            "target_abbreviation",
            "target_type",
            "family_ids",
            "subunit_ids",
            "complex_ids",
        ]
        target_summary = targets[target_cols].drop_duplicates()
        normalized = normalized.merge(target_summary, on="target_id", how="left")

    out = RAW / "iuphar" / "iuphar_ligand_target_interactions_normalized.csv"
    normalized.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(normalized)}")
    return normalized


def fetch_and_normalize_iuphar() -> None:
    target_genes = load_target_genes()
    targets = fetch_iuphar_targets(target_genes)
    interactions = fetch_iuphar_interactions(targets)
    ligand_ids = sorted(
        {
            int(ligand_id)
            for ligand_id in interactions.get("ligand_id", pd.Series(dtype=float)).dropna().unique()
        }
    )
    ligands, structures = fetch_iuphar_ligands(ligand_ids) if ligand_ids else (pd.DataFrame(), pd.DataFrame())
    normalize_iuphar_ligand_target_interactions(targets, interactions, ligands, structures)


def fetch_selleck_l2100_links() -> dict[str, str]:
    log(f"Fetching Selleck L2100 page: {SELLECK_L2100_PAGE}")
    response = requests.get(SELLECK_L2100_PAGE, headers=HEADERS, timeout=120)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")

    links = {}
    for anchor in soup.find_all("a", href=True):
        text = " ".join(anchor.get_text(" ", strip=True).split()).lower()
        href = urljoin(SELLECK_L2100_PAGE, anchor["href"])
        if "download the sdf" in text or href.lower().endswith(".sdf"):
            links["sdf"] = href
        if "download the xlsx" in text or href.lower().endswith(".xlsx"):
            links["xlsx"] = href

    missing = {"sdf", "xlsx"} - set(links)
    if missing:
        raise RuntimeError(f"Missing Selleck L2100 download links: {', '.join(sorted(missing))}")
    return links


def fetch_selleck_l2100_files() -> dict[str, Path]:
    links = fetch_selleck_l2100_links()
    out_dir = RAW / "vendors" / "selleck"
    files = {
        "sdf": download(links["sdf"], out_dir / "L2100_Stem_Cell_Signaling_Compound_Library.sdf"),
        "xlsx": download(links["xlsx"], out_dir / "L2100_Stem_Cell_Signaling_Compound_Library.xlsx"),
    }
    return files


def normalize_selleck_l2100_xlsx(path: Path) -> pd.DataFrame:
    sheets = pd.read_excel(path, sheet_name=None)
    rows = []
    for sheet_name, frame in sheets.items():
        frame = frame.dropna(how="all").copy()
        if frame.empty:
            continue
        frame.columns = [str(column).strip() for column in frame.columns]
        if "Cat" not in frame.columns or "Name" not in frame.columns:
            continue
        frame["source"] = "Selleck"
        frame["source_file"] = path.name
        frame["library_id"] = "L2100"
        frame["library_name"] = "Stem Cell Signaling Compound Library"
        frame["library_expected_compounds"] = SELLECK_L2100_EXPECTED_COMPOUNDS
        frame["source_sheet"] = sheet_name
        rows.append(frame)

    if not rows:
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True)
    normalized_columns = {
        "Cat": "catalog_id",
        "Catalog No.": "catalog_id",
        "Catalog No": "catalog_id",
        "Cat. No.": "catalog_id",
        "Cat.No.": "catalog_id",
        "Name": "compound_name",
        "Compound Name": "compound_name",
        "CAS No.": "cas_number",
        "CAS No": "cas_number",
        "CAS Number": "cas_number",
        "CAS": "cas_number",
        "M.w.": "molecular_weight",
        "Purity": "purity",
        "Target": "target",
        "Targets": "target",
        "Pathway": "pathway",
        "Pathways": "pathway",
        "Information": "description",
        "Description": "description",
        "SMILES": "smiles",
        "Smiles": "smiles",
    }
    combined = combined.rename(columns={key: value for key, value in normalized_columns.items() if key in combined.columns})
    return combined.drop_duplicates()


def fetch_and_normalize_selleck_l2100() -> None:
    files = fetch_selleck_l2100_files()
    normalized = normalize_selleck_l2100_xlsx(files["xlsx"])
    out = RAW / "vendors" / "selleck" / "selleck_l2100_normalized.csv"
    normalized.to_csv(out, index=False)
    log(f"Wrote {out}, rows={len(normalized)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=["all", "lincs", "chembl", "iuphar", "selleck_l2100", "vendors"],
        default="all",
        help="Source to fetch and normalize.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.source in {"all", "lincs"}:
        fetch_and_normalize_lincs()
    if args.source in {"all", "chembl"}:
        fetch_and_normalize_chembl()
    if args.source in {"all", "iuphar"}:
        fetch_and_normalize_iuphar()
    if args.source in {"all", "selleck_l2100", "vendors"}:
        fetch_and_normalize_selleck_l2100()


if __name__ == "__main__":
    main()
