
# Perturbagen Pipeline

Build a target-aware perturbagen library from LINCS/CMap, ChEMBL, IUPHAR/GtoPdb, and vendor libraries.

The main input is `config/targets.yml`. To retarget the library from one differentiation context to another, edit that file and rerun the pipeline.


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Configure Targets

Edit `config/targets.yml`:

```yaml
run_name: beta_cell

pathways:
  Example_Pathway:
    genes: [GENE1, GENE2]
    keywords: [keyword one, keyword two]
```

The same config drives:

- ChEMBL target lookup
- ChEMBL activity and molecule metadata retrieval
- IUPHAR target and ligand interaction retrieval
- master-table pathway labels
- final `priority_score` ranking

For a new use case, replace the pathways, genes, and keywords. For example: hepatocyte, beta cell, cardiomyocyte, or RPE.

`run_name` controls where the archived results are stored. Use simple names such as `hepatocyte`, `beta_cell`, `cardiomyocyte`, or `rpe`.

## Run

Full refresh after changing `targets.yml`:

```bash
.venv/bin/python scripts/run_pipeline.py
```

Override the run folder name from the command line:

```bash
.venv/bin/python scripts/run_pipeline.py --run-name cardiomyocyte
```

Reuse existing raw source files and only rebuild the merged/scored outputs:

```bash
.venv/bin/python scripts/run_pipeline.py --skip-fetch
```

Refresh only target-dependent public sources, while reusing LINCS/vendor files:

```bash
.venv/bin/python scripts/01_fetch_sources.py --source chembl
.venv/bin/python scripts/01_fetch_sources.py --source iuphar
.venv/bin/python scripts/run_pipeline.py --skip-fetch
```

## Outputs

- `data/processed/perturbagen_source_records.csv`: unified source-level records
- `data/processed/perturbagen_master_table.csv`: deduplicated master table before scoring
- `output/perturbagen_master_table.csv`: scored master table
- `output/stem_cell_screen_ready_list.csv`: scored and sorted screening list

The top-level `data/processed/` and `output/` files always represent the latest run.

Each run is also archived under `runs/<run_name>/`:

```text
runs/
  beta_cell/
    config/targets.yml
    data/
      raw/
        chembl/
        iuphar/
      processed/
        perturbagen_source_records.csv
        perturbagens_merged.csv
        perturbagen_master_table.csv
    output/
      perturbagen_master_table.csv
      stem_cell_screen_ready_list.csv
```

This means changing `targets.yml` for a new cell type and rerunning the pipeline will update the latest files, while preserving a named snapshot for that cell type.

## Deduplication

`02_normalize_merge.py` builds a stable `compound_key` using this priority:

1. InChIKey
2. CAS number
3. normalized compound name
4. source ID

## Scoring

`04_score_for_ipsc.py` adds `priority_score` using:

- source support: Selleck/Tocris, LINCS, ChEMBL, IUPHAR
- pathway, target, mechanism, library text matching
- all genes, pathway names, and keywords from `config/targets.yml`

The scoring script keeps the historical filename `stem_cell_screen_ready_list.csv`, but the ranking itself is target-config driven.
