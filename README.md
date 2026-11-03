# HXC Plant Genome Preprocess

Reproducible preprocessing and validation for plant genome annotations. The
pipeline converts NCBI-style GFF3, CDS FASTA, and protein FASTA into a
one-gene-one-transcript WGDI-ready dataset.

## Scope

This repository covers:

- raw file inventory and validation;
- GFF3, CDS, and protein parsing;
- canonical cross-species gene IDs;
- representative transcript selection;
- CDS/protein translation consistency checks;
- organelle and taxon contamination filtering;
- WGDI GFF, lens, CDS, and protein generation;
- manifest and SHA-256 validation;
- a reusable Codex skill package.

BLAST, WGDI execution, synteny analysis, phylogenetic trees, and WGD inference
are outside the current preprocessing scope.

## Repository layout

```text
.
├── 0check_data.sh
├── 1run_blast.sh ... 6run_wgdi.sh
├── run_all_pipeline.ps1
├── plant-genome-resources/
│   ├── scripts/preprocess_wgdi_v2.py
│   ├── scripts/validate_wgdi_dataset.py
│   ├── src/pgci_preprocess/
│   ├── tests/
│   ├── docs/
│   └── resources/
└── skill-packages/plant-genome-wgdi-preprocess/
```

Raw genome files and generated analysis results are intentionally ignored.

## Expected data layout

Place raw downloads and the inventory under:

```text
plant-genome-resources/
├── data/downloads/{show_id}_{Order}_{Family}({Species})/
└── resources/download_inventory.tsv
```

Each species directory must contain one GFF, one CDS FASTA, and one protein
FASTA.

## Standardize and validate

From the repository root:

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --preprocess `
  --refresh-links `
  --expected-count 20
```

On a project that does not need root compatibility links:

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --no-compat-links `
  --expected-count 20
```

The project-specific guide is in
`plant-genome-resources/docs/data-validation.md`.

## Test

```powershell
cd plant-genome-resources
python -m pytest -q
ruff check src tests scripts
mypy src --ignore-missing-imports
```

## Codex skill

The reusable skill package is under:

```text
skill-packages/plant-genome-wgdi-preprocess/
```

Install it by copying that directory into:

```text
%USERPROFILE%\.codex\skills\
```

Example prompt:

```text
Use $plant-genome-wgdi-preprocess to standardize and validate this plant genome dataset.
```

