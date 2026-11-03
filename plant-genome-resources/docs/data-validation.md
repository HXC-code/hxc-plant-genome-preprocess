# Data Validation Guide

## Prerequisites

- Git Bash at `D:\Git\bin\bash.exe`
- Python available as `python`, or set `PYTHON`
- Canonical data root containing `data/downloads` and
  `resources/download_inventory.tsv`

Do not run the `.sh` file directly from Windows PowerShell. Invoke Git Bash.

## Fresh preprocessing

When derived outputs are absent:

```powershell
cd C:\path\to\HXC-Project

& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --preprocess `
  --refresh-links `
  --expected-count 20
```

This runs:

```text
0check_data.sh
→ Python preprocessing
→ manifest validation and SHA-256 checks
→ compatibility-link checks
→ data-validation-report.tsv
```

## Existing output validation

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' --expected-count 20
```

Without legacy `s<show_id>` directories:

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --no-compat-links `
  --expected-count 20
```

## Incremental preprocessing

```powershell
& 'D:\Git\bin\bash.exe' '.\0check_data.sh' `
  --preprocess `
  --species 996 `
  --refresh-links `
  --expected-count 20
```

Incremental runs preserve manifest and summary rows for unselected species.

## Output layers

```text
data/downloads/    immutable raw input
resources/         inventory, manifests, QC, hashes
clean_v2/          representative-transcript mappings and rejection reasons
analysis_v2/       canonical CDS/PEP/GFF/lens files
results/           future BLAST/WGDI analysis output
```

Standardized data shares the project root with raw data but is not written
inside individual raw species directories.

