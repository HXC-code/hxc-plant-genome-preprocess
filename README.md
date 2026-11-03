# Plant Genome WGDI Preprocess

Codex Skill for converting NCBI-style plant genome GFF3, CDS FASTA, and protein
FASTA into a validated one-gene-one-transcript WGDI input bundle.

## Install

Copy this repository into the Codex user skill directory under a folder matching
the skill name:

```powershell
git clone https://github.com/HXC-code/hxc-plant-genome-preprocess.git
Copy-Item `
  .\hxc-plant-genome-preprocess `
  "$HOME\.codex\skills\plant-genome-wgdi-preprocess" `
  -Recurse
```

## Use

```text
Use $plant-genome-wgdi-preprocess to standardize and validate this plant genome dataset.
```

The skill requires a data project containing:

```text
data/downloads/
resources/download_inventory.tsv
```

Run one species first:

```bash
bash scripts/run_pipeline.sh /path/to/plant-genome-resources \
  --species 996 \
  --expected-count 1
```

For a full dataset, omit `--species` and set the expected species count.

## Test

```powershell
python -m pytest -q
ruff check scripts tests
mypy scripts/pgci_preprocess --ignore-missing-imports
```

## Scope

This skill generates and validates WGDI inputs. It does not run BLAST, WGDI
collinearity, Ks, phylogenetic trees, or WGD inference.

