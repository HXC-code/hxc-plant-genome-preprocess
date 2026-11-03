---
name: plant-genome-wgdi-preprocess
description: >-
  Standardize NCBI GFF3, CDS FASTA, and protein FASTA into a reproducible
  one-gene-one-transcript WGDI bundle with unique IDs, lens files, validation
  reports, and manifests. Use when preparing plant genome annotations for
  BLAST/WGDI, validating existing WGDI inputs, or safely reprocessing batches
  of plant genomes. This skill ends before BLAST and WGDI execution.
---

# Plant Genome WGDI Preprocess

Use this skill for stage-1 annotation normalization and stage-2 preprocessing:
raw NCBI-style GFF3/CDS/PEP to validated WGDI-ready files.

## Required input

The target data root must contain:

```text
data/downloads/{show_id}_{Order}_{Family}({Species})/
resources/download_inventory.tsv
```

Each species folder must contain one GFF, one CDS FASTA, and one protein FASTA.
Compressed `.gz` inputs are supported.

Never modify or delete files under `data/downloads`.

## Primary workflow

1. Inspect `resources/download_inventory.tsv` and preserve its exact columns.
2. Run one representative species before a full batch.
3. Read `resources/stage_preprocess_v2_summary.tsv` and per-species QC.
4. Inspect rejection reasons before trusting a `wgdi_ready=yes` result.
5. Run manifest-driven validation as a hard gate.

Use the bundled pipeline:

```bash
bash scripts/run_pipeline.sh /absolute/path/to/plant-genome-resources \
  --species 996 \
  --expected-count 20
```

Without `--species`, all inventory species are processed.

For a project that still uses root `s<show_id>` compatibility links:

```bash
bash scripts/run_pipeline.sh /absolute/path/to/plant-genome-resources \
  --refresh-links \
  --expected-count 20
```

Without `--refresh-links`, only canonical `analysis_v2` data is validated.

## Direct commands

```bash
python scripts/preprocess.py \
  --project-root /absolute/path/to/plant-genome-resources \
  --species 996

python scripts/validate_dataset.py \
  --project-root /absolute/path/to/plant-genome-resources \
  --expected-count 20
```

Use `--skip-hashes` only for a quick structural check. Use the default hash
check for a release or BLAST gate.

## Output contract

```text
clean_v2/{code}/{code}.gene_map.tsv
clean_v2/{code}/{code}.chr_map.tsv
clean_v2/{code}/{code}.rejected.tsv

analysis_v2/{code}/{code}.cds
analysis_v2/{code}/{code}.pep
analysis_v2/{code}/{code}.new.gff
analysis_v2/{code}/{code}.lens
analysis_v2/{code}/{code}.qc.json
```

The four WGDI files must be used as a set:

```text
{code}.cds
{code}.pep
{code}.new.gff
{code}.lens
```

Do not combine a filtered GFF with unfiltered CDS/PEP files.

## Non-negotiable rules

- One gene has one representative transcript.
- Select the longest translation-compatible protein.
- Reject internal stops, non-multiple-of-three CDS, empty sequences,
  invalid alphabets, and translation mismatches.
- Exclude organelles from WGDI input.
- Reject multiple nuclear taxon IDs.
- Preserve source gene, transcript, and protein IDs in mapping files.
- Make gene IDs unique across species.
- Write staged outputs; never mutate raw input.
- Fail closed when the manifest is empty, incomplete, stale, or missing.

## References

- Read [references/workflow.md](references/workflow.md) for the complete stage
  sequence and scaling guidance.
- Read [references/wgdi-input-contract.md](references/wgdi-input-contract.md)
  before changing GFF, lens, CDS, or PEP formatting.
- Read [references/validation-rules.md](references/validation-rules.md) before
  changing acceptance behavior.
- Read [references/pdf-to-pipeline-mapping.md](references/pdf-to-pipeline-mapping.md)
  when comparing this implementation with the original manual `1.py`-`7.py`
  workflow.

## Limits

This skill does not run BLAST, WGDI collinearity, Ks, phylogenetic trees, or
WGD event inference. It prepares and validates their inputs.

