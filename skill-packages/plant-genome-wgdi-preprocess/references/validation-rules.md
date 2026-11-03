# Validation Rules

## Fail-closed checks

- index and hash manifest exist and are non-empty;
- expected species count matches;
- species codes are unique;
- all four files exist and are non-empty;
- GFF has exactly seven columns;
- sequence IDs are integers and gene IDs are unique;
- strand is `+` or `-`;
- coordinates and order are integers;
- coordinates are positive and sorted;
- order is contiguous from 1 per sequence;
- lens has exactly three columns and unique sequence IDs;
- lens length and gene count are positive;
- coordinates do not exceed lens length;
- lens and GFF sequence sets/gene counts agree;
- GFF, CDS, and PEP ID sets are identical;
- FASTA headers are unique;
- FASTA sequences are non-empty and use allowed alphabets;
- gene IDs are globally unique across species;
- QC JSON exists with `status=OK` and matching code;
- SHA-256 of the four files matches the recorded manifest.

## Optional compatibility checks

When legacy BLAST scripts need root `s<show_id>` directories:

- each compatibility file must exist;
- it must be a hardlink to or have identical content with `analysis_v2`;
- unexpected stale `s<digits>` directories fail validation.

Use `--refresh-links` to recreate compatibility files from the canonical index.
Use `--no-compat-links` for manifest-only projects targeting thousands of
species.

## What validation does not prove

The validator verifies preprocessing structure and integrity. It does not prove:

- assembly completeness;
- absence of biological contamination;
- annotation correctness;
- BUSCO completeness;
- suitability for a particular phylogenetic or WGD question.

Those require separate biological QC and scientific review.

