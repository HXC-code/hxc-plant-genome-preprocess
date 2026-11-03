# ADR-0001: Separate Raw Data, Preprocessing, Validation, and Results

**Date**: 2026-11-03
**Status**: accepted
**Deciders**: HXC project team

## Context

The project is moving from 20 test species to thousands. Raw downloads, derived
WGDI inputs, compatibility links, and analysis results were previously mixed
across the HXC project root. A validator that inferred species from `.pep`
files could silently pass an empty or incomplete dataset.

## Decision

Raw data remains immutable under `plant-genome-resources/data/downloads`.
Canonical preprocessing outputs live in `clean_v2` and `analysis_v2`.
`wgdi_v2_index.tsv` is the authoritative species/file manifest. Bash is a thin
orchestration and gate layer; Python owns preprocessing and strict validation.
Root `s<show_id>` files are optional compatibility links, not the source of
truth. Future BLAST, WGDI, synteny, tree, and WGD outputs go under `results/`.

## Alternatives Considered

### Put processed files inside each raw species directory

- **Pros**: One folder per species.
- **Cons**: Mixes immutable source and derived data; regeneration overwrites
  evidence or leaves stale files.
- **Why not**: It breaks provenance and versioned reprocessing.

### Keep root species directories as the only interface

- **Pros**: Existing Bash scripts work without changes.
- **Cons**: Thousands of root directories and hardlinks become stale after
  atomic directory replacement.
- **Why not**: A manifest is safer and scales to thousands of species.

### Implement all preprocessing in Bash

- **Pros**: Fewer runtime dependencies.
- **Cons**: GFF3 relationships, isoforms, translation checks, and atomic output
  become fragile and difficult to test.
- **Why not**: Python already implements and tests the data contract.

## Consequences

### Positive

- Empty, incomplete, stale, and contaminated datasets fail closed.
- Incremental processing can preserve unselected manifest entries.
- Canonical results can be validated without copying sequence data.

### Negative

- Validation requires Python.
- Legacy scripts still need root compatibility links until they consume the
  index directly.

### Risks

- Compatibility links can become stale after reprocessing; `0check_data.sh`
  detects this and `--refresh-links` repairs it.
- Parallel shard runs still need a separate manifest merge step before 3,000
  species are processed concurrently.

