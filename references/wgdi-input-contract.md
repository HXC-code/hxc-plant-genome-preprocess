# WGDI Input Contract

## Simplified GFF

Each line is tab-separated:

```text
sequence_id  gene_id  start  end  strand  order  source_gene_id
```

Requirements:

- `sequence_id` is a simplified integer for this skill's output.
- `gene_id` is globally unique across the processed species set.
- `start` and `end` are positive integers and `start <= end`.
- `strand` is `+` or `-`.
- `order` starts at 1 and is contiguous within each sequence.
- Rows are sorted by sequence ID, start, and order.
- The seventh source ID column is ignored by WGDI but retained for provenance.

This is not normal nine-column GFF3.

## Lens

```text
sequence_id  sequence_length_bp  representative_gene_count
```

Requirements:

- one unique row per sequence;
- positive integer length and gene count;
- gene count equals the number of GFF rows for that sequence;
- all gene ends are within the sequence length.

## CDS FASTA

- Nucleotide FASTA keyed by representative gene ID.
- IDs match the GFF gene IDs exactly.
- Allowed IUPAC nucleotide alphabet is enforced by validation.
- No duplicate or empty records.

## Protein FASTA

- Protein FASTA keyed by representative gene ID.
- IDs match the GFF and CDS IDs exactly.
- No duplicate or empty records.
- Internal stops are rejected during preprocessing.

## BLAST output for the next stage

WGDI expects tab-separated BLAST outfmt 6 with at least:

```text
query_id
subject_id
pident
length
mismatch
gapopen
qstart
qend
sstart
send
evalue
bitscore
```

No header row. Query and subject IDs must exist in the selected WGDI bundle.

