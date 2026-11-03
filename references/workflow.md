# Preprocessing Workflow

## Position in the comparative-genomics pipeline

```text
raw resource download
→ annotation normalization and QC
→ WGDI input bundle
→ BLAST
→ WGDI collinearity/Ks
→ phylogenetic context and WGD analysis
```

This skill covers the first two stages. It does not execute BLAST or WGDI.

## Raw data

The canonical input root is the project folder containing:

```text
data/downloads/
resources/download_inventory.tsv
```

`download_inventory.tsv` must contain:

```text
showid
Species
Order
Family
gff
cds
pep
```

Raw `data/downloads` files are immutable.

## Processing sequence

1. Validate that GFF, CDS, and PEP filenames refer to the same assembly
   accession.
2. Parse GFF regions, genes, transcripts, CDS features, strand, coordinates,
   and protein IDs.
3. Parse CDS/PEP FASTA by protein ID and detect duplicate conflicting IDs.
4. Reject organelle transcripts and non-coding transcripts without a protein.
5. Translate every candidate CDS with the standard genetic code and compare it
   with the protein sequence.
6. Reject internal stops, non-multiple-of-three CDS length, invalid alphabets,
   and incompatible translation.
7. Select one representative transcript per gene:
   longest protein, then longest CDS, then lexical transcript ID.
8. Select WGDI sequences:
   primary chromosomes when present; otherwise apply a minimum gene threshold
   and coverage fallback.
9. Assign globally unique gene IDs and numeric sequence IDs.
10. Write seven-column simplified GFF, three-column lens, CDS, and PEP.
11. Validate the complete manifest and SHA-256 hashes.

## Output layers

```text
clean_v2/{code}/
├── {code}.gene_map.tsv
├── {code}.chr_map.tsv
└── {code}.rejected.tsv

analysis_v2/{code}/
├── {code}.cds
├── {code}.pep
├── {code}.new.gff
├── {code}.lens
├── {code}.all.cds.fa
├── {code}.all.pep.fa
├── {code}.all.gff
├── {code}.all.lens
└── {code}.qc.json
```

The four core files without `.all` are the filtered WGDI set and must be used
together.

## Scaling to thousands of genomes

Process one species at a time and release memory after each output. The pipeline
supports comma-separated show IDs and preserves unselected manifest entries.

For parallel shards, run disjoint show-ID sets in separate workspaces or add a
separate merge step for summary/index files. Do not let multiple processes write
the same summary/index concurrently.

At large scale:

- use `--no-compat-links`;
- keep operation logs and shard IDs;
- validate hashes before freezing a dataset;
- build a comparison graph instead of all-vs-all BLAST.

