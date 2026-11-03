# PDF Manual Steps to Automated Pipeline

The original preprocessing document describes manual `1.py` through `7.py`.
This skill preserves its data intent while replacing fragile manual steps.

| Manual step | Original action | Automated equivalent |
|---|---|---|
| `1.py` | Extract GFF gene/mRNA rows | `gff.py` parses region/gene/transcript/CDS relationships |
| `2.py` | Rename sequence IDs | Numeric sequence IDs generated deterministically |
| `3.py` | Add species prefix and gene order | Globally unique gene IDs generated from species code and order |
| `4.py` | Build lens | Region lengths from GFF plus representative gene counts |
| `5.py` | Extract CDS header IDs | Protein ID extracted from NCBI CDS header attributes |
| `6.py` | Clean CDS header text | ID crosswalk and canonical FASTA headers |
| `7.py` | Replace CDS/PEP sequence names | Atomic writing of filtered CDS/PEP with canonical IDs |

Important differences:

- The manual workflow appends output files and can duplicate data when rerun.
  The automated workflow writes to staging directories and atomically replaces
  completed outputs.
- The manual workflow deletes or edits intermediate files. This skill never
  modifies raw downloads.
- Manual chromosome trimming is replaced by primary chromosome selection and a
  documented gene coverage fallback.
- Gene/mRNA labels are not split with one hard-coded prefix rule; GFF3
  attributes are parsed and URL-decoded.
- CDS/PEP matching is validated by actual translation, not filename order.

