# Isoform Validation Outputs

## Replicate Dot Plots

These plots answer:

`Are the isoform TPM changes visible in individual replicates, or only in the mean?`

Each dot is one replicate. The black line connects condition means.

A strong candidate should usually show dots moving in the same direction across most or all replicates.

## ORF/CDS Tables

These tables answer:

`Could the changed transcript isoform encode a different protein product?`

Useful columns:

- `transcript_type`: protein_coding, retained_intron, nonsense_mediated_decay, etc.
- `protein_id`: Ensembl protein ID, present when a protein product is annotated.
- `cds_bp`: total coding sequence length in base pairs from GENCODE CDS records.
- `protein_length_aa_est`: estimated amino acid length from CDS length.

If two isoforms from the same gene have different protein IDs or different estimated protein lengths, the isoform change may affect the encoded protein.

Protein domain effects require an additional domain database such as UniProt, Pfam, or InterPro. The current outputs are ORF/CDS-level preliminary checks.
