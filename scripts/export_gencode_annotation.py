#!/usr/bin/env python3
"""Export one-row-per-transcript GENCODE metadata for count-based DTU."""

from pathlib import Path

from analyze_transcript_isoforms import load_transcript_annotation


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "metadata" / "gencode_v35_transcripts.tsv.gz"


def main() -> None:
    annotation = load_transcript_annotation()
    if annotation.empty:
        raise RuntimeError("GENCODE v35 annotation could not be loaded")
    if annotation["transcript_id"].duplicated().any():
        raise RuntimeError("Transcript identifiers are not unique")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    annotation.to_csv(OUT, sep="\t", index=False, compression="gzip")
    protein = annotation.query(
        "gene_type == 'protein_coding' and transcript_type == 'protein_coding'"
    )
    print(f"Wrote {len(annotation):,} transcripts to {OUT}")
    print(f"Protein-coding transcripts: {len(protein):,}")
    print(
        "Genes with >=2 protein-coding transcripts: "
        f"{(protein.groupby('gene_id').transcript_id.nunique() >= 2).sum():,}"
    )


if __name__ == "__main__":
    main()
