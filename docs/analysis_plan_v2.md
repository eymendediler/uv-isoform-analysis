# Analysis plan v2: time-dependent protein-coding isoform usage

## Primary biological question

After UV exposure, does the relative usage of protein-coding transcript
isoforms within a gene change over time, and could those transcripts encode
different protein products?

RNA-seq can support a transcript-level isoform switch and its protein-coding
potential. It cannot by itself demonstrate that protein abundance or protein
isoform composition changed; that requires orthogonal validation.

## Locked principles before testing

1. Retain all 12 samples unless a documented technical failure is found.
2. Do not apply a TPM, count, fold-change or delta-usage threshold during the
   Salmon input audit.
3. Preserve deposited `quant.sf.gz` files unchanged.
4. Record every filtering decision and show how conclusions change under
   alternative reasonable filters.
5. Treat time as a four-level categorical variable for the primary omnibus
   test; do not assume a linear response.
6. Separate statistical evidence (FDR) from biological magnitude (delta
   usage), and report both.

## Proposed inferential workflow

1. Audit the 12 deposited Salmon files and sample metadata.
2. Reproduce gene-level time responses as a positive-control analysis.
3. Import Salmon estimated counts and transcript-level uncertainty.
4. Test a gene-level omnibus DTU hypothesis with DRIMSeq, followed by stageR
   transcript confirmation control.
5. Estimate 12, 30 and 60 min versus control contrasts and display complete
   time trajectories. Adjacent-time contrasts are secondary.
6. Repeat the analysis across a predeclared low-expression-filter grid and
   report candidates stable across choices. No single filter becomes primary
   until its retention and sensitivity diagnostics have been reviewed.
7. Use satuRn as a secondary method and SUPPA2 for event-level interpretation,
   not as extra biological replication.
8. Annotate coding candidates by GENCODE transcript/CDS/protein identifiers;
   distinguish transcripts that encode identical proteins from those with
   different CDS or amino-acid products.

## Decisions deliberately left open

- Whether replicate is a valid blocking factor. The model will be compared
  with and without `replicate`; the paired/block interpretation requires
  experimental confirmation from the original laboratory.
- The primary low-expression filter. Candidate grids will be compared before
  selection.
- The delta-usage magnitude used for a “strong candidate” label. Values such
  as 0.05, 0.10 and 0.20 will be summaries, not significance cutoffs.

The primary statistical error target proposed for review is gene-level
FDR <= 0.05 with stage-wise transcript confirmation. This is not yet a final
decision and will be assessed together with power, retained-feature counts and
sensitivity results.

## Exact DRIMSeq filter order

The filter values are not applied simultaneously. DRIMSeq first checks total
gene counts, then removes transcripts that fail the transcript-count rule, and
finally recalculates within-gene proportions among the remaining transcripts
before applying the proportion rule. Genes with fewer than two retained
transcripts are not testable and are removed. An independent implementation of
this sequence reproduced the DRIMSeq retention totals exactly for all three
filter settings.
