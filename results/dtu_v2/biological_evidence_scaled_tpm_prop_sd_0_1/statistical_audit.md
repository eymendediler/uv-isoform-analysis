# Statistical audit of the DTU and gene-expression workflow

## Supported choices

- Salmon quantifications originate from the published GSE268349 RNA-seq data,
  generated as paired-end 150-bp libraries aligned to GRCh38/GENCODE v35.
- tximport `scaledTPM` is appropriate for DRIMSeq DTU because the modeled
  within-gene proportions then correspond to isoform abundance proportions
  without transcript-length scaling.
- The categorical full model `~ time` versus intercept-only `~ 1` is a valid
  four-group omnibus test. It does not assume monotonic change and retains all
  12 biological samples.
- DRIMSeq's Dirichlet-multinomial model accounts for overdispersion in
  transcript proportions among biological replicates.
- stageR at alpha 0.05 is a principled gene-screen/transcript-confirmation
  procedure, but its OFDR claim relies on calibrated input p-values.
- The separate DESeq2 gene-level likelihood-ratio test (`~ time` versus `~ 1`)
  matches both the DESeq2 recommendation and the analysis reported in the
  source paper. tximport supplies original estimated counts plus gene-length
  offsets for this analysis.

## Correction made during this audit

The initial DRIMSeq-to-stageR run omitted the official rnaseqDTU workflow's
non-specific post-hoc filter for transcripts whose across-sample proportion SD
is below 0.10. The workflow reports that DRIMSeq can exceed transcript-level
FDR and that this filter improved FDR/OFDR control in its evaluations.

The filter was applied without using time labels, transcript p-values were set
to 1 where proportion SD < 0.10, and stageR was rerun at alpha 0.05.

Corrected intersection across loose, medium and strict information filters:

- 187 stageR-confirmed transcripts (previous unfiltered sensitivity: 470)
- 111 protein-coding transcripts (previous: 245)
- 12 genes with at least two robust coding transcripts encoding genuinely
  distinct amino-acid sequences (previous: 31)

RUNX1, ERN1, FER, EPC1, SFT2D2, E2F3, ELK4, PDK3, MECP2, NUFIP2, FGF2 and IER3
remain in the 12-gene primary protein-isoform candidate set. FOS and NR4A1 do
not survive this stricter primary definition and must be described only as
unfiltered sensitivity findings for isoform usage; their strong gene-level UV
induction remains valid.

## Remaining limitations

- Three biological replicates per time point limit dispersion estimation and
  power. Very small p-values do not remove short-read assignment uncertainty.
- The experiment harvests separate cultures at each recovery time. The public
  article/GEO descriptions do not establish that replicate labels form valid
  matched longitudinal blocks, so the unblocked categorical model is the
  defensible primary analysis. Blocking should not be claimed without lab
  confirmation.
- Intersection across three expression filters is deliberately conservative
  and is a robustness criterion, not an independently calibrated error rate.
- Different GENCODE protein sequences demonstrate coding potential, not
  protein translation or protein-level switching.
- Isoform-specific junction/read support and an independent DTU method remain
  necessary before laboratory prioritization.
