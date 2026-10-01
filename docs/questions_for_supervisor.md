# Experimental-design and interpretation questions for the supervisor

This document separates points that cannot be resolved from the public paper
and GEO metadata. The accepted statistical design should be revisited after the
original laboratory confirms these details.

## Are samples independent or experimentally blocked?

The paper states that HeLa-S3 cells were grown in 150-mm culture dishes and
samples were collected at 0, 12, 30 and 60 minutes after UV. RNA isolation is
destructive, so the same cells cannot have been measured repeatedly. However,
the public methods and metadata do not establish whether the four samples
labelled `replicate 1` came from one starting culture or processing day and
therefore form an experimental block.

Questions:

1. Were control, 12-, 30- and 60-minute samples within `replicate 1` separate
   dishes split from the same starting culture and processed on the same day?
2. Were replicates 1, 2 and 3 independent culture starts or independent
   experimental days?
3. Were all four time points treated and processed together within a replicate?
4. Are RNA-extraction date, library-preparation batch, operator, kit lot or
   sequencing-lane metadata available?
5. Does the replicate number encode biological matching across time, or is it
   only an arbitrary sample number within each time point?

Decision rule:

- If replicate number denotes a shared starting culture/day, a
  `~ replicate + time` sensitivity model is scientifically justified.
- If replicate numbers are not matched across time, the current `~ time` model
  is appropriate.
- A block term should not be added solely because sample names share a number.

## Biological hypothesis and validation

6. Is the intended primary hypothesis: “During the first 60 minutes after UV,
   within-gene usage of transcripts capable of encoding different protein
   sequences changes as a function of recovery time”?
7. Should transcripts differing only in untranslated regions be included, or
   should the first validation round focus on distinct amino-acid sequences?
8. Which candidates are experimentally accessible? Current primary candidates
   are RUNX1, ERN1, FER, EPC1, SFT2D2, E2F3, ELK4, PDK3, MECP2, NUFIP2, FGF2
   and IER3.
9. Can isoform-specific exon-exon junction RT-qPCR assays be designed for the
   selected transcripts?
10. Is targeted long-read cDNA or amplicon sequencing feasible for a small
    number of genes to confirm full-length exon connectivity?
11. Is protein-level validation feasible when the predicted products differ
    enough in size or possess isoform-specific peptide sequences?

## Statistical decisions to agree in advance

12. Should the primary stageR target remain `alpha = 0.05` overall false
    discovery rate?
13. Should the official-workflow proportion-SD threshold of 0.10 remain the
    primary safeguard while 0.05, 0.075 and 0.10 are reported as sensitivity
    analyses?
14. Should biological magnitude be reported continuously rather than declaring
    an unvalidated delta-usage discovery cutoff?
15. Should 12-vs-0, 30-vs-0 and 60-vs-0 tests remain secondary localization
    analyses with global Benjamini-Hochberg correction across every tested
    transcript-by-time hypothesis?

## Why short-read results remain candidates

Illumina paired-end 150-bp sequencing provides high depth but does not read one
RNA molecule from end to end. Short fragments from exons shared by several
transcripts may be compatible with multiple isoforms; Salmon allocates this
ambiguous evidence probabilistically. The problem is not that long transcripts
are “split into more isoforms.” The problem is that a short fragment may not
identify the complete exon combination from which it originated.

Long-read sequencing can observe much longer or full-length molecules and
therefore provides more direct exon-connectivity evidence. It also has its own
depth, error and protocol biases. Targeted junction assays and long-read data
are complementary validation strategies.

## References

- Source study: https://www.nature.com/articles/s41467-024-55724-7
- Bioconductor rnaseqDTU workflow:
  https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html
- DRIMSeq: https://pmc.ncbi.nlm.nih.gov/articles/PMC5200948/
- stageR: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1277-0
- Long-read isoform analysis/ESPRESSO:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9858503/
- RNA-seq replication, blocking and batch design; open-access article and PDF
  page: https://pmc.ncbi.nlm.nih.gov/articles/PMC2881125/

