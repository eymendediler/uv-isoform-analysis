# Preliminary Isoform Data

## What We Can Say Now

- We built an isoform-level RNA-seq analysis from GEO Salmon `quant.sf.gz` files.
- We mapped transcript isoforms to GENCODE v35 gene names and protein IDs.
- We separated two questions: DTE for expression change and DTU for isoform usage change.
- These are preliminary candidates because the current exact permutation test has limited resolution with 3 vs 3 replicates.

## DTE Protein-Coding Candidate Counts

| time_min | direction | count |
| --- | --- | --- |
| 12 | down | 1113 |
| 12 | up | 64 |
| 30 | down | 574 |
| 30 | up | 98 |
| 60 | down | 879 |
| 60 | up | 237 |

## DTU Protein-Coding Candidate Counts

| time_min | direction | count |
| --- | --- | --- |
| 12 | usage_down | 554 |
| 12 | usage_up | 542 |
| 30 | usage_down | 518 |
| 30 | usage_up | 562 |
| 60 | usage_down | 632 |
| 60 | usage_up | 746 |

## Strong DTE Examples

| gene_name | transcript_name | protein_id | time_min | direction | control_mean_TPM | uv_mean_TPM | log2FC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EGR1 | EGR1-201 | ENSP00000239938.4 | 60 | up | 0.4576 | 265.6 | 8.896 |
| JUN | JUN-201 | ENSP00000360266.2 | 60 | up | 1.079 | 245.5 | 7.702 |
| NR4A1 | NR4A1-209 | ENSP00000449587.1 | 60 | up | 0 | 16.71 | 7.393 |
| EGR1 | EGR1-201 | ENSP00000239938.4 | 30 | up | 0.4576 | 90.96 | 7.351 |
| EGR2 | EGR2-201 | ENSP00000242480.3 | 60 | up | 0 | 12.47 | 6.974 |
| ARC | ARC-201 | ENSP00000349022.2 | 60 | up | 0.1182 | 19.23 | 6.469 |
| NR4A1 | NR4A1-201 | ENSP00000243050.1 | 60 | up | 0.3017 | 35.16 | 6.456 |
| JUN | JUN-201 | ENSP00000360266.2 | 30 | up | 1.079 | 88.39 | 6.23 |

## Strong DTU Examples

| gene_name | transcript_name | protein_id | time_min | direction | control_mean_usage | uv_mean_usage | delta_usage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AL139260.3 | AL139260.3-204 | ENSP00000479064.1 | 60 | usage_up | 0 | 0.9892 | 0.9892 |
| EGR2 | EGR2-201 | ENSP00000242480.3 | 30 | usage_up | 0 | 0.9118 | 0.9118 |
| CXCL8 | CXCL8-201 | ENSP00000306512.3 | 60 | usage_up | 0 | 0.7788 | 0.7788 |
| EGR2 | EGR2-201 | ENSP00000242480.3 | 60 | usage_up | 0 | 0.7512 | 0.7512 |
| HDAC11 | HDAC11-201 | ENSP00000295757.3 | 30 | usage_down | 0.6919 | 0 | -0.6919 |
| SAPCD1 | SAPCD1-208 | ENSP00000411948.2 | 60 | usage_up | 0 | 0.6319 | 0.6319 |
| HDAC11 | HDAC11-213 | ENSP00000395188.2 | 30 | usage_up | 0 | 0.6014 | 0.6014 |
| SPAG1 | SPAG1-202 | ENSP00000373450.3 | 30 | usage_up | 0.2361 | 0.79 | 0.5539 |

## Recommended Next Steps

1. Validate DTE with a count-based RNA-seq method, ideally tximport plus DESeq2 or edgeR.
2. Validate DTU with a specialist method such as DRIMSeq, DEXSeq, satuRn, or IsoformSwitchAnalyzeR.
3. Focus first on UV response genes and DNA damage response genes: JUN, FOS/FOSB, EGR1/EGR2, ATF3, DUSP family, NR4A1, TP53-related pathways.
4. Make gene-specific isoform plots for 10-20 biologically interesting candidates.
5. Check whether isoform changes are coding, NMD, retained intron, or non-coding transcript changes.
