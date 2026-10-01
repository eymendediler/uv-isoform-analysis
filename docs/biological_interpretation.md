# Biological interpretation of primary protein-changing DTU candidates

## Scope and evidence labels

The current results support **time-associated differential transcript usage
(DTU)** after UV in HeLa-S3 cells. They do not yet establish that protein
abundance changes, that both proteins are produced, or that the change is a
causal DNA damage response mechanism.

Evidence is separated into three levels:

1. **Observed here:** statistical DTU and descriptive direction in this dataset.
2. **Biological plausibility:** published gene/pathway connection to stress or
   DNA damage response (DDR).
3. **Direct prior isoform evidence:** a publication tested the relevant protein
   isoforms. Gene-level DDR evidence alone is not labelled isoform validation.

## Strongest interpretable candidates

### RUNX1 — transcript identities confirmed

- **Observed here:** the 453-aa RUNX1-202 share falls from 0.76 at 0 min to
  0.33 at 60 min, while the 250-aa RUNX1-203 share rises from 0.05 to 0.57.
  The 60-vs-0 contrasts have global BH q = 0.0073.
- **Identity validation:** `ENST00000344691.8 / ENSP00000340690.4` maps to
  RefSeq `NM_001001890.3 / NP_001001890.1` and the 453-aa AML1b/RUNX1b protein.
  `ENST00000358356.9 / ENSP00000351123.5` maps to RefSeq
  `NM_001122607.2 / NP_001116079.1`, the reviewed 250-aa AML1a/RUNX1a protein.
- **Known isoform biology:** RUNX1a retains the Runt DNA-binding domain but
  lacks the C-terminal transactivation/regulatory regions of RUNX1b.
- **Interpretation:** the RNA result is therefore a confirmed RUNX1b-to-RUNX1a
  transcript-usage shift, not merely a length-based resemblance. It still does
  not prove a RUNX1b-to-RUNX1a protein-abundance switch.
- **References:** https://www.ncbi.nlm.nih.gov/gene/861,
  https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?DATA=CCDS42922&REQUEST=CCDS,
  and https://pmc.ncbi.nlm.nih.gov/articles/PMC3433167/

### E2F3

- **Observed here:** 465-aa E2F3-201 increases from 0.54 to 0.93 by 60 min;
  334-aa E2F3-202 decreases from 0.46 to 0.07. Both 30- and 60-min contrasts
  pass global BH correction.
- **Known DDR biology:** E2F3a has been shown to be induced after genotoxic
  damage and to support DNA-damage-induced apoptosis; a checkpoint-kinase
  phosphorylation site described in that study is specific to E2F3a rather
  than E2F3b.
- **Interpretation:** this is biologically compatible with isoform-selective
  checkpoint regulation. However, the study used other DNA-damaging agents and
  protein isoform assays; it does not directly prove that our GENCODE 201/202
  shift is UV-driven E2F3a/E2F3b switching.
- **Reference:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2798461/

### ERN1 / IRE1alpha

- **Observed here:** the 977-aa ERN1-201 share decreases from 0.50 to 0.08 at
  60 min; the 70-aa ERN1-208 share rises from 0.40 to 0.92. The 60-vs-0
  contrasts have global BH q = 0.037.
- **Known DDR biology:** full-length IRE1alpha/ERN1 is an ER-stress sensor with
  kinase and endoribonuclease functions. Genotoxic stress can activate its
  regulated IRE1-dependent decay (RIDD) output and influence genome stability,
  cell survival and cell-cycle control.
- **Interpretation:** loss of full-length transcript share and gain of a very
  short predicted protein is intriguing, but there is no direct evidence here
  that the 70-aa product is stable or functional. This candidate needs strong
  junction/read support and protein validation before being called an ERN1
  protein isoform switch.
- **Reference:** https://www.nature.com/articles/s41467-020-15694-y

### EPC1

- **Observed here:** the 813-aa EPC1-202 share falls from 0.44 to 0.05, while
  the 57-aa EPC1-205 share rises from 0.11 to 0.58; 30- and 60-min contrasts
  are significant after global BH correction.
- **Known DDR biology:** EPC1 is a non-catalytic component of the NuA4/TIP60
  acetyltransferase complex. EPC1 recruits MBTD1; this complex participates in
  chromatin signaling and DNA double-strand-break repair-pathway choice.
- **Interpretation:** the gene is directly plausible for chromatin/repair
  biology, but the 57-aa predicted product is unlikely to substitute for the
  characterized full-length EPC1 without direct evidence. A shift toward a
  short coding prediction could represent truncated protein production,
  unstable RNA, or short-read allocation ambiguity.
- **References:** https://pubmed.ncbi.nlm.nih.gov/32209463/ and
  https://pubmed.ncbi.nlm.nih.gov/27153538/

### FER

- **Observed here:** the full-length 822-aa FER-201 share decreases at all UV
  times; 411-aa FER-202 and 163-aa FER-210 shares increase at selected times.
- **Biological interpretation:** FER is a non-receptor tyrosine kinase, but the
  present literature search did not find direct evidence that these exact FER
  isoforms form a UV-DDR switch. The strong reciprocal composition change makes
  it an analytical candidate, not yet a damage-response mechanism.
- **Priority:** validate transcript-specific junctions and ask whether the short
  products retain the kinase domain before assigning functional meaning.

### MECP2

- **Observed here:** 486-aa MECP2-201 decreases sharply at 12 and 30 min while
  172-aa MECP2-216 increases; both time points pass global BH correction.
- **Known isoform biology:** canonical MECP2 alternative splicing produces
  proteins with distinct N termini, and the major historic isoforms have been
  experimentally detected. The 172-aa GENCODE product here is not automatically
  equivalent to those canonical forms.
- **Interpretation:** this supports an early chromatin-regulator transcript
  shift, but published direct UV recruitment evidence for MeCP2 is limited and
  does not validate this exact isoform pair.
- **Reference:** https://pmc.ncbi.nlm.nih.gov/articles/PMC390342/

### PDK3

- **Observed here:** 406-aa PDK3-201 increases at 12 and 30 min while 415-aa
  PDK3-203 decreases; both contrasts pass global BH correction.
- **Biological plausibility:** PDK3 regulates pyruvate dehydrogenase and thus
  mitochondrial carbon flux. Metabolic reprogramming is plausible during acute
  stress, but the search did not identify direct evidence for this precise
  PDK3 isoform pair in UV-DDR.
- **Interpretation:** treat as a metabolism/stress hypothesis rather than a
  known DNA-repair isoform switch.

### NUFIP2

- **Observed here:** 695-aa NUFIP2-201 decreases and 120-aa NUFIP2-202
  increases at 60 min (global BH q = 0.034).
- **Interpretation:** the reciprocal pattern is statistically localized, but
  direct UV-DDR isoform literature was not found. The 120-aa predicted product
  particularly requires validation of junction support and protein stability.

## Omnibus candidates without localized global-BH contrast

- **ELK4:** reciprocal 405/431-aa transcript usage is present descriptively,
  but none of the individual time-vs-0 tests passed the conservative global
  correction. ELK4 is a transcription factor and biologically plausible stress
  responder; the current evidence does not identify one decisive time point.
- **FGF2:** the 288-aa transcript decreases and 155-aa transcript rises most at
  12 min. FGF2 is known to produce low- and high-molecular-weight proteins via
  alternative translation initiation, so transcript length alone cannot be
  directly mapped onto the classic protein isoforms. Reference:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11009566/
- **IER3:** an immediate-early stress-response gene with an early descriptive
  shift, but no individual global-BH-significant time contrast.
- **SFT2D2:** strong descriptive reciprocal usage but no localized contrast at
  the conservative global threshold; little direct UV-DDR isoform evidence was
  identified.

## Is UV-driven isoform switching expected?

Yes at the global mechanistic level, but not yet proven for each candidate.
UV-induced DNA lesions can slow RNA polymerase II elongation and alter
co-transcriptional alternative splicing. Published work has shown widespread UV
splicing changes and functional examples such as BCL-X and ASCC3. Therefore a
time-dependent transcript-usage response is biologically credible.

This does not mean every statistically changing transcript is a functional
DDR isoform. The strongest conclusion is:

> Short-read RNA-seq identifies UV-time-associated, protein-sequence-changing
> transcript-usage candidates. RUNX1 and E2F3 have particularly informative
> prior isoform biology; ERN1 and EPC1 have strong gene/pathway-level DDR
> connections; several short predicted products require direct validation.

Key references:

- UV, RNA polymerase II elongation and alternative splicing:
  https://doi.org/10.1016/j.cell.2009.03.010
- Functional UV-induced ASCC3 coding-to-noncoding isoform response:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5332558/
- Review of reciprocal alternative-splicing/DDR regulation:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7197977/
