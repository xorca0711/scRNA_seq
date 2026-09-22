# Retrospective statistical correction protocol

This is a correction and sensitivity analysis written after the original results
were known. It is not a new preregistration and does not create an independent
validation cohort. Frozen historical trial files remain unchanged.

1. Reconstruct the original G1 myeloid and W1 aMAC animal counts and the original
   G2 compartment donor counts. Cache subtype counts separately. Verify exact
   agreement of animal/donor membership, cell numbers and compartment library
   totals with historical tables. Sum duplicate gene symbols explicitly.
2. Use the official `limma::camera` implementation with residual correlation
   estimated (`inter.gene.cor=NA`), negative correlation disallowed, a two-sided
   directional test and biological samples as design-matrix rows. For the W1
   implementation comparison, retain the original log2(CPM+1) transform, original
   gene filter, original `~ batch + sex + arm` design and Mki67 removal. Also run
   official edgeR TMM normalization and limma voom precision weights. The latter
   is the primary corrected RNA-seq analysis for G1/G2; it changes normalization,
   mean-variance modeling and the gene-set statistic from historical GSEA.
3. G1: resolution minus long-term, original unadjusted animals; all animals with
   batch/sex/genotype adjustment; heterozygotes with batch/sex adjustment; and
   December-round heterozygotes with sex adjustment. Test the complete original
   mouse hallmark/GO family (15–500 expressed genes), then report the frozen
   proliferation/arginine targets. Audit addition of sacrifice day to the design
   and the same-round, male-only arm sizes. No age-adjusted biological persistence
   claim is allowed when time/age cannot be separated from phase and round.
4. G2: IPF minus control, three original compartments. Test the complete eligible
   human hallmark/GO family in discovery, and exactly the 254 frozen discovery
   candidates in the validation cohort. Recompute BH across all tested
   compartments/collections within each cohort; do not obtain discovery FDR by
   retesting only previously significant candidates. Require the historical
   direction and corrected global FDR <0.05 in both cohorts for a corrected
   replicated association. These remain observational compartment associations.
5. Leave one donor out of G2 with the full-sample gene universe fixed, refitting
   TMM/voom and CAMERA. Report direction reversals, p-value ranges and omissions
   below the original three-donor floor. LODO candidate-family FDR is a stability
   diagnostic, not a replacement discovery multiplicity correction.
6. Report deposited subtype cell fractions and their shares of each frozen set's
   raw transcripts per donor. Run subtype-specific CAMERA only where both arms
   retain at least three donors with 50 cells per subtype. These subtype analyses
   are exploratory, selected from the frozen candidate list; they do not validate
   a causal explanation or harmonize different annotation systems.

Limitations retained: G1 injury time and age are inseparable, the human metadata
loaded by G2 do not supply age/sex covariates, both human cohorts are observational,
50-cell selection changes the target donor population, subtype labels differ
between cohorts, and ambient RNA/doublets are not resolved by this correction.
Failure to pass a test establishes no equivalence or biological absence.

After the primary correction was inspected, the review added one explicitly
secondary sensitivity: CAMERA's default fixed inter-gene correlation of 0.01,
using exactly the same counts, designs, gene universes and multiplicity families.
This checks dependence on correlation estimation. It does not replace the primary
estimated-correlation analysis or select whichever setting produces significance.

Reference methods: [limma CAMERA source](https://github.com/bioc/limma/blob/master/R/geneset-camera.R),
[limma package](https://bioconductor.org/packages/release/bioc/html/limma.html),
[edgeR package](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
