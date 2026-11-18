# Statistical methods and study designs

## Splice-chain identity

Each local transcript is identified by sample and local isoform ID. The primary
cross-sample identity is chromosome, strand and the ordered junction chain,
including its splice motifs. The paper hash uses the original serialized
`digest(..., algo="xxhash64")` construction. Retain the readable chain string
and test that distinct strings do not share the same hash.

Transcript starts and ends are not part of this primary identity. Different
local-end variants can therefore contribute to one harmonized chain/sample
count. No-junction records receive sample-specific identifiers and cannot by
themselves establish temporal gain/loss. Missing junctions for a multi-exon
transcript are not treated as evidence that it is single-exon.

The core count matrix sums all local transcript counts before attaching a
canonical annotation. Competing gene labels are ranked by FL support, then
detecting-sample support and a deterministic label tiebreak. Annotation
resolution is audited separately from count conservation. The old descriptive
complexity workflow's inner-join restriction is not a reason to discard counts
from the updated core matrix.

## Paired gain/loss and rarefaction

A pair is evaluable for a gene when both specimens meet the gene-level FL floor.
Within that pair, an isoform is present only when it meets both its FL and
within-gene-fraction thresholds. Below-threshold presence and unavailable gene
evaluation are different states. These are threshold-defined detection events,
not proof of complete biological absence.

A cohort event requires the specified same-direction support and no more than
the allowed opposite-direction events. The study requires at least four
same-direction patients and zero opposite events.

Exact rarefaction samples without replacement from the full integer count
matrix. The primary cohort's shallowest library defines the common target.
Candidate filters are applied after rarefaction. Gene fractions retain all
mapped isoforms in their denominator. Every exclusion uses the same primary
depth target.

Depth robustness requires a full-depth cohort call and recovery of the same
cohort call in the configured fraction of all requested completed repetitions.
It is not a marginal per-patient probability and is not conditioned on genes
remaining evaluable. Conditional detection probabilities can be reported as
separate diagnostics.

## Paired technical properties

Technical summaries retain local transcript information for length and terminus
measurements. FL-weighted and unweighted summaries use their stated denominators.
Missing CAGE, poly(A), chromosome or terminus annotation remains unavailable.
The paired association screen uses within-patient changes; sample-level
correlations are descriptive because the specimens are paired.

Paired property comparisons use signed-rank tests. Technical associations use
Spearman correlations where sufficient nonconstant observations are available.
Global and local correction families are retained explicitly. A nonsignificant
association in a small cohort does not establish the absence of technical
influence.

## Paired remodeling metrics and classes

Retained isoforms meet the pooled support rule. A sample is evaluable when its
gene count, retained count and retained fraction meet the configured criteria;
a pair requires both specimens. The isoform-retention mask is fixed from the
full imported cohort for exclusion comparisons.

Condition compositions are unweighted means of patient/sample proportions.
Deeply sequenced specimens do not receive extra weight in these descriptive
compositions. JSD uses exact zeros and base-2 logarithms. Maximum absolute delta
IF is the largest change between condition-mean fractions. Entropy direction
describes changes in usage diversity and is not an additional selection rule.

The effect rule is JSD at least its threshold **or** maximum absolute delta IF
at least its threshold. Clear dominance requires the configured condition-mean
fraction, margin over the next isoform and cross-patient dominance support.

| Class | Meaning |
|---|---|
| A | Clear dominance in both conditions, with a different dominant chain |
| B | Same clear dominant chain, with verified secondary-isoform remodeling |
| C | Clear dominance in the reference condition only |
| D | Clear dominance in the comparison condition only |
| E | Large redistribution without clear dominance in either condition |
| F | Bootstrap-supported absence of an effect of the prespecified size |
| U | Quantifiable, but neither A–E nor supported F |
| NQ | Insufficient data for the configured assessment |

The B criterion removes the retained dominant chain and checks secondary
composition, including minimum secondary mass. A same-dominant gene does not
automatically become B merely because its overall effect is large.

The paired bootstrap resamples complete subjects, preserving each subject's
two specimens. Effects, driver, entropy, dominance and class are recomputed.
The paper uses 5,000 resamples, 95% intervals and a stability probability of
0.80. Class F requires observed metrics below both effect thresholds and the
upper bootstrap limits also below their thresholds. U is not a stable-negative
class.

Observed class, bootstrap stability and cohort-test FDR are separate outputs.
A cohort P value does not create or erase an observed remodeling class.

## Paired sampling noise and cohort simulations

The split-library check partitions each original count into complementary
binomial halves. It evaluates apparent effects caused by counting noise at
reduced depth. The halves and repeated splits are not new biological specimens.
False-call summaries are stratified by depth and isoform complexity, and the
existing selected-gene overlap audit is retained for stressed strata.

The cohort inference uses finite-count sampling covariance in ILR coordinates
and an additive REML between-patient variance. The variance profile propagates
uncertainty beyond a single plug-in estimate. Null draws retain sample depths
and midpoint context and recompute the same joint effect score. Descriptive
metrics continue to use exact zeros; the finite-count model's pseudocount is a
separate model setting.

The 20,000-draw reference is retained. Eligible reference candidates receive
120,000 independently seeded new draws from their frozen fitted null. Selection
and refinement draws are not pooled. Generalized-Pareto estimates describe
the tail of that same null and are not an independent biological validation.

Known-truth simulations compare the fitted variance model with known-parameter
benchmarks. A variance oracle supplies the true generating variance; a full
oracle additionally supplies true midpoint compositions. Null-only QQ and
false-positive summaries evaluate calibration. A real-data P-value QQ alone
does not establish null calibration. Availability and call rate retain
separate denominators.

The patient count-composition test is a within-patient technical likelihood
comparison. It is interpreted with its patient effect gate and is not a
substitute for biological-replicate cohort inference.

## DTU and stageR

DTU uses the original full-depth integer counts, not the rarefied matrix. The
paired model includes subject and condition. The unpaired model includes
condition only, with an intercept-only biological null. Proportion/count filters are
applied per cohort, and the scheduled gene and real-feature families are then
frozen. Unavailable scheduled tests retain conservative analysis P=1 for family
correction and a distinct availability/reason field. They are not described
as tested DTU-negative genes.

The bundled numerical adapter targets DRIMSeq 1.38.0. It retains the corrected
Hessian/log-determinant handling, observed-reference selection, boundary
likelihood limits, stable retries and precision-grid audit. Unexplained
numerical failures block strict completion. See NOTICE for upstream attribution.

StageR is calculated separately for each configured alpha. The complete gene
screen and scheduled transcript family are retained, including unavailable
P=1 cases. Confirmation values obtained at one alpha are not reused at another.

Directional patient support is a paired isoform-fraction shift agreeing with
the selected cohort driver. It is distinct from a per-patient replicated
DRIMSeq experiment. Structural summaries keep their configured full/pre-filter
or tested-feature denominators, rather than combining them silently.

## Sensitivity and computation

Full named exclusions rerun cohort-dependent analyses. Paired exclusions remove
entire subject pairs; unpaired exclusions remove individual independent
samples. Each retained cohort must satisfy the configured sample/support minima.
No minimum is silently reduced after exclusion. Observed-only leave-one-
patient-out recalculates metrics/classes/drivers and is labelled separately.
Selecting all patients creates separate one-patient omissions, not a joint
exclusion or every possible subset. The support requirement remains unchanged.

Seeds, RNG namespace, numerical settings and source version participate in
stochastic cache identity. Atomic saves and completion records prevent a partial
run from being presented as a completed result. Figures and their styling are
separate from the identity of numerical analysis results.

## Paired ORF-free functional enrichment

The final biological remodeling set retains the locked paired-bootstrap and
patient-direction criteria; a separate strict-class set additionally retains
class stability. Cohort-test nominal/FDR subsets use their frozen testing
family as background, whereas biological/class sets use quantifiable genes.
Patient sets use that patient's scheduled technical-testing family and effect
requirement. Unavailable tests stay in their correction/background families.
The set membership and identifier mapping are written for every comparison.

The paired remodeling GO biological-process analysis uses clusterProfiler enrichGO
with the captured org.Hs.eg.db annotation and its compatible GO.db. Its term
sizes are 10–500. Hallmark uses captured MSigDB memberships with term sizes 5–1,000; KEGG
is an editable additional collection and is not enabled for remodeling in the
study preset. At least eight mapped target genes are required. Both reporting
thresholds are retained. Nominal-P fallback figures state that no pathway passed
the selected pathway FDR; nominal terms are not represented as FDR discoveries.
OrgDb, propagated BP membership, identifier maps and database/package versions
are recorded with the outputs. DTU retains its separate MSigDB GO/matched-
background implementation from the source DTU analysis.


## Unpaired independent-sample analyses

Unpaired mode is selected with comparison.paired=false. Each sample must be
an independent biological unit. An omitted/blank subject ID is filled from its
sample ID; explicitly supplied subject IDs must be unique. At least the configured
minimum samples per condition are required (default three). Repeated specimens
and technical replicates are not made independent by assigning new labels.
The implementation does not model extra batch or repeated-measure covariates.

Unpaired detection results compare condition-level detection/prevalence, with
gene-evaluable denominators retained. Read-support and within-gene-fraction
thresholds distinguish detection from unavailable assessment. Depth control
uses independent rarefaction of complete sample count vectors to the main
cohort target. These results describe condition-associated enrichment or
depletion; they are not within-subject gain/loss events.

Unpaired remodeling summarizes condition-mean isoform fractions and their JSD,
maximum fraction difference, total variation and entropy change. Biological
samples are resampled independently within each condition for percentile bootstrap
intervals at `remodeling.bootstrap.bootstrap_confidence` (default 0.95).
Condition-label permutations of the JSD between group means supply the
independent-sample null, with group sizes preserved. All distinct assignments
are enumerated when their number is no greater than n_permutations; otherwise
the configured Monte Carlo draws use a plus-one correction. FDR correction is
applied within each cohort's eligible testing family. Sample eligibility, retained isoforms and unavailable
metrics remain explicit. Exclusion cohorts retain the main cohort's retained
isoform universe to support comparable effects.

The unpaired `usage_shift` effect label requires JSD at least its configured
threshold **or** maximum absolute delta IF at least its configured threshold.
The downstream unpaired remodeling set additionally requires permutation FDR
below the strictest configured reporting alpha. This effect label alone is
not an inferential discovery or a paired architecture class.

The paired architecture classes A–F, paired class-stability probabilities,
within-patient count-composition test and paired directional-support gates
are not assigned to unpaired results. Unpaired outputs use separate effect
and method labels. They do not establish individual longitudinal trajectories.
Independent-sample DTU retains the full-depth integer-count family and stageR
logic but uses condition-only rather than subject-blocked model designs.
Unpaired QC reports depth, detection and independent-condition comparisons; it
does not calculate within-patient technical-change associations.
Default unpaired resampling uses 1,000 bootstrap and 1,000 permutation draws.
A small group size limits the number of distinct permutations and the attainable
P-value resolution even when a larger draw count is requested. Inspect sample
counts, test availability and the actual number of permutations in the outputs.

## Full leave-one-out sensitivity

sensitivity.mode=none adds no automatic omissions. Keep only the main named
cohort for a run with no exclusions; explicit named cohorts are always included.
selected runs one omission
for every named patient (paired) or sample (unpaired); all runs every eligible
unit omission once. Named joint exclusions use exclude_subjects or
exclude_samples in the appropriate design. Cohort group sizes and absolute
support thresholds are validated before computation.

With N exclusion units, rolling leave-one-out adds N analyses to the main run.
Cohort-dependent resampling, inference and DTU filtering/fitting are repeated;
this is computationally different from an observed-metric influence plot.
Resume keys include design, included units, analysis parameters, seeds and
source content so an incompatible model cannot be reused. The main count
harmonization, depth target and isoform-retention reference remain identifiable.
