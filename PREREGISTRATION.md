# Pre-registration of analysis protocols

Status: DRAFT, pending sign-off. This file is written before any model is
fitted and before the freezing experiments are run. No result reported in the
manuscript may come from a protocol not described here. Any deviation is
recorded in a dated amendment at the foot of this file, with the reason, and is
disclosed in the manuscript.

Version 1.1, drafted 2026-08-12.

Part A governs the freezing experiments that carry the paper. Part B governs
the laboratory spectral library analyses that provide the national-scale
support leg. Part C is common to both.

---

# Part A: the freezing experiments

## A1. Question and quantity

We estimate unfrozen water content as a function of temperature, that is the
soil freezing characteristic curve, from reflectance spectra alone, and we
compare that estimate against a contact reference method measured on the same
specimen at the same time.

## A2. Specimens

Specimens are remolded and prepared at prescribed dry density and gravimetric
water content, so that composition is controlled rather than inherited. The
treatment grid spans moisture, organic content, and salinity, and is fixed in a
configuration file before the first specimen is prepared. Target ranges are set
to cover the values observed in the site archive rather than convenient values.

Intact core material, where used, is treated as a separate and clearly labeled
specimen class and is never pooled with remolded specimens in a single fit.

## A3. The unit of independence

One specimen yields many spectra, because each specimen is measured repeatedly
along a temperature ladder, and at each step with replicate scans. The
specimen is the unit of statistical independence. Every split, fold,
permutation, and bootstrap operates on specimens. A specimen appears in exactly
one fold. Replicate scans never straddle a fold boundary. This rule has no
exceptions, and it is the direct lesson of the replicate structure found in our
earlier work.

## A4. Temperature ladder and hysteresis

The ladder is specified in advance: the temperature set points, the
equilibration criterion that must be met before a spectrum is recorded, and the
dwell time. Freezing and thawing branches are recorded as separate branches and
are reported separately. Any claim about hysteresis rests on the two branches
measured on the same specimen, never on branches pooled across specimens.

Specimen temperature is measured in the specimen. Instrument detector
temperature is not a specimen temperature and is never used as one.

## A5. Reference method

A contact reference measurement is recorded at every temperature step on the
same specimen as the spectrum. The reference method is named in the
configuration file before measurement begins. Its own uncertainty at subzero
temperature is characterized and reported, because it sets the ceiling against
which the optical retrieval is judged.

## A6. Optical protocol and artifacts

The illumination and viewing geometry is fixed and reported. Frost and
condensation on the specimen surface and on any window in the optical path are
the expected failure mode. The mitigation is described in the manuscript, and
the criterion for rejecting a contaminated spectrum is defined before
measurement, not after inspection of the data. The number of rejected spectra
is reported.

Spectral regions are not discarded by default. The water and ice absorption
features between 1300 and 1500 nm and between 1800 and 2100 nm carry the signal
of interest, so any masking applied for display is never applied to the data
entering the retrieval, and any band excluded from the retrieval is justified
by a stated instrument limitation.

## A7. Forward model and retrieval

The forward model is specified as equations before fitting, with every free
parameter named and every fixed constant sourced. Optical constants for ice and
liquid water come from a published compilation, cited, and are not tuned.

The retrieval is evaluated in two modes, both reported: with model parameters
fitted per specimen, and with parameters fixed across specimens. The second is
the honest test of transferability and is reported even when it is worse.

## A8. The retrievable regime

The paper states the temperature range and soil classes over which the
retrieval is trustworthy. That boundary is determined by a pre-specified
criterion applied to the reference comparison, not chosen after inspecting
which range looks best. Performance outside the trustworthy regime is reported
rather than omitted.

## A9. Baselines

The retrieval is compared against a statistical model trained on the same
spectra with the same specimen-level grouping. A physics-based contribution
that cannot beat, or at least match with fewer free parameters, a
straightforward statistical baseline is reported as such.

## A10. External validation

Retrieved curves are compared against the published compilation of measured
freezing curves for soils of comparable texture and organic content. The
comparison is descriptive, since those curves come from other soils, and it is
never presented as validation on held-out samples.

---

# Part B: laboratory spectral library analyses

## B1. Scope and provenance

This part governs analyses using the USDA KSSL mid-infrared spectra distributed
through the Open Soil Spectral Library v1.2, retrieved 2026-08-12 from the
public bucket and held with an md5 manifest and a provenance record on
institutional storage. License CC-BY. No file enters an analysis without
matching its recorded md5.

Counts established by audit before any modeling, to be re-derived by the
analysis code rather than copied by hand:

- KSSL layers carrying a non-empty mid-infrared scan: 76,813.
- Distinct site keys, the pedon grouping unit: 20,716, mean 3.71 layers per
  pedon, maximum 643.
- Layers with both depths: 73,730. With a horizon designation: 59,077. With a
  taxonomic name: 36,487.
- Layers whose taxonomy places them in the Gelisol order: 540 across 109
  pedons.
- Layers with mid-infrared spectra and water retention at both 33 kPa and
  1500 kPa: 17,612 across 3,708 pedons, of which 17,000 also carry clay and
  organic carbon. The Gelisol subset is 116 layers across 35 pedons.

## B2. Grouping

The pedon is the unit of statistical independence. Layers from one pedon share
parent material, climate, landscape position, and frequently one laboratory
run. Every split, fold, permutation, and bootstrap operates on pedons.

## B3. Evaluation

Primary: pedon-grouped ten-fold cross-validation, folds assigned by a seeded
permutation of the pedon list, reported for every property.

Secondary, reported alongside the primary and never in place of it:
leave-region-out by Major Land Resource Area; leave-order-out for the Gelisol
contrast, training on non-Gelisol pedons and testing on Gelisol pedons; and a
random layer-level split reported once, labeled optimistic in the text and in
every table caption, retained only for comparability with a literature that
uses it. The gap between that figure and the primary figure is itself reported
as the leakage magnitude.

## B4. Power

With 109 Gelisol pedons overall and 35 carrying retention data, a null result
is uninformative unless power is quantified. The minimum detectable difference
is computed at the pedon level for each property at 80 percent power before the
contrast is run. Contrasts whose minimum detectable difference exceeds a
scientifically meaningful effect are reported as underpowered and are not
interpreted as evidence of similarity.

## B5. Baselines

Partial least squares regression, Cubist, memory-based learning, and a
one-dimensional convolutional network, each tuned on the same folds under the
same grouping. Beating partial least squares alone is not evidence of
contribution.

---

# Part C: common rules

## C1. Label-noise floor

Reported accuracy is meaningless without the ceiling imposed by the reference
measurement. Repeatability is estimated from replicate determinations where the
records permit, and otherwise from published method precision statements, with
the source stated per property. Every accuracy figure is presented against this
floor. A model that reaches the floor is described as having reached it.

## C2. Hyperparameters and selection

Grids are fixed in a configuration file committed before the first run. No
hyperparameter is chosen after seeing test-fold performance. Selection occurs
only within training folds through nested cross-validation. Every cell of every
grid is reported, including cells that failed to converge.

## C3. Significance

Comparisons use a permutation test that shuffles labels at the grouping level,
specimen in Part A and pedon in Part B, never at the level of an individual
spectrum. The number of permutations is fixed in advance at 10,000. Confidence
intervals come from a bootstrap at the same grouping level.

## C4. Reproducibility

Seeds are fixed and recorded. Package versions are captured per run. Every run
writes fold-level results to institutional storage, retained whether or not the
run reaches the manuscript. Figures are regenerated from those files rather
than from transient state.

## C5. Reporting conventions

Coefficients of determination and effect sizes are reported to two decimal
places. Standard deviations across seeds and root-mean-square errors keep their
native precision. Counts, article numbers, and author lists are verified
against the source record before they are written.

## C6. Claim wording

The compilation of measured freezing curves contains dielectric spectroscopy
among its methods. Any priority claim is therefore worded as the first optical
or infrared determination, and never as the first spectroscopic determination.

## Amendments

None to date.
