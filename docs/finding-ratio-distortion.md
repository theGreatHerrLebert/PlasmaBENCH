# Finding: from_findings distorts simulated A/B ratios (per-sample median normalization)

## Symptom
Stage-2 (YE spike-in superimposed on real plasma) DIA-NN ratios looked wrong: E. coli
appeared *over-separated* (observed A/B ≈ +1.6 vs nominal +1.0) and yeast *compressed*
(≈ +0.67 vs nominal +1.58 in log2 A/B), with large protein-level spread / apparent
bimodality.

## Diagnosis: the simulated BLUEPRINT truth A/B is itself off-nominal
Measured median truth_intensity (= peptide.events × ion.relative_abundance) per species,
B/A, directly from the `synthetic_data.db` blueprints (nominal: yeast 3.0, E. coli 0.5,
human 1.0):

| | yeast B/A | E. coli B/A | human B/A |
|---|---|---|---|
| Stage 1 (human+YE seed) | 2.19 | 0.38 | **0.76** |
| Stage 2 (YE-only seed)  | 2.06 | 0.35 | — (no human in seed) |

A **single factor f** explains every species: yeast = 3.0·f, E. coli = 0.5·f, human = 1.0·f.
Stage-1 human = 0.76 (must be 1.0 by design) is the smoking gun for a uniform A/B rescale.

## Mechanism (confirmed)
TimSim `load_findings` (rustims) computes, per sample:
```
events = total_intensity / median_intensity * upscale_factor * intensity_multiplier
```
where `median_intensity` is the **per-sample median** over that sample's findings
(load_findings.py:469/484/529/543). For a peptide, events_B/events_A =
(int_B/int_A) · (median_A/median_B). The seed int_B/int_A is the nominal ratio, but
`median_A ≠ median_B` because Sample A and B have different intensity distributions
(different spike-in fractions → different medians). So every A/B ratio is scaled by
median(A)/median(B).

Verified directly from the seed CSVs:
- Stage 2: median(seedA)/median(seedB) = **0.707** → predicts yeast 3.0·0.707=2.12 (meas 2.06),
  E. coli 0.5·0.707=0.35 (meas 0.35).
- Stage 1: 0.743 → yeast 2.23 (meas 2.19), E. coli 0.37 (meas 0.38), human 0.74 (meas 0.76).

## Why Stage 1 looked fine but Stage 2 didn't
- The ratio panels human-anchor the *observed* ratios (subtract the human-median log2 B/A).
  The distortion is a **uniform** factor across ALL species incl. human, so anchoring on
  human cancels it exactly → Stage-1 anchored ratios come out ~nominal despite distorted truth.
- Stage 2 uses a YE-only seed (no human in truth) and anchors on the **real plasma** human
  (genuinely 1:1, identical A/B frames). That anchor (≈0) does NOT match the YE seed's
  distortion (×0.707), so the YE ratios are exposed as off-nominal.

## Consequence for the prior interpretation
- The ratio panels compared observed to the **nominal** expected ratios (3.0/0.5). That is the
  wrong reference — it is not what was simulated. Against the *actual* simulated (blueprint)
  truth, E. coli observed (+1.6) ≈ distorted truth (+1.5) → ~accurate; yeast observed (+0.67)
  vs distorted truth (+1.04) → mildly compressed.
- The **bimodality is separate**: a uniform rescale shifts the center, it does not create two
  modes. The compressed-to-1:1 yeast mode appears to be a real plasma-signal-overlap effect on
  top of the truth distortion.

## Proposed fixes (seeking review before locking one in)
1. **Upstream (rustims load_findings):** normalize A and B by a **shared** median (or a fixed
   anchor), not per-sample, so cross-run ratios are preserved. Correct for any A/B-ratio
   benchmark; quietly affects Stage 1 too (just hidden by anchoring).
2. **Scoring-side:** compare observed ratios to the **per-precursor blueprint truth ratio**
   (the actual oracle we already have) instead of nominal. Makes the panels self-consistent
   without changing the sim.
3. **Seed-side:** pre-scale A/B seeds to equal medians before TimSim so load_findings does not
   distort.

## Questions for review
1. Is the root-cause diagnosis correct, and is the median(A)/median(B) math right?
2. Is the "human-anchoring cancels the uniform distortion in Stage 1" argument valid?
3. Is fix #2 (score against blueprint truth) legitimate, or does scoring observed-vs-blueprint
   *mask* a genuine simulator defect that fix #1 should address instead?
4. Any subtlety we're missing — e.g., is the per-sample median normalization intentional/correct
   for some purpose, and would a shared median break something else (absolute event scale,
   detectability)?
5. Does this interact with the imposed-mode design (intensities = Tenzer baseline × species
   fraction) such that the "right" fix differs?
