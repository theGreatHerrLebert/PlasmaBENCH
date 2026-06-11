# Stage 2 — YE spike-in simulated over REAL plasma (design + intensity calibration)

## Goal
Stage 1 simulates all three species and adds noise *sampled from a blank* run. We
observed a low-abundance ratio-bias / detection trend that is *reversed* vs what is seen
on real tTOF plasma (real plasma: dense background → DIA-NN background subtraction →
low-abundance signal too HIGH; our blank-noise SIM → too LOW). Stage 2 tests whether
putting the simulated spike-in on a REAL plasma background reproduces the real-data
direction.

## Mechanism
TimSim `superimpose_on_reference=true` (mutually exclusive with `add_real_data_noise`)
overlays simulated signal onto the frames of a real plasma `.d`. So:
- The seed is **YE-only** (yeast + E. coli; human removed via `--keep-species YEAST ECOLI`).
- The simulated YE signal is added on top of the real plasma frames.
- The blueprint (ground truth) therefore contains **only YE** (verified: a Stage-2 sim
  produced YEAST 18997 + ECOLI 13225 + HUMAN 0).
- Human is the **real plasma background** — no blueprint truth → `truth_scope =
  background_unknown`. Scored empirically only (P1/P2 ratio panels); excluded from the
  ground-truth FDR (P5) and sensitivity (P3), which are restricted to the simulated YE.

## Concrete setup (already built)
- Reference plasma: `PYE_plasma_2022_135` (the *same* plasma used to make the PYE mix),
  runs `G250506_080` and `_081`, RT span 2460 s. Real human plasma, no spike-in.
- A superimposed on run 080, B on run 081 (two real injections as the A/B human backgrounds).
- YE-only seeds from Ute's G-site PYE1 DIA-NN report (run-split A1/B1), native RT ≤2121 s
  (fits 2460 s). A: 24319 yeast + 17538 E.coli precursors; B: 25114 + 17138.
- Search the resulting `.d` with DIA-NN 2.5 (unified cfg, fixed 15 ppm); score with the
  existing P1/P5/P6/P3 panels and `--human-scope background_unknown`.
- A/B oracle ratios (yeast B/A=3, E.coli B/A=0.5) are baked into the seed intensities;
  `proteome_mix=false` is enforced so they are never re-diluted.

## The problem this review is about: YE intensity calibration
In Stage 1 the YE intensity scale only had to be internally consistent (blank = sparse
noise). In superimpose, the real human signal is already in the frames **at its true count
scale**, and we add simulated YE counts on top. Currently YE intensities come from an
imposed Tenzer rank-abundance baseline × an **arbitrary `upscale_factor` (1e5)** — NOT
anchored to the plasma's intensity scale. So the absolute YE level relative to human is a
guess: too high → YE trivially detectable; too low → buried. That directly sets the
low-abundance detectability and ratio behavior we are trying to measure, so it must be
calibrated.

## Proposed calibration
1. **Measure the human background** by re-analyzing the two plasma reference runs with
   DIA-NN 2.5 + human-only fasta (running now). Gives the human per-precursor
   rank-abundance distribution (shape) in each reference run, on the same engine/settings
   we score the SIM with.
2. **Absolute scale**: bridge the reported-quant scale to the raw frame-count scale by
   also measuring the plasma `.d` total frame-intensity sum (the counts superimpose adds
   onto), so the YE *count* total lands at the intended level.
3. **Set YE totals by PYE1 mass fraction relative to measured human:** human = 90% mass,
   so yeast total = (2/90)·human in A and (6/90)·human in B; E.coli = (8/90)·human in A,
   (4/90)·human in B. This preserves the A/B oracle ratios automatically and pins the
   absolute scale to the real background.
4. **Within each species**, distribute by Tenzer rank-abundance (as now) so YE peptides
   span a realistic dynamic range *inside* the human range — low-abundance YE lands near
   the human low tail where detection is genuinely hard.
5. Implement as `build_seed_from_report.py --ratio-mode plasma-calibrated
   --plasma-report <2.5 humanfast report> --plasma-d <plasma.d>`, replacing the arbitrary
   upscale with the mass-fraction-relative-to-human scale.

## Revisions after Codex review (2026-06-09) — ADOPTED

- **Calibration is empirical, not mass-fraction→counts.** Mass fraction sets only the
  *target composition*; MS response differs by species (digestion/ionization/fragmentation).
  Fit the YE intensity scale so simulated **raw MS1/MS2 peak areas** reproduce the *real*
  yeast/E.coli intensity distributions from the real PYE mixed runs; freeze on held-out
  runs before scoring. (Ideally needs the raw PYE-mix `.d`; the PYE DIA-NN report is a
  proxy for the relative YE distribution but only an external-validation target for scale.)
- **No global frame-sum bridge.** Calibrate in matched domains (integrated MS1 isotope
  envelope; MS2 fragment intensity in the DIA windows), stratified by m/z/mobility/RT/
  charge/window; check local saturation. DIA-NN quantity = external validation only.
- **Crossed background design:** A over BOTH 080 and 081, B over BOTH (≥2 plasma
  injections), to separate concentration effect from run-specific interference. The single
  A/080–B/081 pairing is confounded.
- **Trace contamination check (do FIRST):** search unmodified 080/081 against
  human+yeast+E.coli; quantify any real YE at blueprint precursor/fragment/RT locations
  (LC carryover from PYE batch). Exclude/replace contaminated references — else recall is
  inflated and ratios distorted by real non-blueprint YE.
- **Human anchor:** compute from high-confidence, consistently-quantified human precursors
  only; report BOTH anchored and unanchored YE ratios; do not silently stack DIA-NN
  cross-run normalization + post-hoc anchor (anchor could erase the run-level effect Stage 2
  is meant to expose).
- **FDR labeling:** "conditional simulated-YE recall / empirical FDR", NOT global run FDR;
  report DIA-NN q-value separately; consider entrapment/shuffled-YE precursors for errors
  outside the truth set.
- **CRITICAL SCOPE LIMIT:** superposition adds signal AFTER acquisition physics → it does
  NOT reproduce ion suppression, source/fragmentation competition, fill-time, space charge,
  detector saturation. Stage 2 tests whether **dense real background + DIA-NN processing**
  reverse the low-abundance bias — NOT the full real-plasma phenomenon. Claim exactly this.
- Pre-register abundance bins + metrics; hold calibration runs out of evaluation; validate
  against local MS1/MS2 interference distributions vs real PYE, not just ID counts/totals.

## Questions for review
1. Is the mass-fraction-relative-to-human calibration the right model, or should YE be
   calibrated some other way (e.g. match a target ID count, or absolute molar spike)?
2. Report-quant → frame-count bridge: is matching the frame-intensity sum sound, or is
   there a better way to land the YE counts at the right level (DIA-NN normalization,
   MS1 vs MS2, the fact that reported quantity ≠ summed frame counts)?
3. A-over-080 / B-over-081: using two *different* real plasma injections as the A/B human
   backgrounds adds real run-to-run human variation. Since human is unscored
   (background_unknown) and the YE ratio is seed-baked, is this OK, or does it confound
   anything (e.g. normalization, MaxLFQ across runs, the human anchor used by the ratio
   panels)?
4. Is the human anchor still valid? The P1 ratio panels human-anchor on the observed human
   median. In Stage 2 human is REAL plasma (high dynamic range, different per run). Does
   anchoring on it still make sense, or should Stage-2 ratios anchor differently?
5. Any double-counting / contamination risk: the reference plasma IS the plasma used for
   the PYE mix — could it already contain trace yeast/E.coli that would pollute the
   "human-only background" assumption and the YE truth?
6. Validity of ground-truth FDR/sensitivity when human has no truth: P5 restricts to YE
   and treats human as background_unknown (excluded). Is excluding all human findings the
   right move, or does it bias the FDR (e.g. human false IDs that should count)?
7. Anything missing for the Stage-2 result to be scientifically defensible / comparable
   to the real-tTOF observation we're trying to reproduce?
