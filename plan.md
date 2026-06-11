# Plan — `plasmabench/plots/`: LFQbench-style analysis + ground-truth panels

## Goal

Build a reusable plotting/scoring module, `plasmabench/plots/`, that reproduces the
PYE LFQ-Bench paper's analysis panels **in Python** and adds the two panels that
paper structurally cannot produce — **real FDR** and **ground-truth quant
correlation** — because TimSim emits a blueprint and they had no ground truth.

The module must run unchanged across our experimental conditions (it takes a report
path + its blueprint as inputs, hardcodes nothing condition-specific):

- **SIM type:** full SIM ⟷ YE-only SIM over real plasma (coming soon)
- **Engine/version:** DIA-NN 1.8 ⟷ DIA-NN 2.5 on the SIM data

## Two source bodies we are reconciling

1. **Ute's LFQ-Bench R scripts** — `github.com/HanYoo1402/LFQ-Bench-Scripts-for-PYE-Multicenter-Study`
   (preprint DOI 10.21203/rs.3.rs-5618718/v1; source data Zenodo 17018339;
   raw data PXD056598 / JPST003358). 12 scripts, one per figure panel, `tidyverse`/`ggplot2`.
   **No ground truth** — everything is expected-ratio / consensus based.
2. **timsim-bench** (`github.com/theGreatHerrLebert/timsim-bench`, checked out at
   `/scratch/timsim-demo/SUBMISSION/timsim-bench`) — our own Python benchmark code with
   ground-truth-aware FDR/TPR and abundance-binned sensitivity. This is the source for
   the two extra panels.

## What we already have

`scripts/plot_stage1_ratios.py` is our version of Ute's **Fig_5A** (per-species violin)
+ **Fig_5B_C** (LFQbench MA scatter). It works at **precursor level** on
`Precursor.Normalised`, computes `log2(B/A)`, human-anchors (shift so human median = 0),
and reads `EXPECTED_RATIO_B_OVER_A` from `plasmabench/species.py`.

## How Ute's Fig_5A differs from our current script (the deltas to close)

Biology matches exactly (independent confirmation `species.py` is right): same `_YEAS8`
suffix on `Protein.Names`, same "multi-species → exclude" rule, same expected ratios.
Mechanics differ:

| Axis | Ute Fig_5A | Ours (`plot_stage1_ratios.py`) | Decision |
|---|---|---|---|
| Ratio sign | `log2(A/B)` (yeast −log2 3, E.coli +1) | `log2(B/A)` (yeast +log2 3, E.coli −1) | Canonical `log2_b_over_a` + expected B/A stored internally; **one** display API derives `display_ratio = sign · canonical` (default `sign=-1` → paper's `log2(A/B)`), and the expected reference lines/labels/exports flow through the **same** transform. Tests assert yeast & E.coli signs in both conventions (see Decisions) |
| Feature unit | Protein.Group / PG.MaxLFQ | Precursor.Id / Precursor.Normalised | **protein-level (PG.MaxLFQ) primary** for paper comparability; precursor-level a secondary diagnostic (our blueprint truth is precursor-level) |
| Human anchor | none (human=0 by design) | shift human median → 0 | keep ours; document as deviation |
| FDR filter | 4-level: `Q.Value`, `Lib.Q.Value`, `PG.Q.Value`, `Lib.PG.Q.Value` ≤ 0.01 | none (assumes pre-filtered) | **port the 4-level cut** |
| Contaminants | drop `_CONTA` in Protein.Names | none | **port** |
| Replicate filter | ≥3 non-NA per condition | none (median over whatever runs match) | **port as configurable min-replicates** (no-op when 1 rep/sample) |
| MA x-axis | `log2(Mean_A)` | `½·log2(A·B)` (true MA midpoint) | keep ours (better convention) |

**Aggregation:** per-condition replicate aggregation is a **mode**,
`aggregation ∈ {linear_mean, geometric_mean, median}`. Use `linear_mean` for
paper-reproduction parity (Ute averages PG.MaxLFQ in linear space); `geometric_mean`
becomes the benchmark default only *after* a documented comparison. Do **not** claim
geometric mean is the violin "weird grouping" bug Ute's team flagged — Codex review:
insufficient evidence; more plausibly grouping keys / protein-group collapse / NA
handling / mixed runs. Keep that as an open note, not a fix.

## Target panels (port = Ute panel; new = ground-truth-only)

| # | Panel | Source | Notes |
|---|---|---|---|
| P1 | per-species ratio violin + MA scatter | port Fig_5A / Fig_5B_C (already have; apply deltas above) | precursor **and** protein level |
| P2 | ID completeness / counts per species | port Fig_2A / Fig_2B | feeds recall; per-species protein-group counts + run-completeness |
| P3 | missingness / LOD sigmoid | port Fig_5D_E **but** replace its PYE_B-as-proxy-abundance with **real blueprint intensity bins** via adapted `get_sensitivities` (column-select before any `dropna`) — `simulated` species only; all of yeast+E.coli, not just yeast |
| P4 | lower-tertile ratio violin | port Fig_6C_D | minor refinement of P1 |
| **P5** | **simulated-spike-in FDR / TPR** | timsim-bench `fpr_tpr_plots.calculate_fdr` / `calculate_tpr` | **no Ute equivalent** — found-set vs TimSim truth-set, scoped per `(experiment, sample)`, restricted to `simulated` species. Protein scoring: **report the reported group as one unit** — any normalized member ∩ truth = group-TP (primary), strict all-members-in-truth = reported bound. Do **not** use the current member-expansion rule (it dumps non-truth candidates into TP) |
| **P6** | **ground-truth quant correlation** | timsim-bench (simulated `intensity` vs recovered quantity, per species) | **no Ute equivalent** |

## Data model — TWO normalized tables, not one merged table

Codex review (Critical): a single concatenated table that overloads `intensity`
(truth vs observed) by `software` is error-prone and missing columns several panels need.
**Core abstraction = two normalized tables; concatenation is only an interchange format.**

**`truth`** (from the TimSim blueprint, per SIM `.d`):
`experiment, sample, run_id, sequence, sequence_modified, charge, precursor_id,
sim_protein_id, protein_id, protein_names, species, truth_intensity, transmitted,
truth_scope`.

> **Blueprint source contract (verified against a real `synthetic_data.db`).** The
> blueprint is the `synthetic_data.db` SQLite **beside each `.d`** — there is **one DB per
> run** and it has **no run/sample/experiment column** (`frames`/`scans` are acquisition
> geometry only). So `experiment`/`sample`/`run_id` come from the **run manifest** (below),
> not the DB. Within a DB:
> - identity: `peptides ⋈ ions on peptide_id`; precursor key `(sequence_modified, charge)`.
> - `peptides.sequence` is **already UniMod-notated** → it maps to `sequence_modified`;
>   `sequence` (stripped) is derived.
> - `truth_intensity` = `peptides.events × ions.relative_abundance` (events carries the
>   dilution; relative_abundance splits a peptide across charges, sums to 1).
> - **`transmitted` flag (REQUIRED).** `peptides ⋈ ions` over-counts: not every ion reaches
>   `fragment_ions` (verified: 40,184 ions vs 39,815 in `fragment_ions` for sample A → 369
>   never transmitted). Untransmitted precursors are **unfindable**, so counting them makes
>   them false negatives in P3/P5. Add `transmitted = ion ∈ fragment_ions`; **transmitted
>   truth is the primary identification denominator** (timsim-bench intersects fragment_ions
>   before scoring). Optionally report all-generated sensitivity as a secondary curve.
> - **protein id split.** `peptides.protein_id` is an **internal integer** → keep as
>   `sim_protein_id`. `protein_id` (for P5 joins vs DIA-NN accessions) is a **normalized
>   accession parsed from `peptides.protein`** (the FASTA header, e.g.
>   `sp|Q08211|DHX9_HUMAN` → `Q08211`); header parsing must be tested. `protein_names` =
>   the full header → feeds `assign_species`. `peptides.decoy` available for exclusion.

### Run manifest (the join authority)
Because the blueprint carries no run id and DIA-NN `Run` stems don't match `.d`/dir names,
**a manifest is the single source of truth**: it maps each physical `.d` + its
`synthetic_data.db` → canonical `{experiment, sample, run_id}`. Both `truth` and
`observations` resolve `(experiment, sample, run_id)` **through the manifest**, which
**fails on any missing or duplicate mapping**. Downstream: **P5 = scoped set comparison**
(found-set vs truth-set per `(experiment, sample[, run])`), **never an inner join**; **P6 =
explicitly validated many-to-one precursor join**.

**`observations`** (from the engine report, per condition):
`experiment, sample, run_id/filename, software, sequence, sequence_modified, charge,
precursor_id, protein_id, protein_names, protein_group, species, truth_scope,
observed_intensity, pg_maxlfq, q_value, lib_q_value, pg_q_value, lib_pg_q_value` (the four
q-values gate the LFQ-Bench filter; carry whichever the engine emits). `truth_scope` is
joined in from the truth/condition (so human observations are tagged `background_unknown`
for filtering).

Panels are **panel-specific joins** of these two on `(experiment, sample[, run_id])`.
`table_processing.process_diann_table` is the starting point for `observations` but does
**not** currently carry `run_id`, `Lib.Q.Value`, `Lib.PG.Q.Value`, or a `species` column —
the adapted loader must add them.

### Truth coverage / `truth_scope` (resolves the YE-only-real-plasma condition)
Every truth/observation row carries `truth_scope`:
- `simulated` — yeast & E. coli (blueprint has exact truth).
- `background_unknown` — human in the YE-only-SIM-over-real-plasma condition (real plasma,
  **no** blueprint truth). In the full-SIM condition, human is `simulated`.

Consequences, enforced in the metric layer:
- **P5 is "simulated-spike-in FDR/TPR"** restricted to rows confidently assigned to
  yeast/E. coli. Do **not** label it "real FDR" in the YE-only condition. In full-SIM it
  is genuine global FDR.
- **P3 sensitivity** computed only over `simulated` species; human sensitivity is
  *unavailable* in the YE-only condition (state this).
- Human (`background_unknown`) rows **still feed P1/P2/P4** as empirical observations.
- Cross-species / ambiguous findings → **`unscorable`, never FP** (matches `species.py`
  multi-species→None). Documented limitation: an error mis-assigned to human escapes the
  spike-in FDR.

### Three integration fixes when lifting timsim-bench code
1. **Species labeling.** `table_labeling.py` hardcodes `.contains("YEAST")`. Our yeast is
   EC1118 = `_YEAS8`. **All species assignment routes through `plasmabench/species.py`**
   (handles `_YEAS8`, multi-species→None), populating the `species` column. `protein_id`
   alone is bare accessions — assignment needs `protein_names` (suffix-bearing).
2. **Sample A/B detection.** `quant_utils.load_and_prep` infers A/B from `001.d`; our SIM
   filenames use `sampleA`/`sampleB`. Parameterize the A/B rule.
3. **Scope truth-set by `(experiment, sample[, run])`.** The lifted `analyze_datasets`
   forms ONE TimSim set across all supplied data — if A and B blueprints differ, a find in
   A gets credited for existing in B. FDR/TPR/sensitivity/P6 joins must scope per
   `(experiment, sample[, run])` before combining summary statistics.
4. **`get_sensitivities` `dropna()` trap.** It calls unrestricted `dropna()` — the exact
   silent-data-loss trap our own findings contract warns about (one null optional column
   deletes rows). **Select required columns first**, and set the truth key explicitly
   (precursor truth → `(sequence_modified, charge)`, not unmod sequence+charge).

## Proposed module layout

```
plasmabench/plots/
  __init__.py
  loader.py        # report+blueprint → canonical merged long table (the contract above)
  ratios.py        # P1, P4  (violin + MA + tertile), log2(A/B) display
  identifications.py # P2
  sensitivity.py   # P3  (wraps get_sensitivities with blueprint intensity)
  fdr.py           # P5  (wraps calculate_fdr/tpr + group-credit)
  quant.py         # P6
  style.py         # shared house style (lift FONT_* + species colors from existing script)
scripts/
  plot_stage1_ratios.py  # becomes a thin CLI over plasmabench.plots.ratios (back-compat)
```

## timsim-bench reuse strategy — DECISION NEEDED

timsim-bench is a research repo (notebooks, not cleanly packaged) and is **not** currently
a PlasmaBENCH submodule (we have `rustims`, `evident`). Options:
- **(A) Vendor** the small function set (`calculate_fdr`, `calculate_tpr`,
  group-credit protein rule, `create_sets`, `get_sensitivities`, `process_diann_table`)
  into `plasmabench/plots/`, keep-in-sync. **Matches existing precedent** — CLAUDE.md
  already notes `peptides.py` is "kept in sync with timsim-bench/benchmark/helpers.py".
- (B) Add timsim-bench as a third submodule and import.
- (C) pip-install timsim-bench.

**DECISION: (A) vendor as an *adapted port*** (Codex). A submodule would drag in
timsim-bench's notebook/research deps and unstable APIs. We must edit the lifted code
anyway (species, A/B, scoping, `dropna`, group-credit). Requirements: **record the source
commit** of each ported function in a header comment, add **focused parity tests** against
the original on a fixed fixture, and call it an *adapted port* — not vaguely "keep in
sync." Update the CLAUDE.md note that currently describes peptides.py that way.

## Replicate generation — controlled-jitter precision (resolves open Q1)

Stage 1 sims one `.d` per sample, so the paper's CV/precision panels (Fig_4) and the
≥N-replicate filter looked impossible. They're not — TimSim's **`from_existing`** route is
built for exactly this. Recipe:

1. **Base sim** = the normal Stage-1 `from_findings` run → SIM `.d` + its synthetic-
   experiment DB (peptides/ions/proteins; carries `retention_time_*`, `inv_mobility_*`,
   `events`).
2. **Each replicate** = a `from_existing` run with `existing_path` = the base DB and the
   three jitter knobs set:
   - `rt_variation_std` (seconds) → Gaussian on the RT column (`add_normal_noise`)
   - `ion_mobility_variation_std` → Gaussian on the IM column (`add_normal_noise`)
   - `intensity_variation_std` → **log1p-space** Gaussian on `events`
     (`add_log_noise_variation`); ≈ fractional CV at high abundance, compressed at low.
   **Identity and species are identical** across replicates; **truth intensity is NOT** —
   it is deliberately re-jittered each replicate (that's the point). Each replicate's own
   `synthetic_data.db` records its own jittered `events`, i.e. **per-replicate truth**.

**What this measures (and doesn't).** Replicate CV here reflects *our injected* dispersion,
not instrument/biological variance — a **controlled precision floor**, NOT the paper's
real-world CV (**label it as such, never as instrument CV**). Two subtleties from review:
- **Configured sigma ≠ injected CV.** `add_log_noise_variation` is Gaussian in **log1p**
  space with a **hard clip at the base run's max** → linear-space noise is not mean-zero,
  the most-abundant peptide can only move *down*, and low intensities don't follow the
  high-abundance lognormal approx. So **P7 plots recovered CV vs the *realized* truth CV
  computed from each replicate's blueprint**, with configured sigma as the experimental
  setting only — not as the x-axis truth.
- **A/B dilution is preserved only *in distribution*, not exactly per replicate.** Fine for
  the quant oracle because each replicate is scored against *its own* blueprint truth; the
  expected A/B ratio holds across replicates, not within one.

**Dilution / A/B ratio is carried in `events` — guard it on two fronts.** The spike-in
ratio lives in the `events` column: `from_findings` sets `events ∝ seed intensity`
(`load_findings.py:484`), the base sim persists it, and `from_existing` loads `events`
(`simulator.py:906`) and only re-jitters it (`:964`).
- **Config guards (render-configs must assert all):** `from_existing=True`,
  `from_findings=False`, **`proteome_mix=False`** (else the dilution loop
  `simulator.py:970–979` overwrites `events = total_events × dilution_factor`),
  **`phospho_mode=False`** (it mutates peptide/ion *identity* post-load,
  `simulator.py:981`), **`re_scale_rt=False`**.
- **Post-generation invariant assertions (catch future TimSim changes config checks
  miss):** each replicate's blueprint must have **identical modified precursor keys** and
  **unchanged species membership** vs the base, and we **record the realized event ratios**
  rather than assume them.

Caveats to honor in the runner:
- `from_existing` and `from_findings` are **mutually exclusive** (simulator raises) — base
  is `from_findings`, replicates are `from_existing`. Two-step, not one run.
- Keep `re_scale_rt` **off** for replicates (it re-predicts deterministic RT → no jitter);
  jitter must come from `rt_variation_std`.
- Noise draws from the **unseeded global RNG** (`sample_seed=41` only seeds peptide
  subsampling, not noise) → separate replicate processes differ for free, but exact jitter
  is **not reproducible**. If we need reproducible replicates, that's a small upstream
  rustims enhancement (add a noise seed) — flag, don't block on it.

## In scope, reframed
- **P7 (new): controlled-jitter precision** — port Fig_4A_C/4B_D (RT-CV / protein-CV
  boxplots) over `from_existing` replicates, plotted as recovered-CV vs injected-CV. The
  one panel set that genuinely needs the replicate machinery above.

## Out of scope (for now)
- Fig_3C_D / Fig_3E_F (PYE1-vs-PYE3/9 dilution series) — not our PYE1 A/B design.
- The scorer's numeric report tables (`stage1-eval-spec`) — this plan is the **plotting/
  metric layer**; the scorer can consume the same `loader.py` tables later.

## Decisions (resolved in Codex review, 2026-06-09)

- **Data model:** two normalized tables (`observations` + `truth`), panel-specific joins;
  concatenation only as interchange. (was Q1/Q2)
- **Sign:** single display API `display_ratio = sign · canonical_log2_b_over_a`; expected
  lines/labels/exports flow through it; tests assert both conventions. (Q1)
- **Aggregation:** mode `{linear_mean, geometric_mean, median}`; linear for paper-repro;
  geometric not claimed as the R bug. (Q2)
- **timsim-bench reuse:** vendor as *adapted port* — record source commit, add parity
  tests. (Q3)
- **Primary unit:** **protein-level (PG.MaxLFQ) primary**, precursor-level secondary
  diagnostic. (Q4)
- **Protein FDR scoring:** any-member-∩-truth = group-TP primary, strict = reported bound;
  no member-expansion. (Q5)
- **YE-only-real-plasma:** `truth_scope` tagging; truth metrics restricted to `simulated`
  species; P5 labelled "simulated-spike-in FDR"; human enters P1/P2/P4 only; cross-species
  → `unscorable`. (Q6)

## Decisions (resolved in Codex review v2, 2026-06-09)
- **P6 quant-correlation metric:** primary **Spearman(log truth, log observed)** per
  `(run, species)`; secondary **log-Pearson after a documented per-run scale
  normalization**. (Spearman also sidesteps the `events`-vs-reported-quantity unit
  mismatch.) Decide once here so scores are comparable across panels.
- **Run identity:** explicit **manifest** (not path-inference) — see "Run manifest" above.
- **Truth universe:** `transmitted` truth (∩ `fragment_ions`) is the primary denominator.
- **Replicate truth:** per-replicate (re-jittered), scored against each replicate's own
  blueprint; P7 x-axis = realized truth CV, not configured sigma.

## Smallest shippable slice
`loader.py` + **precursor-level P6** on one A/B report, with the full contract enforced:
1. **manifest-based run resolution** (fail on missing/duplicate);
2. **modified-sequence canonicalization** (blueprint UniMod ↔ engine notation via
   `plasmabench/peptides.py`);
3. **decoy exclusion** (`peptides.decoy`);
4. fields: `transmitted`, `sim_protein_id`, normalized-accession `protein_id`;
5. **uniqueness assertion** on `(experiment, sample, run_id, sequence_modified, charge)` for
   truth; **explicit duplicate-collapse policy** for observations;
6. **tests against the existing Sample-A blueprint + its DIA-NN parquet** (real fixtures).

Then layer P1 (protein-level ratio violin) as the second panel. All other panels are
pure functions of the two tables.
```
