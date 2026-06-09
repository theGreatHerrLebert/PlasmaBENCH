# PlasmaBENCH Stage 1 — simulated PYE1 DIA benchmark

Source: [`plasmabench-stage1/`](plasmabench-stage1/) (run definitions in [`../simulations/stage1/`](../simulations/stage1/))

## Problem

We want to compare DIA search engines (DIA-NN, FragPipe, …) on TimsTOF
mixed-proteome plasma data, but real data has no exact ground truth: which
peptides are truly present, from which organism, at what between-sample ratio,
is unknown. The mixed-proteome plasma study of Ute Distler and colleagues
(*Nat Commun* 2025) establishes the spike-in design; PlasmaBENCH reproduces it
*in silico* so the truth is exact.

## Trust Strategy

**Validation via simulated ground truth.** TimSim seeds from the real Sample A
/ Sample B reports (`from_findings`) and emits a blueprint alongside each
simulated `.d`. Identification and quantification calls from each engine are
scored against (a) that blueprint and (b) the known PYE1 ratios.

This is not a substitute for real-data benchmarking — it is the layer where
exact numerical truth exists. Stage 2 (SIM overlaid on a real plasma `.d`) and
the experimental runs are the real-data complement.

## Evidence

> **Status:** the evaluator (`make stage1-eval`) is **not yet implemented** (stub,
> exits 64). The contract below is the declared target; see the scorer
> requirements note before building it.

| Aspect | Value |
|---|---|
| Oracle | TimSim peptide/precursor blueprints per sample (the final `synthetic_data.db`, **intersected** A∩B for ratios) + known B/A ratios (human 1.0, yeast 3.0, E. coli 0.5). The raw blueprint ratios carry a single global offset from TimSim's per-file median normalization — the scorer **must human-reference-normalize** (estimate one scale from human precursors, apply to all species) and report raw + normalized. |
| Tolerance | Per metric: species-resolved recall, true FDR, and human-anchored log2-ratio recovery error. Numeric thresholds set per release. |
| Command | `make stage1-eval` (scores per-tool output against the blueprint + ratio map) |
| Artifact | `results/stage1/metrics.json` plus per-tool tables under `results/stage1/<tool>/` |

Inputs:
- Sample A / Sample B reports (DIA-NN / FragPipe), provided by Ute, converted to
  TimSim seeds by `scripts/build_seed_from_report.py`.
- Background noise injected from a real TimsTOF DIA-PASEF blank acquired with the
  study method.
- Organism assignment via `plasmabench/species.py` (UniProt `_HUMAN`/`_YEAS8`/`_ECOLI`;
  the yeast spike-in is the EC1118 strain, tag `_YEAS8`, not the canonical `_YEAST`).

## Assumptions

- TimSim reproduces the TimsTOF DIA-PASEF aspects that matter for the scored
  metrics (RT, IMS, fragmentation, intensity scaling). Plasma matrix effects
  beyond injected noise are out of scope for stage 1.
- The seed reports are representative enough to make differential engine
  behaviour visible.
- Organism assignment is unambiguous for scored peptides; shared/ambiguous
  peptides are excluded from the ratio oracle.

## Failure Modes

- **Seed-engine circularity** — seeding from one engine and scoring that same
  engine inflates its self-recall. Score every engine against the union
  blueprint and read *cross-engine* differences, not absolute self-recall.
- **Simulator overfit** — tools tuned on synthetic data may improve stage 1
  without improving real-data behaviour; inspect stage 2 / experimental runs
  before release-grade claims.
- **Validation theater** — sweeps without a pre-declared decision rule give
  volume, not signal. Each sweep references a claim and a tolerance.
- **Convention gaps** — differing FDR / mass-tolerance / modification
  conventions across tools must be documented in `configs/<tool>/NOTES.md` so
  disagreement is attributed correctly.
- **Ratio leakage** — a misassigned organism pollutes two species' ratio
  distributions at once; audit the species map (the seed builder prints it)
  before scoring.
