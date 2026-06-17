# Reproducing PlasmaBENCH — step-by-step guide

This is the end-to-end *Anleitung* for reproducing the PlasmaBENCH results from
the archived data. It covers three depths of reproduction, from cheapest to most
complete:

1. **Replay the figures & metrics** from our published DIA-NN reports — no
   search engine needed (~minutes).
2. **Re-search the simulated `.d`s** with your *own* DIA-NN / FragPipe, then
   re-score (~hours, needs the engines).
3. **Regenerate the simulations** from the seeds with TimSim, then search and
   score (~a day, needs a GPU + the TimSim build).

> **You must bring your own search engines.** DIA-NN and FragPipe are *not*
> redistributed here — they have their own licenses. Download/install them
> yourself (DIA-NN 1.8.1, 2.5.0, 2.6.0 are the three version-axis points;
> FragPipe for the cross-engine arm). `scripts/get_diann.sh` helps fetch DIA-NN
> builds you are licensed to use; FragPipe install is manual. None of the
> archived data contains engine binaries or predicted spectral libraries.

---

## 0. Prerequisites

| Need | For which tier | Notes |
|---|---|---|
| `git` (with submodules) | all | clone recursively |
| Python 3.12 | all | analysis venv |
| ~50 GB free disk | tiers 1–2 | the archive unpacks to ~44 GB |
| DIA-NN 1.8.1 / 2.5.0 / 2.6.0 | tiers 2–3 | **your own** licensed binaries |
| FragPipe | tier 2 (cross-engine) | **your own** install |
| NVIDIA GPU + CUDA | tier 3 only | TimSim signal rendering |
| Docker (optional) | any | hermetic analysis/replay image |

---

## 1. Get the code

```bash
git clone --recurse-submodules https://github.com/theGreatHerrLebert/PlasmaBENCH.git
cd PlasmaBENCH
# or, if you already cloned without --recurse-submodules:
make submodules        # = git submodule update --init --recursive
```

The submodules are `rustims/` (the TimSim engine, pinned commit) and `evident/`
(the claim-manifest validator + trust patterns). Both are required.

---

## 2. Get the data (from Zenodo)

The raw inputs, simulated `.d`s, blueprints, and DIA-NN reports are **not in git**
(they are 10s of GB and gitignored). Download them from the Zenodo record:

> **Zenodo DOI:** [`10.5281/zenodo.20733913`](https://doi.org/10.5281/zenodo.20733913) (concept DOI — always resolves to the latest version)

The archive ships as three tarballs that unpack **into the repo root**, recreating
the exact directory layout the scripts expect:

```bash
# from the repo root, after downloading the .tar.zst files into ./_zenodo_dl/
zstd -d --stdout _zenodo_dl/plasmabench-data-raw.tar.zst    | tar -xf - 
zstd -d --stdout _zenodo_dl/plasmabench-simulations.tar.zst | tar -xf -
zstd -d --stdout _zenodo_dl/plasmabench-results.tar.zst     | tar -xf -

# verify integrity against the shipped checksums
sha256sum -c _zenodo_dl/SHA256SUMS
```

What each tarball contains:

| Tarball | Unpacks to | Contents |
|---|---|---|
| `plasmabench-data-raw` | `data/raw/` | search FASTA, the blank `.d` (noise source), the real G PYE1 `.d`s, FragPipe DDA inputs |
| `plasmabench-simulations` | `simulations/` | simulated Sample A/B `.d`s, **`synthetic_data.db` (the TimSim blueprint truth)**, seeds, rendered configs |
| `plasmabench-results` | `results/` | our DIA-NN 1.8 / 2.5 / 2.6 reports (`report.tsv`/`.parquet`, stats, `.quant`) — **predicted spectral libraries (`*.predicted.speclib`, ~20 GB) are deliberately excluded** |

After unpacking, the tree matches what `CLAUDE.md` and `evident.yaml` describe.

---

## 3. Set up the environment

Pick one:

```bash
# (a) Analysis only — enough for tier 1 and scoring. No TimSim, no GPU.
make venv && source .venv/bin/activate

# (b) From source with TimSim — needed for tier 3 (regenerating simulations).
#     Builds imspy_connector + imspy-* + timsim from the pinned rustims submodule.
make venv-source && source .venv/bin/activate

# (c) Docker analysis image — hermetic replay, code pinned, data mounted.
#     (used in tier 1 below; nothing to install locally)
```

---

## 4. Tier 1 — replay figures & metrics (no search engine)

Each EVIDENT claim has one turnkey `make` target. With the data in place:

```bash
make stage1-eval        # Stage-1 species-resolved recall / FDR / ratio recovery
make stage2-eval        # Stage-2 spike-in ratios on real plasma (sim-QC gate)
make quant-eval         # DIA-NN cross-run normalization → ratio IQR inflation
make real-vs-sim-eval   # real-vs-SIM ratio compression capstone
make lowinput-eval      # low-input per-member under-quantification
make idoverlap-eval     # cross-engine real-data ID overlap (1.8 vs 2.5 / 2.6)
```

Outputs land under `results/figures/`, `results/stage1/`, `results/stage2/`.
The mapping of claim → command → artifact is below (§7).

**Hermetic alternative** (runs the same target inside the pinned analysis image
with your data bind-mounted):

```bash
make eval-docker EVAL=lowinput-eval
# override data roots if they live elsewhere:
PLB_OUT=/path/to/simulations PLB_DATA=/path/to/data make eval-docker EVAL=stage1-eval
```

Validate the claim manifest itself:

```bash
python evident/workflow/validate_manifest.py evident.yaml
```

---

## 5. Tier 2 — re-search the simulated `.d`s with your own engines

This re-derives the DIA-NN reports under `results/` instead of trusting ours.

**DIA-NN** (the version-comparison axis — repeat per version you have):

```bash
# args: <version> <fasta> <out-dir> <.d ...>
bash configs/dia-nn/run.sh 2.5 \
  data/raw/fasta/plasma_yeast_ecoli.fasta \
  results/diann-2.5-AB \
  simulations/stage1/dia/sampleA_imposed/*/*.d \
  simulations/stage1/dia/sampleB_imposed/*/*.d
```

The unified search config lives in `configs/dia-nn/` (only the binary differs
across versions — that is the whole point of the comparison). Document any
untunable knob in `configs/<tool>/NOTES.md`.

**FragPipe** (cross-engine arm): run your install with the workflow in
`configs/fragpipe/` (see its `NOTES.md`) over the same simulated `.d`s.

Then re-run the tier-1 `make *-eval` targets, which now read *your* reports.

> **Cross-engine integrity:** scoring an engine on data seeded from its own
> report inflates self-recall. Read *cross-engine differences* against the union
> blueprint, never absolute self-recall. The three-species design is the control
> (a yeast peptide quantified at the *E. coli* ratio is a detectable error). See
> `evident.yaml` → `failure_modes`.

---

## 6. Tier 3 — regenerate the simulations from scratch (needs TimSim + GPU)

Only if you want to rebuild the simulated `.d`s rather than download them.
Requires the `make venv-source` environment (§3b).

```bash
# 6a. Build a seed from a DIA-NN/FragPipe report (one per sample; --run-prefix
#     splits a combined A+B report — without it the A/B ratio is destroyed).
python scripts/build_seed_from_report.py --engine diann \
  --report <A+B report.tsv> --sample A --run-prefix A1 \
  --ratio-mode imposed --rt-max 1480 \
  --out simulations/stage1/seeds/seed_sampleA_imposed.csv
# ...repeat for sample B.

# 6b. Render configs (fills findings_path/save_path/experiment_name).
python scripts/render_configs.py \
  --base simulations/stage1/configs/base-dia.toml \
  --seeds-dir simulations/stage1/seeds --out-dir simulations/stage1/dia \
  --samples A B --seed-suffix _imposed --path-mode host

# 6c. Simulate (emits the .d + synthetic_data.db blueprint).
timsim simulations/stage1/dia/sampleA_imposed/config.toml
timsim simulations/stage1/dia/sampleB_imposed/config.toml

# 6d. then re-search (§5) and re-score (§4).
```

The whole Stage-1 flow is also wrapped:

```bash
ENGINE=diann REPORT_A=.../A/report.tsv REPORT_B=.../B/report.tsv bash scripts/run_stage1.sh
DRY_RUN=1 bash scripts/run_stage1.sh    # plan only
```

See `CLAUDE.md` (“Stage-1 data flow”, “the findings contract”) and
`simulations/stage1/README.md` for the scientific rationale and the two
silent-data-loss traps the seed builder defends against.

---

## 7. Claim → command → artifact

| EVIDENT claim (`evident.yaml`) | Command | Artifact |
|---|---|---|
| `plasmabench-stage1-simulated-truth` | `make stage1-eval` | `results/stage1/metrics.json` |
| `plasmabench-stage2-superimpose-truth` | `make stage2-eval` | `results/stage2/metrics.json` |
| `plasmabench-diann-crossrun-normalization` | `make quant-eval` | `results/figures/quant_18v25_normalization.png` |
| `plasmabench-real-vs-sim-validation` | `make real-vs-sim-eval` | `results/figures/real_vs_sim_compression.png` |
| `plasmabench-lowinput-overseparation` | `make lowinput-eval` | `results/figures/sim_lowabund_underquant.png` |
| `plasmabench-ideoverlap-real` | `make idoverlap-eval` | `results/figures/real_overlap_18_vs_25.png` |

Each claim’s `evidence.command` field in `evident.yaml` is the authoritative
source for its replay command.

---

## 8. Troubleshooting

- **`make stage1-eval` exits 64** — the `help` text still describes the old stub;
  the real recipe runs `scripts/stage1_eval.py`. If it errors on missing inputs,
  confirm the `simulations/` and `results/` tarballs are unpacked.
- **DIA-NN reports 0 IDs on the SIM** — known gotcha: DIA-NN auto mass-accuracy
  can pick `--mass-acc 0`. Use a fixed 15 ppm (the config does). The simulation
  is fine; do **not** re-simulate.
- **Yeast peptides missing from the FASTA** — use the EC1118 strain proteome
  (organism 643680), the one Ute uses; entry names are `_YEAS8`, not `_YEAST`.
- **Ratios look halved (~0.7×)** — the historical `from_findings` per-sample
  median normalization bug. Fixed upstream (rustims PR #407, shared
  `reference_median`); use the `dia-fix`/`dia-refmed` simulation trees.

---

## 9. Citation & provenance

PlasmaBENCH simulates the mixed-proteome plasma benchmark of Distler et al.,
*Nat Commun* 2025, DOI [10.1038/s41467-025-64501-z](https://doi.org/10.1038/s41467-025-64501-z)
(data PXD056598). The simulator is TimSim from the
[rustims](https://github.com/theGreatHerrLebert/rustims) project. See
`zenodo/ARCHIVE_README.md` for the data-record citation and the exact pinned
versions.
