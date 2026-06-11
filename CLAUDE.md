# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PlasmaBENCH benchmarks DIA software (DIA-NN, FragPipe, …) on TimsTOF **mixed-proteome plasma** data, using **TimSim**-generated synthetic ground truth. The design ("PYE1": Plasma · Yeast · E. coli) is a three-species ratio test à la LFQbench/HYE: human plasma background stays 1:1 between Sample A and B, while yeast (2%→6%, expected B/A = 3.0) and E. coli (8%→4%, expected B/A = 0.5) move in opposite directions. Because TimSim emits a blueprint alongside every simulated `.d`, species-resolved recall, FDR, and ratio-recovery can be scored against exact truth.

Read `README.md` and `simulations/stage1/README.md` first — they carry the scientific rationale this file does not repeat.

## The one architectural decision that governs everything: `from_findings`, not `proteome_mix`

TimSim has two simulation routes. **PlasmaBENCH uses `from_findings` exclusively** for Stage 1:

- **`from_findings`** (this repo): peptides come from Ute's *real* DIA-NN/FragPipe reports for Sample A and B. The A/B spike-in ratios are **baked into the report intensities + organism membership** of the seed CSV. Config: `from_findings = true`, `[proteome_mix] proteome_mix = false`, and **no `[digestion]` block** (peptides come from the seed, not a FASTA digest).
- **`proteome_mix` + dilution CSV** (NOT used here): in-silico FASTA digest with ratios from a `dilution_factors.csv`. This is the fully-synthetic route. **Never enable `proteome_mix` on this path** — it would override the baked-in ratios.

The `fasta_path` in configs is *ignored* by `from_findings` but kept valid because the downstream search-engine step (DIA-NN/FragPipe) needs it.

### The findings contract (and its two silent-data-loss traps)

The seed CSV has 6 columns: `protein, sequence, intensity` (required) + `charge, rt, im` (optional). TimSim's `load_findings` feeds this to a Rust parser. Two failure modes that `build_seed_from_report.py` defends against — preserve these invariants in any edit:

1. **Optional columns must be fully populated or omitted entirely — never partially NaN/blank.** TimSim runs `dropna()` over *every present column*, so one all-blank `rt`/`im` column deletes every row. The builder only writes an optional column if `.notna().all()`.
2. **`sequence` must parse in the Rust `PeptideSequence` parser or it *aborts the whole process*.** Sequences must be UniMod-notated (`PEP[UNIMOD:4]TIDE`, optional leading N-term tag). The builder validates with `is_valid_timsim_sequence` and drops bad rows *before* they reach Rust.

## Python package (`plasmabench/`)

Not pip-installable; scripts add the repo root to `sys.path`. Four modules, each load-bearing for correctness:

- **`species.py`** — organism assignment by UniProt entry-name suffix. Tags live in `Protein.Names` (`ALBU_HUMAN`), **not** `Protein.Group` (bare accessions). Yeast appears as `_YEAS8` (strain-specific, the one Ute uses) — **not** `_YEAST`; both map to YEAST. A peptide matching >1 organism returns `None` and is **excluded from the ratio oracle** (cross-species leakage is the worst failure mode here). `PYE1_COMPOSITION` and `EXPECTED_RATIO_B_OVER_A` are the quant oracle.
- **`peptides.py`** — engine notation → UniMod. DIA-NN: `(UniMod:n)` → `[UNIMOD:n]`. FragPipe: stripped `Peptide` + separate `Assigned Modifications` column, mass→UniMod via sagepy's `mass_to_mod` (lazily imported so the DIA-NN path needs no sagepy). Kept in sync with `timsim-bench/benchmark/helpers.py`.
- **`intensities.py`** — only used in `--ratio-mode imposed`. Replaces reported (ratio-compressed) intensities with a Tenzer rank-abundance baseline × per-species dilution. Baseline is keyed on a **`crc32(sequence|charge)`** seed (NOT Python's salted `hash()`) so Sample A and B assign the *same* baseline to a shared peptide → exact fold change. Ported verbatim from rustims `simulate_proteins.py`; keep it in sync.
- **`utility.py`** — generic table/sequence helpers.

### observed vs imposed ratio mode

`build_seed_from_report.py --ratio-mode`:
- **`observed`** (default): use the report's real (DIA-NN-compressed) intensities; oracle = the seeded ratio (not nominal — DIA-NN compresses ratios, so the truth is what was actually seeded).
- **`imposed`**: discard found intensities, assign Tenzer-baseline × species-fraction so A/B is exactly nominal (1.0/3.0/0.5) — a clean accuracy oracle. The current Stage-1 runs (`results/diann-2.5-AB*`) use imposed seeds (`seed_sample*_imposed.csv`).

## Stage-1 data flow

```
DIA-NN/FragPipe report (combined A+B, long-format)
  → scripts/build_seed_from_report.py  (one run per sample, --run-prefix/--run-contains to split conditions)
  → simulations/stage1/seeds/seed_sample{A,B}[_imposed].csv
  → scripts/render_configs.py  (fills findings_path/save_path/experiment_name into base-dia.toml)
  → simulations/stage1/dia/sample{A,B}/config.toml
  → timsim <config.toml>  → SIM .d + blueprint
  → configs/dia-nn/run.sh  (search the SIM .d's)  → results/
```

Reports are **long-format with both conditions × 6 replicates concatenated**. `--run-prefix A1`/`B1` (G-site) or `--run-contains Slot2-11`/`Slot2-12` (L-site) selects one condition; without it A and B merge and the ratio is destroyed. Per-precursor intensity = **median across replicates** (not sum — sum biases toward frequently-detected peptides). `protein` aggregates as the *union* of protein strings (so cross-species ambiguity is detected, not order-dependent).

### The Stage-1 recipe (current, locked)

G-site found IDs simulated over the **L-method DIA blank** `G250513_005_Slot2-2_1_17993.d` (the only blank available). The G IDs elute to ~2121s (35 min) but the L blank spans only ~1500s (25 min), so seeds are RT-rescaled (`build_seed_from_report.py --rt-max 1480`) to fit. RT is not a scored quantity, so uniform compression is safe. The blueprint truth is exact regardless of this mismatch.

## Common commands

```bash
# Environment — timsim is NOT on PyPI; build from the pinned rustims submodule.
make venv-source        # builds imspy_connector (maturin) + imspy-* + timsim into .venv
source .venv/bin/activate
make venv               # lightweight: dev tools only, NO timsim (seed building, analysis)

# Docker alternative (once the from-findings image is published to GHCR):
make pull && make shell # GPU shell with timsim on PATH; mounts /work /data /sim
make timsim-help

# Stage 1 end-to-end:
ENGINE=diann REPORT_A=.../A/report.tsv REPORT_B=.../B/report.tsv bash scripts/run_stage1.sh
DRY_RUN=1 bash scripts/run_stage1.sh                          # plan only

# Build one seed manually (note --run-prefix to split a combined report):
python scripts/build_seed_from_report.py --engine diann \
  --report .../report.tsv --sample A --run-prefix A1 \
  --ratio-mode imposed --rt-max 1480 \
  --out simulations/stage1/seeds/seed_sampleA_imposed.csv

# Render configs from a host venv (rewrites container /work,/data,/sim → host paths):
python scripts/render_configs.py --base simulations/stage1/configs/base-dia.toml \
  --seeds-dir simulations/stage1/seeds --out-dir simulations/stage1/dia \
  --samples A B --seed-suffix _imposed --path-mode host

# Search the SIM .d's (two DIA-NN versions = the version-comparison axis):
bash configs/dia-nn/run.sh 2.5 data/raw/fasta/plasma_yeast_ecoli.fasta results/diann-2.5 \
  simulations/stage1/dia/sampleA_imposed/*.d simulations/stage1/dia/sampleB_imposed/*.d

# Validate the EVIDENT claim manifest:
python evident/workflow/validate_manifest.py evident.yaml

# Run a module doctest (peptides/species carry doctests):
python -m doctest plasmabench/peptides.py -v
```

There is **no pytest suite and no linter config** — the tests that exist are doctests in `plasmabench/peptides.py` and `plasmabench/species.py`. `make stage1-eval` is a deliberate stub that exits 64; the scorer (species-resolved ratio recovery, FDR, recall vs blueprint) is **not yet implemented** — it's the next major piece of work (see `TODO.md`, memory `stage1-eval-spec`).

## Paths and portability

Configs reference container-fixed roots so they stay machine-portable: `/work`=repo, `/data`=`data/` (override `PLB_DATA`), `/sim`=`simulations/` (override `PLB_OUT`). `render_configs.py --path-mode host` rewrites these to host paths for from-source venv runs (where `/work` etc. don't exist). Raw inputs live under `data/raw/{dia/blanks,dia/plasma,fasta}/` and are gitignored.

## Submodules

- **`rustims/`** — tracks branch `main` (the from_findings 6-column contract **plus** the shared `reference_median` fix are now merged to main via PR #407; `rustims:latest` still differs, so PlasmaBENCH owns its own image). Currently pinned at `2748248`. Bumping the pin: see `docs/venv-setup.md`. The Makefile `IMAGE_DIGEST` is built+pushed automatically by `.github/workflows/docker-publish.yml` when a `rustims`-pointer change lands on `main`; read the pushed `sha256` from that run into `IMAGE_DIGEST` (until then it stays `UNPINNED` and `make pull-pinned` errors by design).
- **`evident/`** — read-only reference; provides the claim-manifest validator and trust-pattern docs.

Clone with `git clone --recurse-submodules` or `git submodule update --init --recursive` (`make submodules`).

## Cross-engine integrity (when scoring)

The benchmark's credibility rests on avoiding circularity: seeding from one engine's report and then scoring that same engine inflates its self-recall. Score every engine against the *union* blueprint and read cross-engine *differences*, not absolute self-recall. The three-species design is itself a control — a yeast peptide quantified at the E. coli ratio is a detectable error. Search settings must be unified across tools; document any untunable knob in `configs/<tool>/NOTES.md`. See `evident.yaml` `failure_modes` for the full list.
