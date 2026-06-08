# Stage 1 — fully simulated PYE1 (Sample A + Sample B)

Seed TimSim from Ute's real DIA-NN / FragPipe reports and simulate DIA-PASEF
runs for both PYE1 samples. The simulated `.d` files contain exactly the
peptides that were detected, at their real relative intensities; the A/B
spike-in ratios are therefore intrinsic to the seed.

## The two routes (why we use `from_findings`, not `proteome_mix`)

| | `from_findings` (THIS stage) | `proteome_mix` + dilution CSV |
|---|---|---|
| Source of peptides | real Sample A/B reports | in-silico digest of FASTAs |
| Source of A/B ratios | baked into report intensities + organism membership | `dilution_factors.csv` per FASTA |
| Reflects "what was detected" | yes | no (idealized) |
| Config | `from_findings = true`, `[proteome_mix] proteome_mix = false` | `from_findings = false`, `[proteome_mix] proteome_mix = true` |

We take the first row. The expected B/A ratios (human 1.0, yeast 3.0, E. coli
0.5) emerge from the data, and the benchmark question becomes: do *other*
search engines, run on the SIM, recover them?

## Inputs (TODO — from Ute)

- **Sample A / Sample B reports** — DIA-NN `report.tsv`/`report.parquet` and/or
  FragPipe `psm.tsv`/`ion.tsv`. Converted to seeds by
  `scripts/build_seed_from_report.py` (one run per sample).
- **TimsTOF DIA blank `.d`** — `data/raw/dia/blanks/BLANK_PLACEHOLDER.d`, for
  `add_real_data_noise`. Wire the real filename into `configs/base-dia.toml`.
- **FASTA** — `data/raw/fasta/plasma_yeast_ecoli.fasta` (human plasma + YEAST +
  ECOLI). Ignored by from_findings but referenced for downstream tooling and
  for the search-engine configs.

## Build the seeds

```bash
python scripts/build_seed_from_report.py --engine diann \
  --report /path/to/sampleA/report.tsv --sample A \
  --out simulations/stage1/seeds/seed_sampleA.csv
python scripts/build_seed_from_report.py --engine diann \
  --report /path/to/sampleB/report.tsv --sample B \
  --out simulations/stage1/seeds/seed_sampleB.csv
```

The builder prints a **species audit** (HUMAN/YEAST/ECOLI row counts). Check it
against the expected composition before locking the seeds — a wrong organism
map silently corrupts the ratio oracle.

## Run

```bash
# Build timsim from the pinned submodule (until the Docker image carries from_findings):
make venv-source && source .venv/bin/activate

# One-shot: build seeds (if REPORT_* given) → render configs → run both samples.
ENGINE=diann \
  REPORT_A=/path/to/A/report.tsv REPORT_B=/path/to/B/report.tsv \
  bash scripts/run_stage1.sh

# Plan only:
DRY_RUN=1 bash scripts/run_stage1.sh
```

Equivalent step-by-step:

```bash
python scripts/render_configs.py --base simulations/stage1/configs/base-dia.toml \
  --seeds-dir simulations/stage1/seeds --out-dir simulations/stage1/dia \
  --samples A B --path-mode host
for s in A B; do timsim simulations/stage1/dia/sample$s/config.toml; done
```

## TimSim version

Needs the **`from_findings`** mode (6-column CSV: `protein, sequence, intensity`
required + `charge, rt, im` optional), which is on rustims
`feature/simulate-from-findings`, not yet on main — so the upstream
`rustims:latest` Docker image lacks it. The `rustims/` submodule is pinned to
that branch; build from source via `make venv-source` until the image is
rebuilt.

## Memory note

Large from_findings seeds are RAM-hungry (MHCbench's 54k-peptide seed peaked
>20 GB and needed a big-memory host). Smoke-test on a `head` of the seed before
committing to a full run.

## TODO before the runs
- [ ] Reports + blank `.d` received and wired into `configs/base-dia.toml`.
- [ ] `build_seed_from_report.py` validated against the actual report columns.
- [ ] `gradient_length` and DIA window settings matched to the study method.
- [ ] Species audit reviewed for both samples.
