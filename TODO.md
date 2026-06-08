# TODO

Open work for PlasmaBENCH. Mirror into the GitHub Project board so each item has an owner.

## Now (Stage 1 — simulated PYE1)

- [ ] **Get the reports from Ute** — Sample A and Sample B, both **DIA-NN** (`report.tsv`/`report.parquet`) and **FragPipe**. These seed the simulation. *Owner: distleu.*
- [ ] **Get a TimsTOF DIA blank `.d`** acquired with the study's DIA method, for noise injection (`add_real_data_noise`). A plasma-only `.d` would additionally unlock stage 2 (`superimpose_on_reference`). *Owner: distleu.*
- [ ] **Confirm the PYE1 ratios** — A = 90 % plasma / 2 % yeast / 8 % E. coli, B = 90 / 6 / 4. Expected B/A: human 1.0, yeast 3.0, E. coli 0.5.
- [ ] **Validate the seed builder** — `scripts/build_seed_from_report.py` handles DIA-NN `(UniMod:N)` and FragPipe mass-bracket notation, but the exact column names / mod syntax are version-dependent. Check against the real files once they arrive; adjust `plasmabench/peptides.py` if needed.
- [ ] **Build the seeds** — one `seed_sample{A,B}.csv` per sample; verify the HUMAN/YEAST/ECOLI split (`plasmabench/species.py`) matches the expected composition before locking.
- [ ] **Pin the rustims branch** — `feature/simulate-from-findings` carries the from_findings 6-column contract; not yet on main. Pin the `rustims/` submodule at the branch tip and rebuild the Docker image; the upstream `rustims:latest` does not have it.
- [ ] **Run TimSim DIA** for both samples; capture rustims commit + Docker digest per run.
- [ ] **Memory sizing** — like MHCbench, large from_findings seeds are RAM-hungry. Smoke-test on a head of the seed before the full run.

## Software setup
- [ ] Unify search params across DIA-NN / FragPipe (and later Spectronaut / MaxQuant); document untunable ones in `configs/<tool>/NOTES.md`.
- [ ] Confirm decoy export per tool.
- [ ] Capture tool license/version in `configs/<tool>/VERSION`.

## Stage 2 — SIM on real plasma (planned)
- [ ] Locate/acquire a real plasma-only `.d` (study method); wire `superimpose_on_reference`.
- [ ] Compare semi-synthetic vs. fully-synthetic vs. real experimental runs.

## Repo / infra
- [ ] Add collaborators; decide visibility.
- [ ] Decide raw-data hosting.
- [ ] CI: config lint + EVIDENT manifest validation.
