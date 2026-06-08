# PlasmaBENCH

Benchmark suite for comparing DIA software on TimsTOF **mixed-proteome plasma** data, with **TimSim**-generated synthetic ground-truth datasets seeded from real search-engine results.

The name nods to LFQbench / HYE — same three-species ratio idea, applied to a human-plasma background with yeast and *E. coli* spike-in.

## Goal

Reproduce, *in silico*, the mixed-proteome plasma benchmark of Ute Distler and colleagues (*Nat Commun* 2025, DOI [10.1038/s41467-025-64501-z](https://doi.org/10.1038/s41467-025-64501-z)) on the TimsTOF, so that DIA pipelines can be scored where the truth — which peptides are present, from which organism, at what ratio between samples — is exactly known. Because TimSim emits a blueprint alongside every simulated `.d`, we can quantify species-resolved FDR, recall, and quantification (ratio-recovery) accuracy.

This is the plasma/quant counterpart to [MHCbench](https://github.com/theGreatHerrLebert/MHCbench) (immunopeptidome recall/FDR). Same TimSim `from_findings` engine, same EVIDENT trust scaffolding, different sample type and oracle.

## The PYE1 design (Plasma · Yeast · E. coli)

Two samples with a fixed human-plasma background and complementary spike-in ratios. The plasma stays 1:1 while the two spike-in organisms move in opposite directions — the classic three-species ratio test.

| Species | Sample A | Sample B | Expected B/A | log2(B/A) |
|---|---|---|---|---|
| Human (plasma) | 90 % | 90 % | **1.00** | 0.00 |
| Yeast | 2 % | 6 % | **3.00** | +1.585 |
| *E. coli* | 8 % | 4 % | **0.50** | −1.00 |

A pipeline that quantifies well should place each species' peptide population tightly around its expected log2 ratio, with no cross-species leakage.

## How it works (TimSim `from_findings`)

Ute's real DIA-NN / FragPipe peptide reports for Sample A and Sample B are converted into TimSim's 6-column **findings** contract and simulated directly — so the simulated runs contain *exactly the peptides that were really detected*, at their real relative intensities. The A/B ratios are therefore **baked into the report intensities + organism membership**; we do **not** use TimSim's `[proteome_mix]` / dilution-factor path (that route is for fully synthetic FASTA digestion). See `simulations/stage1/README.md` for the contrast.

```
DIA-NN / FragPipe report (Sample A)  ─┐
DIA-NN / FragPipe report (Sample B)  ─┤→ scripts/build_seed_from_report.py
                                       │     (→ UniMod sequences, species tag,
                                       │      protein/sequence/intensity[/charge/rt/im])
                                       └──→ simulations/stage1/seeds/seed_sample{A,B}.csv
                                                  │
                                                  └──→ timsim <rendered config.toml> → SIM .d + blueprint
```

The findings contract (from rustims `feature/simulate-from-findings`): columns `protein, sequence, intensity` (required) + `charge, rt, im` (optional). **Optional columns are emitted as real NaN or omitted entirely — never empty strings** (TimSim's `load_findings` runs `dropna` over all present columns, so a single all-blank `rt`/`im` column would silently delete every row).

## Stages

### Stage 1 — fully simulated PYE1 (in flight)
Seed TimSim from the real Sample A / Sample B reports and simulate **DIA-PASEF** runs (optionally N replicates with light jitter). Background noise injected from a real TimsTOF blank acquired with the same DIA method. Oracle: TimSim blueprint + the known B/A species ratios above.

### Stage 2 — SIM spiked onto real plasma (planned)
Use TimSim's `superimpose_on_reference` to overlay the simulated spike-in signal onto a **real plasma-only `.d`**, giving a semi-synthetic sample with a realistic matrix but known spike-in truth. Compare against the real experimental TimsTOF runs from the *Nat Commun* study.

### Out of scope (for now)
Orbitrap / Thermo `.raw` — TimSim emits TimsTOF Bruker `.d` only. The study's Orbitrap arm is not reproducible here; the experimental Orbitrap data stand on their own.

## Software under test

DIA-NN and FragPipe to start (both already in Ute's report set), with room for Spectronaut / MaxQuant. Settings must be unified across pipelines (mass tolerance, peptide length, modifications, FDR target); parameters that cannot be set to the unified value are documented per tool in `configs/<tool>/NOTES.md`. Where possible, configure each tool to also output decoys for entrapment-style FDR checks (the three-species design already gives a natural cross-organism control: a yeast peptide quantified at the *E. coli* ratio is an error).

## Repository layout

```
PlasmaBENCH/
├── README.md            – this file
├── PROJECT.md           – kanban mirror
├── TODO.md              – open work, owners
├── evident.yaml         – EVIDENT claim manifest (validated in CI)
├── evident/             – submodule: theGreatHerrLebert/evident (read-only reference)
├── rustims/             – submodule: theGreatHerrLebert/rustims (pinned; from_findings build path)
├── plasmabench/         – Python helpers (report→UniMod, species assignment, metrics)
├── cases/               – per-stage EVIDENT narratives
├── configs/             – per-software DIA search configs (DIA-NN / FragPipe / …)
├── docker/              – compose for the timsim runtime
├── data/                – raw + intermediate (gitignored; fileshare links instead)
├── scripts/             – seed builder, config renderer, run orchestration
├── simulations/         – timsim run definitions and outputs
└── results/             – figures, tables, manuscript-bound outputs
```

Clone with submodules:

```bash
git clone --recurse-submodules git@github.com:theGreatHerrLebert/PlasmaBENCH.git
# or, after a plain clone:
git submodule update --init --recursive
```

## Running TimSim

Stage 1 needs the **`from_findings`** TimSim mode (6-column CSV input), which lives on the rustims `feature/simulate-from-findings` branch — not yet on `main`, so the upstream `ghcr.io/thegreatherrlebert/rustims:latest` image does **not** contain it. Two paths:

1. **Docker (preferred once published).** PlasmaBENCH builds its own runtime image from the pinned `rustims/` submodule and pushes to GHCR (`.github/workflows/docker-publish.yml`):

   ```bash
   make pull            # ghcr.io/thegreatherrlebert/plasmabench-timsim:from-findings
   make shell           # interactive GPU shell with timsim on PATH
   make timsim-help
   ```

2. **From source** until the image is published:

   ```bash
   make venv-source     # builds timsim from the pinned rustims submodule into .venv
   source .venv/bin/activate
   ```

   Details in [`docs/venv-setup.md`](docs/venv-setup.md).

### Data layout and per-machine paths

Configs reference inputs/outputs through container-fixed roots, so they stay portable:

| Container path | Host default | Override |
|---|---|---|
| `/work` | `<repo>` | n/a (always the repo) |
| `/data` | `<repo>/data` | `PLB_DATA` |
| `/sim`  | `<repo>/simulations` | `PLB_OUT` |

Raw inputs live under `data/raw/`:

```
data/raw/
├── dia/blanks/      – TimsTOF DIA-PASEF blank(s) for noise injection
├── dia/plasma/      – real plasma-only .d (stage 2 superimpose target)
└── fasta/           – plasma + yeast + E. coli reference FASTA(s)
```

To run on another machine, drop the same files at the same repo-relative paths under `data/raw/`, or set `PLB_DATA=/abs/path` in `docker/.env` (see [`docker/.env.example`](docker/.env.example)).

## Collaborators

| Slack | Name | Role | GitHub |
|---|---|---|---|
| dateschn | David Teschner | TimSim / simulation | theGreatHerrLebert |
| distleu | Ute Distler | study + reports (Sample A/B), real data | _TBD_ |
| tenzer | Stefan Tenzer | PI | _TBD_ |

ProteoBench is using the same study's dataset for its plasma module — a natural point of contact for cross-validation.

## Related work

- **Mixed-proteome plasma benchmark** — Ute Distler and colleagues, *Nat Commun* 2025, DOI [10.1038/s41467-025-64501-z](https://doi.org/10.1038/s41467-025-64501-z). The study PlasmaBENCH reproduces *in silico*; source of the Sample A/B definitions and the real DIA-NN / FragPipe reports.
- **TimSim preprint** — Teschner et al., *Complete Simulation of timsTOF PASEF Raw Datasets with Timsim*, [Research Square 2026](https://www.researchsquare.com/article/rs-9032301/v1). The simulator used for ground truth.
- **MHCbench** — sibling benchmark on JY/Raji immunopeptidome data; same TimSim `from_findings` engine and EVIDENT scaffolding.
- **LFQbench** — Navarro et al., *Nat Biotechnol* 2016, [10.1038/nbt.3685](https://doi.org/10.1038/nbt.3685) — the original three-species ratio benchmark this design follows.
