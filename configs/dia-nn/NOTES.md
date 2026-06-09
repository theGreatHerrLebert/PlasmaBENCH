# DIA-NN — config notes

Software-under-test runs **two DIA-NN versions** (1.8.1 and 2.5.0) on the same
SIM `.d` files with the **identical** unified settings file
[`diann-pye1.cfg`](diann-pye1.cfg). The only intended differences are each
version's built-in defaults — that contrast *is* the experiment.

- Settings: `diann-pye1.cfg` (validated to parse cleanly on **both** binaries —
  no unknown options on 1.8.1 or 2.5.0).
- Runner: `bash run.sh <1.8|2.5> <fasta> <out_dir> <runA.d> <runB.d>`
- Versions: see [`VERSION`](VERSION). Binaries are gitignored under `tools/`;
  fetch with `bash scripts/get_diann.sh`.

## Unified settings (diann-pye1.cfg)

| Setting | Value | Flag | Rationale |
|---|---|---|---|
| Mode | library-free, predicted lib | `--fasta-search --gen-spec-lib --predictor` | matches Ute's library-free DIA-NN run (report has `Predicted.*` columns) |
| Enzyme | Trypsin (no cut before P) | `--cut K*,R*,!*P` `--met-excision` | real trypsin does not cut K\|P / R\|P; `K*,R*` (cut before P) split ~3.8% of seed peptides out of the library. `!*P` + N-term Met excision match the seed's real-trypsin peptides and the known-good HeLa 2.5 config |
| Missed cleavages | 1 | `--missed-cleavages 1` | covers 100% of seed peptides (max 1 MC under trypsin/P) — **TODO(log): confirm vs Ute** |
| Peptide length | 7–30 | `--min/max-pep-len` | matches seed builder bounds |
| Precursor charge | 1–4 | `--min/max-pr-charge` | seeds contain charges 1–4 |
| Precursor m/z | 300–1800 | `--min/max-pr-mz` | DIA windows cover 300–1165; range left wide |
| Fragment m/z | 200–1800 | `--min/max-fr-mz` | DIA-NN default |
| Fixed mod | Carbamidomethyl (C) | `--unimod4` | plasma is alkylated; the only mod in Ute's report (`UniMod:4`) |
| Variable mods | **none** | `--var-mods 0` | the SIM ground truth contains **no** variable mods (seed has only fixed CAM), so searching Ox(M)/Ac would only add search space. Keeps recall/FDR clean and identical across tools. **TODO(log): confirm Ute used none** |
| FDR | 1% precursor + protein | `--qvalue 0.01 --matrix-qvalue 0.01` | standard |
| MBR / cross-run | on | `--reanalyse` | A+B analysed together → quant matrix; needs ≥2 `.d` (pass A and B) |
| Mass accuracy | **fixed 15/15 ppm** | `--mass-acc 15 --mass-acc-ms1 15` | **MUST be fixed, not auto.** `--mass-acc 0` (auto) → **0 IDs** on the SIM (see "Zero-IDs gotcha" below); 15 ppm → ~27k precursors. DIA-NN's own timsTOF advice is 10–15 ppm |
| Scan window | auto | _(omit `--window`)_ | inferred from data |

## Zero-IDs gotcha (2026-06-09) — auto mass accuracy kills the SIM search

The first stage-1 run returned **0 IDs at 1% FDR** on both samples despite the
SIM being clean (1e10 MS1/MS2 signal, clean isotope envelopes, realistic b/y
fragment ladders, Gaussian elution, 100% precursor coverage, mass calibrated to
3.6 ppm). Root cause was **`--mass-acc 0 --mass-acc-ms1 0` (auto)**: on the
ultra-clean simulated m/z, DIA-NN's auto-calibration collapsed to a window where
targets could not separate from decoys.

Proof (same prebuilt `--lib`, same `--reanalyse`, same A+B):
- `--mass-acc 0  --mass-acc-ms1 0`  → **0** IDs
- `--mass-acc 15 --mass-acc-ms1 15` → **27,152 / 27,444** precursors, ~5,150 protein
  groups/sample, RT-prediction accuracy 0.011–0.017 (was 0.229).

NOT the cause (each independently ruled out): intensity / dynamic range (the SIM
signal is fine), `--reanalyse`/MBR (works once mass-acc is fixed), the
"separate pipeline step" advisory (benign — the known-good HeLa 2.5 SIM run prints
it too and yields IDs), and library coverage (our digest reproduces 96.2% of seed
peptides). Reference known-good run:
`/scratch/timsim-demo/TIMSIM-HeLa-RUSTW-SMOKE/diann/` (DIA-NN 2.5.0, one-step).

`--smart-profiling` was dropped to mirror that known-good config; the isolation
test searched the prebuilt `--lib`, so re-validate the one-step cfg end-to-end on
the next full run.

## 1.8.1 vs 2.5.0 — inherent differences (NOT flag-controlled; this is the point)

- **Quantification:** 2.x quantifies with **QuantUMS** by default; 1.8.1 uses the
  legacy MaxLFQ-style path. Expect this to dominate the quant (ratio-accuracy)
  differences between the two on identical input.
- **DL models:** 2.5.0 ships newer retention-time / spectrum / charge predictors
  than 1.8.1. Affects the in-silico library and thus IDs.
- **Output format:** 2.x defaults to `.parquet`; we pass `--out report.tsv` so
  `build_seed_from_report.py --engine diann` (and the metrics code) read both the
  same way. `report.parquet` is also readable by the seed builder if preferred.
- `--predictor` is explicit-on in 1.8.1 and effectively default in 2.x; passed in
  both for an identical command line (2.5.0 accepts it without warning).
- **No `--gen-spec-lib` / no `--window`** in the unified config, matching the proven
  invocation in rustims `validate/diann_executor.py` (used by MHCbench): just
  `--fasta-search --predictor`. Two 2.5.0 notes:
  - `--window 0` warns "scan window radius should be a positive integer" — removing
    `--window` (let 2.5 auto-infer) clears it. 1.8.1 read `0` as auto.
  - 2.5.0 prints "the in silico-predicted library must be generated in a separate
    pipeline step ... now without activating FASTA digest" for ANY one-step
    library-free run (it persists with or without `--gen-spec-lib`). It is an
    **advisory, non-fatal** — FASTA digest still runs (it generates ~5.8M precursors
    / a 30,767-protein library and searches). MHCbench's one-step runs hit the same
    advisory and produced valid results. If recall ever looks degraded, switch 2.5 to
    the explicit two-step (generate `.speclib`, then search with `--lib`).
  - `--gen-spec-lib`/`--out-lib` dropped because we don't need the library written
    out (MHCbench doesn't either); valid on 1.8.1 too, so the config stays unified.

## Open items

- **TODO(Ute's DIA-NN log):** not yet in the repo (still on Seafile `DIA_NN/`).
  Once available, reconcile missed-cleavages, variable mods, mass-accuracy, and
  the exact 1.8.x patch version against the table above.
- **FASTA: built** — `data/raw/fasta/plasma_yeast_ecoli.fasta` (30,818 seqs:
  human reviewed 20,431 / EC1118 yeast 5,984 `_YEAS8` / E. coli K12 4,403),
  reproducible via `scripts/get_fasta.sh`. QC'd at **99.1% seed-peptide coverage**
  (`scripts/qc_fasta_coverage.py`). Ute's exact FASTA would close the last <1% but
  is no longer a blocker.
- Both runs are blocked on a SIM `.d`, which needs a working `timsim`
  (from_findings build).

## Seed relevance

`build_seed_from_report.py --engine diann` reads `report.tsv` / `report.parquet`.
Columns consumed: `Protein.Names` (species via `_HUMAN`/`_YEAS8`/`_ECOLI`) with
`Protein.Group`/`Protein.Ids` fallback, `Modified.Sequence` (`(UniMod:n)`),
`Precursor.Charge`, `Precursor.Quantity`, `RT` (min→s), `IM` (1/K0), `Run`.
Validated against Ute's actual G/L exports.
