# NEXT_STEPS

Picking up stage 1 after the 2026-06-09 session. Status is **unblocked**: the
SIM is good and DIA-NN 2.5 now produces real IDs. See
`configs/dia-nn/NOTES.md` ("Zero-IDs gotcha") and memory `diann-zero-ids-mass-acc`
for the full debugging trail.

## Where things stand

- **SIM `.d` files are validated** (`simulations/stage1/dia/sample{A,B}_imposed/…/*.d`).
  Don't re-simulate — frames, isotopes, fragments, elution, coverage all check out.
  The 6-order intensity spread is fine; the intensity model needs **no** change.
- **DIA-NN 2.5 config is fixed and validated end-to-end** →
  `results/diann-2.5-AB-v3/report.parquet`: A 27,104 / B 27,133 precursors,
  ~5,130 protein groups/sample. Root cause of the earlier 0 IDs was
  `--mass-acc 0` (auto); now fixed to `--mass-acc 15 --mass-acc-ms1 15`.
- **First species-resolved sanity check passed** (imposed seed, true ratios
  human 1.0 / yeast 3.0 / ecoli 0.5): recovered E. coli 0.47 (≈exact), yeast 1.85
  (right direction, compressed), human 0.76 (anchor skewed). Direction correct for
  all three; the compression/anchor offset is expected DIA-NN normalization
  behaviour, to be handled by the scorer.

### Discard these (stale, 0-ID, pre-fix)
- `results/diann-2.5-AB/`, `results/diann-2.5-AB-v2/` — produced with `--mass-acc 0`.
  (`results/diag-diann-2.5-*` are throwaway diagnostics; keep or delete as you like.)

## Next steps (in order)

1. **DIA-NN 1.8.1 run — completes the version axis.**
   ```bash
   THREADS=14 bash configs/dia-nn/run.sh 1.8 data/raw/fasta/plasma_yeast_ecoli.fasta \
     results/diann-1.8-AB \
     simulations/stage1/dia/sampleA_imposed/PLB-S1-DIA-sampleA-imposed/PLB-S1-DIA-sampleA-imposed.d \
     simulations/stage1/dia/sampleB_imposed/PLB-S1-DIA-sampleB-imposed/PLB-S1-DIA-sampleB-imposed.d
   ```
   1.8.1 honours `--out report.tsv` (2.5 writes `report.parquet`; `run.sh` now
   accepts either). Confirm the fixed config parses cleanly on 1.8.1 too — it
   shares `diann-pye1.cfg`. Watch the auto-vs-fixed mass-acc behaviour: 1.8 may
   tolerate auto, but keep 15/15 for a unified command line.

2. **Build the stage-1 scorer** (`make stage1-eval` is still a stub that exits 64).
   This is the main remaining piece. Per memory `stage1-eval-spec`:
   - truth = TimSim blueprint ∩ across engines; exclude singletons.
   - species-resolved recall + true FDR (cross-organism leakage = error).
   - quant: **human-anchor normalization** before B/A ratios — this is what
     corrects the human-0.76 / yeast-1.85 skew seen above (global norm is fooled
     by yeast+ecoli changers outnumbering the human anchor ~5:1).
   - report each engine against the **union** blueprint (avoid seed-engine
     circularity); read cross-engine differences, not absolute self-recall.
   - Writes `results/stage1/metrics.json`; wire it into `make stage1-eval` and the
     EVIDENT manifest (`evident.yaml`).

3. **FragPipe** as the second engine (already in Ute's report set). Add
   `configs/fragpipe/` (config + `VERSION` + `NOTES.md`), unify search params with
   DIA-NN, run on the same SIM `.d`, feed into the scorer.

4. **Cross-check the `observed`-mode seeds** (`seed_sample{A,B}.csv`, no `_imposed`).
   The validated run used the imposed (clean-ratio) seeds; the observed seeds carry
   the real DIA-NN-compressed intensities and are the realism arm. Same DIA-NN run,
   compare ratio recovery against the imposed run.

## Gotchas to remember

- **Never `--mass-acc 0` on the SIM** — auto-calibration collapses to 0 IDs on the
  ultra-clean simulated m/z. Fixed 10–15 ppm always.
- The DIA-NN 2.5 *"in silico-predicted library must be generated in a separate
  pipeline step"* warning is **benign** — the known-good HeLa SIM run prints it too.
  Do not "fix" it by going two-step.
- DIA-NN 2.5 main report is always `report.parquet` regardless of `--out` extension.
- Known-good reference run for comparison:
  `/scratch/timsim-demo/TIMSIM-HeLa-RUSTW-SMOKE/diann/` (DIA-NN 2.5.0).
