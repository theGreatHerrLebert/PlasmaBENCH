# Stage 2 — YE spike-in simulated *over real plasma* (superimpose)

Stage 1 simulates all three species and injects noise **sampled from a blank**
(`add_real_data_noise=true`). Stage 2 instead **superimposes** the simulated
yeast + E. coli signal onto a **real PYE-plasma `.d`** (`superimpose_on_reference=true`,
the two modes are mutually exclusive). The human plasma background is then *real*, not
simulated — so:

- the seed is **YE-only** (human peptides removed; human comes from the reference run);
- human carries **no blueprint truth** → score with `--human-scope background_unknown`
  (P5/P3 restrict to the simulated spike-ins; human only feeds the empirical P1/P2/P4);
- this is the condition that tests whether the reversed low-abundance quant trend is a
  blank-noise artifact (see memory `stage1-noise-blank-not-plasma`).

## Status

**Reference plasma is DOWNLOADED** (2026-06-09): the plasma used for the PYE mix is
`data/raw/dia/plasma/PYE_plasma_2022_135/` — runs `G250506_080_Slot2-16_1_17973.d` and
`…_081_…_17974.d` (2 of 6 reps fetched; rest available from the Seafile share). RT span
**0.9–2460 s (41 min)**, longer than the blank — so the seed IDs (elute to ~2121 s) **fit
natively**: set `--gradient-length 2460` and YE seeds need **no aggressive RT compression**
(unlike Stage-1's `--rt-max 1480` for the blank).

Still pending: **YE-only seeds** `seeds/seed_sample{A,B}_ye.csv` — build with
`build_seed_from_report.py --keep-species YEAST ECOLI` from the same G-site PYE1 report used
for Stage-1 seeds (drop `--rt-max`, or use `--rt-max 2400` to stay inside the 2460 s span).

## Re-render once the plasma .d + YE seeds exist

```bash
python scripts/render_configs.py \
  --base simulations/stage1/configs/base-dia.toml \
  --seeds-dir simulations/stage2/seeds --out-dir simulations/stage2/dia \
  --samples A B --seed-suffix _ye --path-mode host \
  --experiment-prefix PLB-S2-DIA \
  --superimpose \
  --reference-path /data/raw/dia/plasma/<PYE_plasma>.d \
  --gradient-length <plasma_span_seconds>
```

Then `timsim <config.toml>` per sample, search with `configs/dia-nn/run.sh`, and run the
same panels with `--human-scope background_unknown`.
