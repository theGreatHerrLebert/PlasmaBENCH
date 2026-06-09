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

`dia/sample{A,B}_ye/config.toml` are **scaffolding** rendered with placeholders. Two
things are still pending before they can run:

1. **The real PYE-plasma `.d`** (Ute is uploading it; goes under
   `data/raw/dia/plasma/`). The configs currently point at
   `…/plasma/PYE_PLASMA_PLACEHOLDER.d`.
2. **YE-only seeds** `seeds/seed_sample{A,B}_ye.csv` (human excluded). Not yet built.

Also confirm `gradient_length` matches the plasma run's RT span (Ute's plasma gradient
differs slightly from the blank — she says that's fine; set `--gradient-length` to match).

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
