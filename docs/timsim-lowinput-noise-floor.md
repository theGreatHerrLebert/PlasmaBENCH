# TimSim gap: low-input synthetic peaks lack an additive noise floor → ratio over-separation

*Ready to file as a rustims issue. Validated against real experimental data (PlasmaBENCH PYE1).*

## Summary
TimSim (`from_findings`) reproduces real PYE1 A/B spike-in ratios well at moderate-to-high
abundance, but **spuriously over-separates fold changes at low abundance** — simulated
low-abundance precursors get ratios *more extreme* than nominal, an effect absent in real data.
Root cause: the **lower-intensity member of each A/B pair is under-quantified at low abundance**
because synthetic peaks have no realistic additive noise floor; the small peak collapses below
where DIA-NN quantifies it faithfully, inflating the fold change.

## Evidence (PlasmaBENCH, DIA-NN 1.8, raw Precursor.Quantity)
Validated against Ute's real G250506 PYE1 runs (same batch as the Stage-2 reference plasma),
searched with the identical config (only the binary differs across DIA-NN 1.8/2.5/2.6).

1. **Real (matched 1-vs-1) is flat across abundance** — E. coli ≈ +0.95, yeast ≈ −1.40 (A/B,
   log2) at every abundance octile. SIM **over-separates at low abundance** — E. coli +1.4→+1.8,
   yeast −2.3→−2.4 at the lowest octile, converging to nominal only at high abundance. Engine-
   independent (1.8/2.5/2.6), both blank (Stage 1) and real-plasma (Stage 2) backgrounds.

2. **It is not the seed.** The blueprint truth ratio is exactly nominal and *perfectly flat*
   across abundance (imposed mode: the ratio is baked in independent of abundance). Binning by
   *true* (blueprint) abundance — not reported — the observed ratio still over-separates, so it
   is not a reported-abundance / complete-case-selection artifact. It is not MBR (a `--reanalyse`
   no-MBR control leaves it intact).

3. **Per-member decomposition pins it.** Quant recovery `log2(observed/truth)` vs true-abundance
   octile, split by the higher/lower member of each A/B pair:

   | low → high octile | higher member | lower member | over-sep (hi−lo) |
   |---|---|---|---|
   | E. coli | +0.04 → −0.07 (flat) | **−0.36** → −0.06 | +0.40 → −0.01 |
   | yeast   | +0.11 → −0.07 (flat) | **−0.51** → −0.06 | +0.62 → −0.02 |

   The higher member is recovered accurately at all abundances; the lower member is under-
   quantified by ~0.4–0.6 log2 (≈ 1.3–1.5×) in the bottom octile. The higher−lower gap **equals**
   the measured over-separation to the decimal — it fully accounts for the effect. Member-
   symmetric (holds where A is the big side *and* where B is) → it is about small-vs-big peak,
   not A/B or species.

## Why real data doesn't show it
Real low-abundance peaks ride on an **additive noise floor / co-eluting interference** that
*lifts the small peak*, keeping the A/B ratio accurate. A clue: Stage 2 (simulated spike-in
superimposed on *real* plasma) over-separates **less** than Stage 1 (blank) — the real
background partially supplies the missing floor. TimSim's purely synthetic peaks have no such
additive baseline, so the small peak falls below the faithful-quant regime.

## Proposed fix
Add a **realistic additive intensity floor / low-count noise baseline** to simulated precursor
peaks (and/or intensity-dependent quant jitter that mimics the detector at low ion counts), so
small peaks retain quantifiable signal rather than collapsing below DIA-NN's quant limit. The
PlasmaBENCH analysis above is the validation target: after the fix, the SIM per-member recovery
curves (and the ratio-vs-abundance curve) should flatten to match real data.

## Reproduce
PlasmaBENCH: `scripts/plot_real_vs_sim_compression.py` (real-vs-SIM ratio vs abundance),
`scripts/plot_lowabund_underquant.py` (per-member recovery decomposition). Real-data loader:
`plasmabench.plots.loader.load_observations_real`.
