# TimSim gap: deterministic sub-1.0 ion-count collapse → low-abundance ratio over-separation

*Ready to file as a rustims issue. Validated against real experimental data (PlasmaBENCH PYE1);
the code-level mechanism below was confirmed by reading the rendering path (Codex review).*

## Summary
TimSim (`from_findings`) reproduces real PYE1 A/B spike-in ratios well at moderate-to-high
abundance, but **spuriously over-separates fold changes at low abundance** — simulated
low-abundance precursors get ratios *more extreme* than nominal, an effect absent in real data.
Root cause: the **lower-intensity member of each A/B pair is under-quantified at low abundance**
because the synthetic-signal renderer uses **deterministic fractional intensities with a hard
`< 1.0` per-pixel threshold + rounding** (no count sampling). The dimmer member loses
proportionally more of its sub-unit per-pixel contributions before they are ever written, so its
quantity drops and the fold change inflates.

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

## Code-level mechanism (confirmed against the rendering path)
1. **Integer event scaling, truncated, floored at 1**: `events = max(1, int(intensity / median *
   upscale_factor * multiplier))` (`load_findings.py:504/565`). Dim findings lose fractional
   events and eventually collapse to the same 1-event floor (matters below `upscale_factor`,
   default 100000).
2. **Deterministic multidimensional dilution**: events are split across frame × scan × charge ×
   isotope/fragment as `frame_abundance * scan_abundance * ion_abundance * total_events`
   (`precursor.rs:125/133`, `dia.rs:531/568`) — no Poisson/binomial count draw.
3. **Hard intensity floor of 1.0 + rounding + uint32**: per-pixel contributions below an
   intensity of `1.0` are **discarded** (`filter_ranged(…, 1.0, 1e9, …)` `precursor.rs:166`;
   `intensity_min.unwrap_or(1.0)` `dia.rs:504/772`), survivors rounded (`dia.rs:321/400`) and
   cast to uint32 (`tdf.py:238`). Note this is an **intensity** floor in the simulator's own
   units — *not* an ion count (there is no fractional ion), and it sits **far below the real
   timsTOF detector minimum (~11)**, so it is both the wrong unit and the wrong magnitude.

Per-peak intensity is proportional to `total_events`, so a fainter peptide has all per-pixel
contributions scaled down and more of them cross below the fixed `1.0` threshold; for DIA MS2 the
filter is applied **per fragment spectrum BEFORE frame aggregation** (`dia.rs:663`), so two sub-1.0
contributions at the same coordinate cannot even combine to survive. The under-recovery of the
dimmer member is therefore the **interaction of abundance-proportional fractional spreading with a
nonlinear, early per-contribution threshold** — *not* the floor value alone (spreading by itself
conserves intensity; 1.0 is permissive, below the real ~11). **This is a code-supported renderer-side
cause that is *capable* of producing the measured per-member under-quantification — it is NOT proven
to be the sole or quantitative driver; DIA-NN's own integration is not excluded.** NB: raising the
floor toward ~11 as a standalone change clips MORE and worsens it; and naive Poisson with sub-1
pixel means yields mostly zeros, not ~11 — so neither is a fix on its own.

**Why the existing noise mechanisms don't fix it (incl. real-data noise — observed):**
- `add_uniform_noise` is `abundance + abundance*noise` *renormalized to preserve total*
  (`utility.py:56`) — multiplicative chromatographic-shape jitter, NOT an additive floor; creates
  no signal where abundance is 0; defaults off.
- Reference-noise / real-data superimpose is genuinely additive, but it is a **separate background
  population added at its OWN m/z** — the noise step samples unrelated peaks from real frames and
  the superimpose adds whole real frames (`add_noise_from_real_data.py:110/128`). It therefore does
  not add intensity at the erased peptide's specific (m/z, RT, IM) coordinate, so it cannot
  reconstitute that peak's clipped quantity — regardless of ordering (additive signal is not
  order-dependent; the fact that it runs after the `<1.0` filter at `assemble_frames.py:132/143` is
  incidental). It only adds nearby interference that DIA-NN may partly integrate → a *less* extreme
  ratio (mild compression), which is why adding real-data noise **mildly reduced but did not remove**
  the over-separation (Stage-2 over-separates less than blank Stage-1).

## Proposed fix (candidates — to be settled by the audit below)
The bias lives in the **spreading × early-threshold interaction**, so the levers are:
1. **Relocate the threshold** — aggregate per-coordinate contributions FIRST, then apply a single
   intensity threshold (so sub-1.0 fragments that belong to the same peak can combine). This alone
   may remove much of the bias, with no count model.
2. **Discrete ion allocation + a gain model** — allocate integer ions to pixels and apply an
   explicit ion→recorded-intensity (detector gain) mapping with a properly-placed detector
   threshold. Plain `Poisson(expected_pixel_intensity)` is NOT sufficient (sub-1 means → mostly
   zeros), and adding it while keeping an ~11 floor could make loss worse.
Do **not** simply raise the 1.0 floor toward ~11 as a standalone change — it clips more and worsens
the dim member. **Decisive experiment:** audit the *rendered* signal (rendered-frame intensity vs
blueprint truth, per member, vs true abundance) to confirm the renderer — not DIA-NN integration —
is the locus, then ablate threshold-order / count-model and re-score against the per-member-recovery
target. Only then is the fix (and which lever) established.
Do NOT add events in `load_findings` (that changes blueprint truth and imposes compression) or a
constant to the EMG abundance (renormalized → adds no counts). A realistic additive low-count
floor is one candidate term within this count model, but the deterministic `< 1.0` truncation is
the demonstrated primary mechanism. The PlasmaBENCH per-member recovery + ratio-vs-abundance
curves are the validation target: after the fix they should flatten to match real data.

## Reproduce
PlasmaBENCH: `scripts/plot_real_vs_sim_compression.py` (real-vs-SIM ratio vs abundance),
`scripts/plot_lowabund_underquant.py` (per-member recovery decomposition). Real-data loader:
`plasmabench.plots.loader.load_observations_real`.
