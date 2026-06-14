# PlasmaBENCH claims audit (EVIDENT integration)

Every substantive claim in `results/report.html`, with how the data supports it. Support grades:

- **GT** — *ground-truth validated*: scored against the TimSim blueprint or an exact computation on the reports.
- **EXP** — *experimentally validated*: checked against real experimental data (Ute's G PYE1 runs).
- **CTRL** — *strong, control-backed*: not a single ground truth but confirmed by explicit controls (shared-set, matched-structure, adversarial/Codex review).
- **CAV** — *caveated / limitation*: true as stated but bounded; or an honestly-registered limitation.

Cross-engine integrity rule throughout: seed from one engine, score vs the **union** blueprint, read cross-engine *differences*, not absolute self-recall.

| # | Claim | Grade | Backing evidence / artifact | Caveat / failure mode |
|---|---|---|---|---|
| §2 | On **raw** precursor quantity, A/B ratio recovery agrees across DIA-NN 1.8/2.5/2.6 at ion & protein level (E. coli ≈ +1.0, yeast ≈ −1.6) | **GT** | `build_ratio_table` vs blueprint nominal; figs `stage1_diann{18,25}_ratios_{protein,ion}` | Protein uses raw Σ-precursor (not PG.MaxLFQ) — a deliberate choice (§8) |
| §2 | PG.MaxLFQ **compresses** the protein ratio (+0.91 vs ion +1.05) | **GT** | direct computation, both engines | — |
| §3 | Within-run quant correlation ρ ≈ 0.92–0.96, a tie across versions | **GT** | Spearman vs blueprint truth intensity; figs `stage1_diann{18,25}_quant` | rank-fidelity only; not a ratio statement |
| §4 | 1.8 true FDR ~2% precursor / **~7.5% protein**; 2.5/2.6 hold precursor <0.6%, protein ~1.4%; recall comparable | **GT** | `score_stage1_fdr.py` vs blueprint, blank-subtracted | protein FDR is group-credit (any/strict reported); ~79% of 1.8's is human |
| §4 | 2.6 ≈ 2.5 on FDR/recall (2.6 did not improve FDR) | **GT** | same scorer, 2.6 reports | — |
| §5/§9b | Cross-engine ID-overlap counts (SIM and real) | **GT** | exact set ops on q<0.01 IDs; Venn figs | protein "unique" on real not truth-labellable |
| §6 | Detection sensitivity & ratio bias rise/settle with **true** abundance; 2.5 not less sensitive at matched true abundance | **GT** | `build_true_abundance` (true-abundance bins) | reported-abundance view is the confound it corrects |
| §7 | Stage-2 (SIM spike-in on real plasma) recovers nominal ratios after the `from_findings` distortion fix | **GT** | sim-QC: blueprint A/B nominal (yeast ~3.0, ecoli ~0.5); rustims PR #407 | — |
| §7 | from_findings per-sample-median normalization distorted A/B (~0.7×) — **found & fixed** | **GT** | reproduced factor from seeds; fixed + re-sim QC; PR #407 merged | — |
| §7 | Stage-2 intensity calibration to real-PYE level | **CAV** | report-proxy calibration; E. coli 0.63× vs 0.53× target | proxy, not raw MS1/MS2 peak areas |
| §8 | 2.5/2.6's cross-run **normalization** inflates A/B ratio spread ~2–9× (`.Normalised`); **raw `Precursor.Quantity` recovers it** | **CTRL** | shared-precursor intersection (gap persists on identical ions), common-column set, residual-vs-truth; `diag_quant_{iqr,residual,method}.py`; Codex-reviewed | not MaxLFQ (protein rollup) nor QuantUMS extraction — both `.Normalised` bad, both raw fine |
| §8 | 2.6 inherits 2.5's normalization regression (not fixed) | **CTRL** + **EXP** | same on SIM and on real data | — |
| §8 | The raw-quant decision generalizes to real data | **EXP** | real reports show the same `.Normalised` spread; raw recovers | — |
| §9 | **SIM over-separates spike-in ratios at low abundance; matched real data is flat** | **CTRL** + **EXP** | matched 1-vs-1 (real = mean of 6 injection pairs); vs **true** abundance; human control flat; Codex-reviewed (corrected from an earlier 6-rep-averaging over-claim) | x-axis is reported-abundance **rank**, not matched absolute input; complete-case selection controlled only by same-pipeline comparison |
| §9 | Our unified 2.5/2.6 reproduce Ute's own report numbers (ecoli +0.89–0.90 / yeast −1.32–1.33) | **EXP** | `load_observations_real` on both reports | Ute's settings differ (LabL); cross-check, not identical pipeline |
| §9a | Cause = the **lower-intensity member of each pair is under-quantified at low abundance**; higher−lower gap = the over-separation | **GT** | per-member recovery log2(obs/truth) vs true abundance; gap matches measured over-sep to the decimal; member-symmetric; `plot_lowabund_underquant.py` | — |
| §9a | Mechanism = deterministic **sub-1.0 ion-count collapse** (truncation + `<1.0` per-pixel threshold + rounding, no count sampling); real-data noise injected *after* the threshold so it can't help | **CTRL** | Codex review of the rustims rendering path (load_findings.py:504/565, precursor.rs:165, dia.rs:504/673/731, assemble_frames.py:132/143) | the proposed fix is **unverified** — no after-fix SIM yet |
| §10 | Superposition limit: Stage 2 does not reproduce ion suppression / fill-time / detector saturation | **CAV** | by construction | scoped limitation |

## Headline supported conclusions
1. **TimSim is validated against real ground truth** — at matched abundance the ratio behaviour, ID counts, and version differences reproduce reality (our 2.5/2.6 = Ute's numbers). **[EXP]**
2. **DIA-NN 2.5/2.6's cross-run normalization inflates ratio spread; quantify from raw — holds on SIM *and* real.** **[CTRL+EXP]**
3. **One concrete simulator gap, mechanism pinned:** SIM over-separates spike-in ratios at low abundance because low-intensity synthetic peaks are erased by deterministic sub-1.0 ion-count thresholding; real data is flat. Fix proposed (Rust per-pixel count model), **not yet verified**. **[CTRL+EXP, fix pending]**

## Not-yet-supported / open (honest register)
- The §9a **fix** (per-pixel count model) is proposed but unvalidated — needs an after-fix SIM re-run scored against this same per-member-recovery target.
- Stage-2 calibration is a **report-proxy**, not raw peak areas.
- The low-abundance over-separation is established vs reported-abundance **rank**; an absolute-input-matched axis would strengthen it further.
