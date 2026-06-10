# Finding: DIA-NN 2.5 has worse two-sample ratio precision than 1.8 (quant)

## Question
Does DIA-NN 2.5 perform *worse* than 1.8 on quant? 2.5 is the more sensitive identifier
(it is why we use it for the one-shot blank-background ID step), but the user observed that
its quantitative ratios look noisier. This separates "quant" into two distinct axes and
tests each on both stages, then controls for the two confounds an independent reviewer
raised (precursor-set coverage, and human-anchoring).

## Two axes of quant (do not conflate)
- **Within-sample abundance recovery** — Spearman correlation of observed precursor intensity
  vs blueprint truth intensity, within a single run. Measures rank-abundance fidelity.
- **Cross-sample ratio precision** — IQR (interquartile range) of the per-precursor log2(A/B)
  fold change, per species. Measures how *tight* the fold-change estimates are. For a
  ratio-test benchmark (PYE1), this is the headline quant metric.

Ratio orientation: tables below are in the **A/B (paper) display convention**. Expected
values: human 0, yeast −log2 3 ≈ −1.585, E. coli +1. IQR is orientation-invariant (a sign
flip and constant shift do not change spread), so the precision comparison is unaffected by
convention.

## Results — axis 1: within-sample abundance correlation (Spearman; higher = better)
| stage | engine | yeast | ecoli |
|---|---|---|---|
| Stage 2 (fixed) | 1.8 | 0.938 | 0.944 |
| Stage 2 (fixed) | 2.5 | 0.915 | 0.919 |
| Stage 1 (blank) | 1.8 | 0.932 | 0.934 |
| Stage 1 (blank) | 2.5 | 0.923 | 0.931 |

→ 2.5 is **marginally but consistently lower** (~0.92 vs ~0.94, lower in all four cells).
Not a "tie" — a small, repeatable deficit — but the gap is minor relative to axis 2.

## Results — axis 2: A/B ratio precision (IQR of log2 fold change; lower = tighter = better)

Two computations are reported. **(a)** is the protein-level / linear-mean ratio table the
report figures use. **(b)** is the ion-level *control* that answers the coverage confound:
for each species it compares the two engines' full precursor sets **and** the subset of
precursors *both* engines quantified (`scripts/diag_quant_iqr.py`).

**(a) Protein-level, linear-mean (report figures):**
| stage | engine | yeast | ecoli |
|---|---|---|---|
| Stage 2 (fixed) | 1.8 | **0.21** | **0.11** |
| Stage 2 (fixed) | 2.5 | 1.04 | 0.91 |
| Stage 1 (blank) | 1.8 | **0.15** | **0.11** |
| Stage 1 (blank) | 2.5 | 0.28 | 0.74 |

**(b) Ion-level, full vs shared-precursor intersection** (n = precursors per species):
| stage | species | engine | n_full | IQR_full | n_shared | IQR_shared |
|---|---|---|---|---|---|---|
| Stage 1 | ECOLI | 1.8 | 7495 | 0.215 | 7235 | 0.209 |
| Stage 1 | ECOLI | 2.5 | 7797 | 0.448 | 7235 | **0.414** |
| Stage 1 | YEAST | 1.8 | 9682 | 0.226 | 9336 | 0.220 |
| Stage 1 | YEAST | 2.5 | 10154 | 0.429 | 9336 | **0.386** |
| Stage 2 | ECOLI | 1.8 | 7218 | 0.191 | 6922 | 0.181 |
| Stage 2 | ECOLI | 2.5 | 7498 | 1.642 | 6922 | **1.632** |
| Stage 2 | YEAST | 1.8 | 9933 | 0.351 | 9517 | 0.335 |
| Stage 2 | YEAST | 2.5 | 10509 | 1.614 | 9517 | **1.604** |

→ 2.5's ratio dispersion is **2× wider (Stage 1) to ~9× wider (Stage 2)** than 1.8's, in
both stages. n is thousands of precursors per cell, so this is not a small-sample effect.

## Confounds controlled (the two an independent reviewer flagged)

1. **Coverage / not-like-for-like — RULED OUT.** 2.5 reports more precursors (560–990 extra
   per species), and the worry was that those extra marginal IDs, not worse quant on shared
   analytes, drive the spread. On the **shared precursor intersection** (identical ions both
   engines quantified) the gap **persists almost entirely**: Stage-2 ecoli 2.5 goes 1.642 →
   1.632 when restricted to shared (vs 1.8's 0.181); yeast 1.614 → 1.604. The extra precursors
   move IQR by <0.02. The wider spread is 2.5 quantifying the *same* ions less precisely, not
   a coverage artifact.

2. **Human-anchor — DOES NOT AFFECT IQR.** Anchoring subtracts a single scalar (the human
   median log2 B/A) from every feature (`ratios.py:145`); replicates are already collapsed to
   one ratio per feature before this. A constant shift cannot change within-species IQR.
   Empirically confirmed: Stage-2 2.5's anchor is only −0.040 (tiny) yet its IQR is 1.6; the
   spread is intrinsic to the unshifted ratios, not anchor-induced. (Anchoring *would* matter
   for IQR only if ratios were pooled across run-pairs with different per-pair offsets — they
   are not; each table is one A-vs-B comparison.) Separately, Stage-1 2.5 has a large anchor
   (−1.67 ⇒ human B/A median ≈ +1.67), meaning 2.5's human ratio is itself badly off-center —
   a related accuracy issue, but orthogonal to the precision finding.

## Results — axis 2 accuracy companion: residual vs blueprint truth

IQR captures spread but is blind to bias. With exact blueprint truth we score the
per-precursor residual `e = observed_log2(A/B) − truth_log2(A/B)` (both human-anchored, so the
global intensity-scale convention is removed and what remains is species-specific quant error;
`scripts/diag_quant_residual.py`, ion level). `bias = median(e)` (accuracy), `IQR/MAD` (precision),
`med|e|` (combined):

| stage | species | engine | n | bias | IQR | MAD | med\|e\| |
|---|---|---|---|---|---|---|---|
| Stage 1 | ECOLI | 1.8 | 7460 | −0.032 | 0.214 | 0.087 | **0.092** |
| Stage 1 | ECOLI | 2.5 | 7750 | −0.076 | 0.445 | 0.195 | 0.215 |
| Stage 1 | YEAST | 1.8 | 9622 | +0.044 | 0.226 | 0.087 | **0.096** |
| Stage 1 | YEAST | 2.5 | 10081 | +0.065 | 0.425 | 0.204 | 0.211 |
| Stage 2 | ECOLI | 1.8 | 7155 | −0.009 | 0.187 | 0.094 | **0.095** |
| Stage 2 | ECOLI | 2.5 | 7396 | −0.338 | 1.635 | 0.564 | 0.613 |
| Stage 2 | YEAST | 1.8 | 9843 | +0.049 | 0.348 | 0.112 | **0.119** |
| Stage 2 | YEAST | 2.5 | 10351 | −0.183 | 1.606 | 0.664 | 0.655 |

Two things the residual shows that IQR alone did not:
- **1.8 is accurate *and* precise everywhere** — bias within ±0.05 log2, combined med|e| ≈
  0.09–0.12 log2 (≈ 7–9 % ratio error) across both stages and species.
- **2.5 in Stage 2 is also *biased*, not merely noisy** — fold changes compress toward 1:1
  (ecoli bias −0.34, yeast −0.18 log2), on top of the ~9× wider spread. On the blank background
  (Stage 1) 2.5's *accuracy* is fine (bias −0.08…+0.07) but its precision is already ~2× worse;
  it is the real-plasma background that degrades 2.5's accuracy as well. So 2.5's combined error
  (med|e|) is ~2× 1.8 in Stage 1 and ~6× in Stage 2.

## Interpretation
- The finding survives both confound controls and **recurs in both stages** (blank-background
  SIM and YE-on-real-plasma), with the from_findings median distortion already fixed. So it is
  not the Stage-2 pairing, not coverage, not anchoring, not the seed distortion. The most
  defensible statement: *on matched precursors, 2.5 estimates the A/B fold change with
  substantially more spread than 1.8, across both stage configurations.*
- **Mechanism — ISOLATED to DIA-NN 2.5's cross-run NORMALIZATION** (`diag_quant_method.py`,
  re-reads the same report under different quantity columns, no re-running DIA-NN):

  | engine | quant column | S2 ecoli IQR | S2 ecoli bias | S2 yeast IQR | S2 yeast bias |
  |---|---|---|---|---|---|
  | 2.5 | `Precursor.Normalised` (used) | 1.64 | −0.34 | 1.61 | −0.18 |
  | 2.5 | `Precursor.Quantity` (raw) | **0.23** | −0.02 | **0.30** | +0.04 |
  | 2.5 | `Ms1.Area` (raw) | **0.22** | +0.03 | **0.31** | −0.04 |
  | 1.8 | `Precursor.Normalised` (used) | 0.19 | −0.01 | 0.35 | +0.05 |
  | 1.8 | `Precursor.Quantity` (raw) | 0.16 | −0.02 | 0.21 | +0.03 |

  The cause is **not** MaxLFQ (a protein-level aggregation; our endpoint is precursor-level) and
  **not** the QuantUMS peak extraction (the *raw* `Precursor.Quantity` is fine). It is
  specifically 2.5's **cross-run normalization** layer: switching from the `.Normalised` to the
  raw column — same report, no re-search — collapses 2.5's Stage-2 IQR from 1.64 to 0.23 and
  removes the compression bias, landing right next to 1.8 (0.16/0.21). Both `.Normalised`
  columns (MS2 and MS1) are bad and both raw columns are fine, confirming the normalization
  step, not extraction. The asymmetry is the proof: 1.8's normalization is benign (0.19→0.16),
  2.5's is a ~7× regression (1.64→0.23). The normalization also *compresses* ratios toward 1:1
  (the Stage-2 bias), consistent with a normalization that assumes most signal is unchanged
  between runs and partially "corrects away" the genuine spike-in fold changes — worse on the
  complex real-plasma background (Stage 2) than the blank (Stage 1).

## Bottom line (the three-way split)
| axis | winner |
|---|---|
| Identification / sensitivity | **2.5** (more IDs — why we use it for blank ID) |
| Within-run abundance ranking | 1.8 (marginal, consistent) |
| Two-sample ratio precision | **1.8** (clearly tighter; survives shared-precursor + anchor controls) |

So "2.5 is worse on quant" is **true for ratio precision** — the metric that matters most for a
ratio-test benchmark — and the shared-precursor control makes it robust, not a coverage or
anchoring artifact. But it is **fixable, and the fix is selectable**: the regression lives
entirely in 2.5's cross-run normalization, so quantifying from raw `Precursor.Quantity` /
`Ms1.Area` (or re-running 2.5 with `--no-norm`) restores 1.8-level ratio precision.

## Actionable: choosing the quant method
- **Free, no re-search:** score 2.5 from `Precursor.Quantity` (raw) instead of `Precursor.Normalised`.
  Our loader already supports it (`load_observations(..., quant_col="Precursor.Quantity")`).
- **At search time:** re-run 2.5 with `--no-norm` (disables the cross-run/RT normalization);
  its `.Normalised` column then equals raw.
- **Cross-engine fairness:** we currently feed both engines `Precursor.Normalised`. That is
  benign for 1.8 but penalizes 2.5. The clean apples-to-apples comparison is **raw quantity for
  both** (isolates peak extraction from each engine's normalization choice); on raw, 2.5
  (0.23/0.30) ≈ 1.8 (0.16/0.21). Report the default-normalized numbers too, since they reflect
  what a user gets out-of-the-box.

## Still open / recommended before publishing
1. ~~Accuracy companion to precision.~~ **DONE** — residual-vs-truth table above
   (`diag_quant_residual.py`). Confirms 1.8 accurate+precise; reveals 2.5 also ratio-compressed
   in Stage 2.
2. **Uncertainty.** Add bootstrap CIs on the IQRs/biases, **resampling proteins/peptides** (not
   precursors — same-protein precursors are dependent; precursor bootstrap pseudoreplicates).
3. ~~Mechanism ablation.~~ **DONE** — isolated to 2.5's cross-run normalization
   (`diag_quant_method.py`); raw quantity recovers 1.8-level precision. Optional confirmation:
   re-run 2.5 with `--no-norm` and check `.Normalised` == raw.
4. **Distribution shape.** Show the residual/ratio distributions, not just IQR — to expose any
   multimodality or intensity dependence the scalar hides.
