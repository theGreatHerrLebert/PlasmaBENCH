# PYE1: simulated vs real — what we found

**TL;DR — TimSim reproduces your real PYE1 data well almost everywhere; the one real gap is
that the simulation over-separates spike-in ratios at *low* abundance, which we have narrowed to
a specific part of the simulator (with candidate fixes and a test to confirm it).**

Every finding below is backed by a registered, one-command-replayable claim in our EVIDENT
manifest (`evident.yaml`); the claim id is in parentheses.

---

### 1. The simulation matches reality at normal abundance.
At moderate-to-high abundance the simulated A/B ratios, the ID counts, and the differences
between DIA-NN versions all line up with your real G PYE1 runs — in fact our re-search with
2.5/2.6 reproduces your own reported ratios (E. coli +0.89–0.90, yeast −1.32–1.33 in log₂).
*(plasmabench-real-vs-sim-validation)*

### 2. The one real deviation: low-abundance over-separation.
For the *faintest* peptides the simulation reports fold-changes that are **too extreme**
(spike-in ratios blown outward), whereas your real data stays flat there — so the compression
behaviour you noticed really is different between simulated and real, and it lives at low input.
*(plasmabench-real-vs-sim-validation, report §9)*

### 3. We have narrowed where the gap comes from — with a fix to test.
The faint partner of each A/B pair is the one that's under-measured at low signal. The strongest
code-level *candidate* cause (not yet proven the only one) is in how the simulator renders faint
peaks: it spreads each peptide's signal across many pixels and applies an intensity threshold to
each fragment *before* combining them, so for a faint peptide a large fraction of that thinly-spread
signal falls under the threshold and is dropped — and the real-data noise we add can't recover it,
because it is separate background at its own m/z (it lands elsewhere, not on the faint peptide's
peak). Your real data simply does not show this — it stays flat. The fix is to change how faint
signal is rendered (combine a peak's pieces *before* thresholding, and/or model discrete ions with
detector gain) — **not** raising the intensity floor, which would drop more; we'd confirm it with a
rendered-signal audit. *(plasmabench-lowinput-overseparation — candidate cause + fix, not yet built)*

### 4. DIA-NN 2.6 behaves essentially like 2.5 on quantification.
2.6 did **not** change the quant behaviour relative to 2.5 — both are more sensitive than 1.8,
but their default normalisation makes the A/B ratios noticeably noisier; reading the **raw**
intensity column instead brings 2.5/2.6 back in line with 1.8.
*(plasmabench-diann-crossrun-normalization)*

### 5. ID overlap, 1.8 vs 2.x (the thing you asked about).
On the real runs, 2.5/2.6 identify ~23–25% more precursors than 1.8 and recover ~92% of 1.8's
IDs plus a large extra set — but 1.8 still keeps a small unique fraction (not a strict subset),
and 2.5 vs 2.6 are nearly identical.
*(plasmabench-ideoverlap-real)*

### 6. On the simulated data (where we have ground truth), 2.5/2.6 are far better FDR-calibrated.
At the same 1% setting, 1.8's true protein error rate is ~7.5% vs ~1.4% for 2.5/2.6, at
comparable recall — but this can only be measured on the simulation, since real data has no
ground truth.
*(plasmabench-stage1-simulated-truth, plasmabench-stage2-superimpose-truth)*

---

**Reproducibility.** Each finding has a one-command replay (`make <claim>-eval`) that re-derives
its figure/metric from your raw data, runnable in a Docker image with the data mounted — so any
of these can be re-checked independently.
