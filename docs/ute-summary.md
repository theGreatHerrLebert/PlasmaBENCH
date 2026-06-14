# PYE1: simulated vs real — what we found

**TL;DR — TimSim reproduces your real PYE1 data well almost everywhere; the one real gap is
that the simulation over-separates spike-in ratios at *low* abundance, and we traced that to a
specific, fixable simulator behaviour.**

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

### 3. We found the cause — and it is fixable.
The smaller partner of each A/B pair is under-measured at low signal because the simulator
spreads each peptide's signal across many pixels and clips them at an intensity floor of **1**
in its own arbitrary units — far below the real timsTOF detector minimum (~11) — and builds the
peaks deterministically **without real ion-count statistics**, so a faint peptide's thinly-spread
signal is erased; the real-data noise we add comes too late in the pipeline to compensate.
The fix is to model the real detector (a realistic intensity floor + ion-count statistics)
instead of the deterministic clip. *(plasmabench-lowinput-overseparation — proposed, not yet built)*

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
