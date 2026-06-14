"""Assemble the PlasmaBENCH results into one self-contained HTML report.

Embeds the figures from results/figures/ as base64 (portable single file) and
weaves in the narrative + key-findings tables. Regenerate after new results:

    python scripts/build_report.py --out results/report.html
"""
from __future__ import annotations

import argparse
import base64
import html
from pathlib import Path

FIGDIR = Path("results/figures")


def img(name: str, caption: str) -> str:
    p = FIGDIR / name
    if not p.exists():
        return f'<p class="missing">[missing figure: {html.escape(name)}]</p>'
    b64 = base64.b64encode(p.read_bytes()).decode()
    return (f'<figure><img src="data:image/png;base64,{b64}" alt="{html.escape(name)}"/>'
            f'<figcaption>{html.escape(caption)}</figcaption></figure>')


def table(headers, rows) -> str:
    h = "".join(f"<th>{html.escape(str(c))}</th>" for c in headers)
    r = "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>"
                for row in rows)
    return f'<table><thead><tr>{h}</tr></thead><tbody>{r}</tbody></table>'


CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:1100px;
margin:0 auto;padding:2rem;color:#1a1a1a;line-height:1.55}
h1{border-bottom:3px solid #4E79A7;padding-bottom:.3rem}
h2{margin-top:2.4rem;border-bottom:1px solid #ddd;padding-bottom:.2rem;color:#2a4d69}
h3{margin-top:1.6rem;color:#444}
figure{margin:1.2rem 0;text-align:center}
img{max-width:100%;border:1px solid #e0e0e0;border-radius:6px}
figcaption{font-size:.86rem;color:#666;margin-top:.4rem;text-align:left}
table{border-collapse:collapse;margin:1rem 0;font-size:.9rem}
th,td{border:1px solid #ccc;padding:.35rem .7rem;text-align:left}
th{background:#f0f4f8}
.key{background:#f7faf7;border-left:4px solid #59A14F;padding:.6rem 1rem;margin:1rem 0}
.caveat{background:#fdf6f0;border-left:4px solid #E15759;padding:.6rem 1rem;margin:1rem 0}
.missing{color:#b00;font-style:italic}
code{background:#f3f3f3;padding:.1rem .3rem;border-radius:3px;font-size:.88em}
small{color:#888}
"""


def build() -> str:
    S = []
    S.append('<!doctype html><html lang="en"><head><meta charset="utf-8">'
             '<meta name="viewport" content="width=device-width, initial-scale=1">'
             f'<title>PlasmaBENCH Report</title><style>{CSS}</style></head><body>')
    S.append("<h1>PlasmaBENCH — TimSim-grounded benchmark of DIA software on mixed-proteome plasma</h1>")
    S.append("<p><small>Generated from <code>results/figures/</code>. Self-contained (figures embedded).</small></p>")

    # ---- Overview ----
    S.append("<h2>1. Design</h2>")
    S.append("""<p>PYE1 is a three-species ratio test (à la LFQbench/HYE): human <b>plasma</b>
    background stays 1:1 between Sample A and B, while <b>yeast</b> moves 2%→6% (expected
    B/A = 3.0) and <b>E. coli</b> 8%→4% (B/A = 0.5) in opposite directions. <b>TimSim</b>
    generates synthetic TimsTOF DIA <code>.d</code> from real DIA-NN findings
    (<code>from_findings</code>), emitting an exact ground-truth <i>blueprint</i> alongside
    each run — so recall, FDR, and ratio recovery are scored against truth, not consensus.
    Software-under-test: <b>DIA-NN 1.8.1 vs 2.5</b> (the version Ute used vs current).</p>""")
    S.append("""<p>Two simulation regimes: <b>Stage 1</b> simulates all three species and injects
    noise sampled from a real blank run; <b>Stage 2</b> superimposes only the yeast+E.coli
    spike-in onto a <i>real</i> plasma <code>.d</code> (human = real background, no blueprint
    truth). Every analysis panel is built on a two-table model (truth from the blueprint,
    observations from the engine report) and was twice code-reviewed (incl. an independent
    Codex review of the design and the code).</p>""")

    # ---- Stage 1: ratios ----
    S.append("<h2>2. Stage 1 — ratio recovery (LFQbench)</h2>")
    S.append(table(["level (DIA-NN 2.5)", "E. coli A/B", "Yeast A/B", "n"],
                   [["protein (raw Σ-precursor)", "+1.03", "−1.63", "4,855"],
                    ["peptide", "+1.05", "−1.65", "20,227"],
                    ["ion", "+1.05", "−1.66", "21,533"]]))
    S.append('<div class="key">Expected A/B: E. coli +1.0, yeast −1.58 (=−log₂3), human 0. '
             'All levels agree (E. coli ≈ +1.05, yeast ≈ −1.65) because the benchmark quantifies '
             'from <b>raw precursor intensity</b> (<code>Precursor.Quantity</code>; protein = sum '
             'over the group), not DIA-NN\'s normalized <code>PG.MaxLFQ</code>/'
             '<code>Precursor.Normalised</code>. With MaxLFQ the protein rollup would read +0.91 '
             '(compressed) — and for 2.5 the compression is severe. <b>Why we use raw: §8.</b> '
             '<i>DIA-NN 2.6, run through the identical grid, matches 2.5 on the raw axis '
             '(protein +1.03/−1.63, ion +1.05/−1.65) — see §8.</i></div>')
    S.append(img("stage1_diann25_ratios_protein.png", "DIA-NN 2.5, protein level (raw Σ-precursor). Plasma orange, yeast blue, E. coli green."))
    S.append(img("stage1_diann25_ratios_ion.png", "DIA-NN 2.5, ion level — agrees with protein on raw quantity."))
    S.append(img("stage1_diann18_ratios_protein.png", "DIA-NN 1.8, protein level (raw Σ-precursor) — +1.02/−1.61, matches 2.5 on raw quantity."))
    S.append(img("stage1_diann18_ratios_ion.png", "DIA-NN 1.8, ion level — +1.03/−1.63, matches 2.5 ion (+1.05/−1.66)."))

    # ---- Stage 1: quant ----
    S.append("<h2>3. Stage 1 — quant correlation vs truth (ground-truth only)</h2>")
    S.append("<p>Simulated truth intensity vs recovered quantity, per species. All three versions "
             "(1.8/2.5/2.6) reach Spearman ~0.92–0.96 — strong rank fidelity and essentially a tie "
             "(within-run abundance ranking, the first axis of §8). The paper cannot produce this "
             "panel (no truth).</p>")
    S.append(img("stage1_diann25_quant.png", "DIA-NN 2.5 precursor quant correlation. ρ 0.92–0.96 per species."))
    S.append(img("stage1_diann18_quant.png", "DIA-NN 1.8 precursor quant correlation. ρ 0.93–0.96 per species — comparable to 2.5."))

    # ---- Stage 1: FDR ----
    S.append("<h2>4. Stage 1 — empirical FDR / recall (ground-truth only)</h2>")
    S.append('<div class="key">Headline: at the same nominal 1% q-value, DIA-NN 1.8\'s <b>true</b> '
             'error rate is ~2% at precursor and <b>~7.5% at protein</b>, while 2.5 holds precursor '
             '&lt;0.5%. Recall is comparable (ion 67% vs 66%, protein 88.6% vs 90%) — so at low '
             'abundance 2.5 detects as many precursors yet is better-calibrated. This independently '
             'reproduces the known "1.8 FDR control is poor" observation. Protein FDR is group-credit; '
             '"any-member" / "strict" definitions reported separately.</div>')
    S.append(table(["level", "1.8 FDR", "2.5 FDR", "2.6 FDR", "recall 1.8 / 2.5 / 2.6"],
                   [["ion", "1.9%", "0.49%", "0.57%", "66% / 67% / 67%"],
                    ["peptide", "2.0%", "0.51%", "0.59%", "67% / 70% / 69%"],
                    ["protein (group-credit: any / strict)", "7.5% / 7.8%", "1.44% / 1.65%", "1.42% / 1.63%", "90% / 88.6% / 88%"]]))
    S.append('<div class="key">DIA-NN <b>2.6 ≈ 2.5</b> on FDR/recall too — well-calibrated '
             '(ion ~0.57%, protein group-credit ~1.4%) and far below 1.8\'s 1.9%/7.5%. So 2.6 inherits '
             '2.5\'s better FDR control along with its quant-normalization behaviour (§8).</div>')
    S.append(img("stage1_diann25_fdr_blanksub.png", "DIA-NN 2.5 FDR/TPR, blank-subtracted. Protein group-credit FDR over the 1% line (inference inflation)."))
    S.append(img("stage1_diann26_fdr_blanksub.png", "DIA-NN 2.6 FDR/TPR, blank-subtracted — matches 2.5 (ion ~0.57%, protein ~1.4%)."))
    S.append(img("stage1_diann18_fdr.png", "DIA-NN 1.8 — FDR far over 1% (ion ~1.9%, protein group-credit ~7.5%)."))
    S.append("<p><small>Blank subtraction: the noise-source blank yields 1,039 IDs (1,035 human, "
             "4 yeast); subtracting blank-explained IDs from the FP side barely moves FDR "
             "(ion 0.53→0.49%), i.e. <b>most SIM false positives are not explained by IDs seen in "
             "the blank</b> — not a proof that every remaining FP is genuine, but blank leakage is "
             "not the driver.</small></p>")

    # ---- Stage 1: overlap ----
    S.append("<h2>5. Stage 1 — cross-engine ID overlap</h2>")
    S.append(img("stage1_overlap_18_vs_25.png", "1.8 vs 2.5 ID overlap (q<0.01, A+B pooled). Precursor/peptide ~90% shared; "
             "protein asymmetric — 1.8 reports 647 unique proteins vs 2.5's 90. Given 1.8's higher protein FDR (~7.5%), the "
             "extra IDs are plausibly enriched for false positives, but confirming that needs per-protein truth labels (not done here)."))

    # ---- Stage 1: abundance ----
    S.append("<h2>6. Stage 1 — sensitivity &amp; ratio bias vs TRUE abundance</h2>")
    S.append("<p>Binning by the blueprint's true abundance (shared across engines) removes the "
             "reported-intensity-scale confound. In the low-abundance transition zone DIA-NN 2.5 "
             "detects as many or more precursors as 1.8 (top row) and is better-calibrated (bottom "
             "row) — so the reported-abundance view that made 2.5 <i>look</i> less sensitive is an "
             "artifact; at matched true abundance 2.5 is not less sensitive.</p>")
    S.append(img("stage1_trueabund_ion_18_vs_25.png", "True-abundance: detection sensitivity (top) + ratio bias (bottom), 1.8 vs 2.5. Human flat = control."))
    S.append(img("stage1_ratiobias_ion_18_vs_25.png", "Ratio bias vs each engine's own reported abundance — within-engine diagnostic; 1.8 fans wider at low abundance."))

    # ---- Stage 2 ----
    S.append("<h2>7. Stage 2 — spike-in over REAL plasma (superimpose)</h2>")
    S.append("""<p>The simulated yeast+E.coli is overlaid on a real PYE-plasma <code>.d</code>
    (the actual plasma used for the mix). Contamination check: searching the plasma vs
    human+yeast+E.coli finds only ~0.4% yeast/E.coli — at the level expected from FDR-level
    false IDs and with no coherent high-confidence YE population, so no meaningful carryover
    was detected (low detectable contamination, not a proof of zero). YE intensity is
    calibrated (report-proxy) to sit at the real-PYE level relative to human:
    <b>yeast 0.40×, E. coli 0.63×</b> the human median (A/B oracle ratios preserved). Human is
    real background with no truth (<code>truth_scope=background_unknown</code>): scored
    empirically only; FDR/sensitivity restricted to the simulated spike-ins.</p>""")
    S.append('<div class="key">A simulator bug surfaced here and was fixed (§7a). Post-fix the '
             'blueprint A/B truth is nominal again (yeast 3.01, E. coli 0.51), so everything below is '
             'scored against correct truth.</div>')

    S.append("<h3>7a. A from_findings ratio-distortion bug — found &amp; fixed</h3>")
    S.append("""<p>Initial Stage-2 ratios looked badly off (E. coli "over-separated", yeast
    compressed, large spread). Root cause: TimSim <code>load_findings</code> scaled events by each
    sample's <i>own</i> median intensity, so for an A/B experiment the differing medians silently
    rescaled every cross-sample ratio by median(A)/median(B) ≈ 0.7 — the <i>simulated truth itself</i>
    was off-nominal. Stage 1 was affected too (human truth B/A came out 0.74 instead of 1.0 — the
    smoking gun), but hidden because the ratio panels human-anchor, which cancels a uniform factor.
    Fixed via a shared reference median (rustims PR #407, merged); sim-QC confirms blueprint A/B back
    to nominal. All Stage-2 numbers below use the corrected re-sim.</p>""")

    S.append("<h3>7b. Ratio recovery (corrected truth, raw quantity) — 1.8 vs 2.5</h3>")
    S.append(table(["protein A/B (expected ecoli +1.0, yeast −1.58)", "E. coli", "yeast", "human anchor (080 / 081)"],
                   [["DIA-NN 1.8", "+1.01", "−1.61", "0.00 / 0.00 (stable)"],
                    ["DIA-NN 2.5", "+1.02", "−1.62", "0.00 / 0.00 (stable)"]]))
    S.append('<div class="key">On real plasma, with <b>raw precursor quantity</b> both engines recover '
             'the spike-in ratios accurately and reproducibly (E. coli ≈ +1.0, yeast ≈ −1.6 on both '
             'backgrounds; human anchor ≈ 0). An earlier version of this section reported 2.5 as '
             '<i>compressed and unstable</i> (yeast −1.11, anchor swinging −0.31→−1.51). That was '
             '<b>not a real 2.5 defect</b> — it was DIA-NN 2.5\'s cross-run normalization '
             '(<code>PG.MaxLFQ</code>/<code>Precursor.Normalised</code>) misbehaving on a spike-in with '
             'controlled loading; quantifying from raw <code>Precursor.Quantity</code> removes it entirely '
             'and 2.5 lands on top of 1.8. <b>This is the headline of §8.</b></div>')
    S.append(img("stage2fix_18_ratios_protein.png", "DIA-NN 1.8, Stage-2 corrected, protein (raw Σ-precursor) — yeast −1.61, E. coli +1.01."))
    S.append(img("stage2fix_25_ratios_protein.png", "DIA-NN 2.5, Stage-2 corrected, protein (raw Σ-precursor) — yeast −1.62, E. coli +1.02; matches 1.8 (the MaxLFQ compression is gone)."))
    S.append(img("stage2fix_18v25_ratiobias.png", "Ratio bias vs abundance on real plasma, 1.8 vs 2.5 (corrected truth, raw quantity)."))

    S.append("<h3>7c. FDR / recall on real plasma (simulated spike-in)</h3>")
    S.append(table(["", "ion FDR", "ion recall", "protein FDR", "protein recall"],
                   [["DIA-NN 2.5", "0.83%", "62%", "0.94%", "90%"],
                    ["DIA-NN 2.6", "0.81%", "63%", "0.94%", "90%"],
                    ["DIA-NN 1.8", "0.88%", "59%", "1.73%", "89%"]]))
    S.append('<div class="caveat"><b>Scope — read FDR like-for-like.</b> Stage 2 scores FDR only '
             'over the simulated <b>spike-in (YE)</b>: human is real plasma with no blueprint '
             '(<code>background_unknown</code>), so human false discoveries are invisible to the '
             'scorer. Stage 1\'s headline protein FDR (1.8 ~7.5%) is over <i>all</i> species and is '
             '<b>~79% human</b> false groups — so it is NOT comparable to Stage 2. Restricting Stage 1 '
             'to YE-only gives protein FDR <b>2.07% (1.8)</b>, essentially equal to Stage 2\'s 1.73%. '
             '<b>So the spike-in FDR is consistent across blank-SIM and real-plasma</b> — DIA-NN is not '
             '"better on real plasma"; we just cannot see the human-background errors there. The 7.5% '
             'all-species figure remains a valid, separate statement about 1.8\'s human protein-inference '
             'inflation.</div>')
    S.append(img("stage2fix_25_fdr.png", "DIA-NN 2.5 Stage-2 simulated-spike-in FDR (left) &amp; recall (right), corrected, human excluded."))
    S.append(img("stage2fix_26_fdr.png", "DIA-NN 2.6 Stage-2 simulated-spike-in FDR (left) &amp; recall (right) — matches 2.5 (ion 0.81%, protein 0.94%)."))
    S.append(img("stage2fix_18_fdr.png", "DIA-NN 1.8 Stage-2 simulated-spike-in FDR (left) &amp; recall (right) — protein FDR 1.73% vs 2.5's 0.94%; recall comparable."))

    S.append("<h3>7d. Bias &amp; sensitivity vs TRUE abundance — blank vs real plasma</h3>")
    S.append(img("stage1_vs_stage2fix_trueabund.png", "Stage 1 (blank) vs Stage 2 (real plasma, corrected truth) — ratio bias & detection sensitivity vs true abundance. Human panel Stage-2-empty (no truth)."))

    # ---- 1.8 vs 2.5 quant ----
    S.append("<h2>8. DIA-NN 1.8 vs 2.5 vs 2.6 quant — the cross-run normalization (2.6 did not fix it)</h2>")
    S.append("""<p>2.5/2.6 are the more sensitive identifiers (why we use them for the one-shot
    blank/plasma background ID), but their <i>quantitative ratios</i> looked far noisier than 1.8
    throughout the drafts above. Splitting "quant" into two axes and controlling the confounds an
    independent review raised (precursor-set coverage, human-anchoring) isolates the cause cleanly.
    <b>DIA-NN 2.6 (released 2026-06-10) was run through the identical grid</b> (same SIM <code>.d</code>s,
    same unified <code>diann-pye1.cfg</code> — only the binary differs) as a third point on the
    version axis.</p>""")
    S.append(table(["axis", "metric", "DIA-NN 1.8", "DIA-NN 2.5", "DIA-NN 2.6", "verdict"],
                   [["within-run abundance", "Spearman ρ (truth)", "0.93–0.94", "0.92–0.93", "0.92–0.93",
                     "≈ tie across all three"],
                    ["cross-sample ratio precision", "IQR log₂(A/B), default Normalised, Stage 2",
                     "0.19–0.35", "1.6", "1.3–1.4", "1.8 far tighter; 2.6 ≈ 2.5"]]))
    S.append('<div class="key"><b>The spread is the 2.5/2.6 cross-run normalization, and it is '
             'selectable.</b> Re-reading the <i>same</i> report under different quantity columns '
             '(no re-search): on Stage 2, 2.5\'s ion-level ratio IQR is 1.64 with '
             '<code>Precursor.Normalised</code> but <b>0.23 with raw <code>Precursor.Quantity</code></b> '
             '(and 0.22 with <code>Ms1.Area</code>) — landing next to 1.8\'s 0.18; <b>2.6 behaves the '
             'same</b> (1.39 → 0.24). This is not a coverage artifact: <code>Precursor.Normalised</code> '
             'is a strictly-positive rescale of <code>Precursor.Quantity</code> (identical precursor '
             'support), and on the common feature set quantified under <i>every</i> column the IQRs are '
             'unchanged. The same switch removes a <b>compression bias</b> (2.5/2.6 Stage-2 ratios pulled '
             '~0.1–0.34 log₂ toward 1:1). It is <i>not</i> MaxLFQ (a protein rollup — our endpoint is '
             'precursor-level) and <i>not</i> the QuantUMS peak extraction (raw is fine): both '
             '<code>.Normalised</code> columns are bad, both raw columns are fine. 1.8\'s normalization '
             'is benign (0.19→0.16); 2.5\'s and 2.6\'s are ~6–7× regressions. <b>2.6 did not fix this</b> '
             '(marginally milder on real plasma than 2.5, but the same behavior). The benchmark therefore '
             'quantifies from raw precursor intensity everywhere (§2/§7), which is why all three versions '
             'agree above.</div>')
    S.append(table(["quantity column (Stage 2)", "1.8 ecoli/yeast IQR", "2.5 ecoli/yeast IQR", "2.6 ecoli/yeast IQR"],
                   [["Precursor.Normalised (default)", "0.19 / 0.35", "1.64 / 1.61", "1.39 / 1.27"],
                    ["Precursor.Quantity (raw)", "0.16 / 0.21", "0.23 / 0.30", "0.24 / 0.29"],
                    ["Ms1.Area (raw)", "0.15 / 0.21", "0.22 / 0.31", "0.22 / 0.32"]]))
    S.append(img("quant_18v25_iqr.png", "A/B ratio precision (IQR) on raw Precursor.Quantity — the benchmark standard. 1.8 ≈ 2.5 ≈ 2.6 in both stages."))
    S.append(img("quant_18v25_normalization.png", "Mechanism: IQR by quantity column, Stage 2. The spread lives in the .Normalised columns for BOTH 2.5 and 2.6 (not 1.8); raw Quantity / Ms1.Area recover 1.8 levels."))
    S.append("""<p><small>Accuracy companion (residual vs blueprint truth, human-anchored): on raw
    quantity all three versions are accurate and precise (bias ≤ 0.05 log₂, median|error| ≈ 0.09–0.13);
    2.5 and 2.6 on the default Normalised column are both noisier <i>and</i> ratio-compressed on real
    plasma. Practical knobs: score from <code>Precursor.Quantity</code> (free, no re-search;
    <code>load_observations(quant_col=…)</code>) or re-run with <code>--no-norm</code>. <b>Version-axis
    takeaway: 2.6 ≈ 2.5 on quant — the normalization regression carries forward.</b> Full write-up +
    reproduction: <code>docs/finding-diann-1.8-vs-2.5-quant.md</code>,
    <code>scripts/diag_quant_{iqr,residual,method}.py</code>.</small></p>""")

    # ---- Real vs SIM validation ----
    S.append("<h2>9. Real vs SIM — validation against experimental data</h2>")
    S.append("""<p>Ute provided the <b>real experimental PYE1 runs</b> behind this simulated method
    (the G250506 G-site acquisitions — the <i>same batch</i> as the Stage-2 reference plasma): 12
    runs, 6 replicates × Sample A/B. We searched them with DIA-NN 1.8/2.5/2.6 under the <i>same
    unified <code>diann-pye1.cfg</code></i> used for SIM (only the binary differs) and read RAW
    Precursor.Quantity. Real PYE1 is itself a ratio test (nominal A/B yeast 3.0, E.&nbsp;coli 0.5,
    human 1.0), so the ratio-vs-abundance behaviour compares directly to SIM. No blueprint exists on
    real data, so FDR is not scorable — but the <b>ratio axis</b> (the one that matters here) is.</p>""")
    S.append('<div class="key"><b>Cross-check:</b> our unified 2.5/2.6 reproduce Ute\'s own report '
             'almost exactly (global A/B ecoli +0.89–0.90 / yeast −1.32–1.33 vs her +0.89/−1.33, '
             '~27–28k IDs/run) — so the comparison is sound, not a settings artifact.</div>')
    S.append(img("real_vs_sim_compression.png",
                 "A/B ratio vs abundance octile (1=lowest…8=highest), DIA-NN 1.8, raw quant, "
                 "human-anchored. REAL is the mean of the 6 single-A-vs-single-B injection pairs "
                 "(matched to SIM's 1-vs-1 structure). REAL stays flat across abundance; SIM "
                 "over-separates hard at low abundance; human (control) flat."))
    S.append(table(["lowest-abundance octile (A/B), matched 1-vs-1", "E. coli", "yeast"],
                   [["REAL G PYE1 (mean of 6 pairs)", "+0.98  (flat)", "−1.39  (flat)"],
                    ["SIM Stage 2 (real plasma)", "+1.44", "−2.29"],
                    ["SIM Stage 1 (blank)", "+1.83", "−2.40"],
                    ["nominal", "+1.00", "−1.58"]]))
    S.append('<div class="key"><b>The real-vs-SIM gap: TimSim produces a spurious abundance-dependent '
             'ratio <i>over-separation</i> at low input that real data does not.</b> On the matched '
             '1-vs-1 comparison (single A vs single B run, the same structure as SIM), <b>real ratios '
             'are flat across abundance</b> (E.&nbsp;coli ≈ +0.95, yeast ≈ −1.40 at every octile — a '
             'uniform mild global compression vs nominal, but no abundance dependence), while <b>SIM '
             'blows outward at low abundance</b> (E.&nbsp;coli +1.4–1.8, yeast −2.3–2.4 at the lowest '
             'octile), converging to nominal only at high abundance. Same pipeline, same 1-vs-1 '
             'structure, opposite low-end behaviour → the over-separation is a property of the '
             '<b>simulated signal</b>, not the analysis. Engine-independent (1.8/2.5/2.6 agree, per §8) '
             'and present on both blank (Stage 1) and real-plasma (Stage 2) backgrounds. Most plausibly '
             'the simulator under-models the <b>noise floor / ion statistics at low input</b> (low-count '
             'ions get extreme fold-changes instead of regressing toward nominal) — though this analysis '
             'establishes the <i>observed</i> over-separation, not the exact mechanism.</div>')
    S.append('<div class="caveat"><b>What this is NOT (corrected after review):</b> an earlier version '
             'reported real data <i>compressing</i> at low abundance. That was an artifact of averaging '
             'all 6 real replicates per condition (<code>min_replicates=1</code>, linear-mean stabilises '
             'low-abundance ratios); the <b>matched 1-vs-1</b> real curve above is flat. The robust claim '
             'is the <i>SIM over-separation</i>, not real compression. Further caveats: the x-axis is '
             'reported-abundance <i>rank</i> (octile), not a matched absolute input across datasets '
             '(different intensity scales) — so the claim is the within-dataset shape; complete-case '
             'selection near the detection limit can bias low-octile ratios (controlled here only by the '
             'same-pipeline comparison); human-anchoring is a constant shift (shape-invariant). '
             'Reproduce: <code>scripts/plot_real_vs_sim_compression.py</code>.</div>')

    S.append("<h3>9a. Cause of the over-separation — the low-intensity peak is under-quantified</h3>")
    S.append("""<p>The gap can be localised, because SIM has blueprint truth. Binning the <b>truth</b>
    ratio by true abundance, it is exactly nominal and perfectly flat (the seed is abundance-independent
    by construction) — so the over-separation is introduced <i>after</i> the seed, in TimSim's rendering
    + DIA-NN's quant. Decomposing it per A/B-pair member (quant recovery
    <code>log2(observed/truth)</code> vs true abundance):</p>""")
    S.append(table(["recovery log2(obs/truth), lowest→highest octile", "higher member (big peak)", "lower member (small peak)"],
                   [["E. coli", "+0.04 → −0.07  (flat)", "−0.36 → −0.06"],
                    ["yeast", "+0.11 → −0.07  (flat)", "−0.51 → −0.06"]]))
    S.append(img("sim_lowabund_underquant.png",
                 "Per-member quant recovery vs true abundance, SIM Stage 2. The big peak is accurate "
                 "everywhere; the small peak is under-quantified at low abundance; the shaded gap = the "
                 "over-separation."))
    S.append('<div class="key"><b>Mechanism, pinned:</b> the <i>higher</i>-intensity member of each '
             'pair is recovered accurately at all abundances, while the <i>lower</i> member is '
             'under-quantified at low abundance (~0.4 E.&nbsp;coli / ~0.6 yeast log2 in the bottom octile). '
             'The higher−lower gap <b>equals</b> the over-separation to the decimal — it fully accounts '
             'for it — and is member-symmetric (holds where A is the big side and where B is), so it is '
             'about small-vs-big peak, not A/B or species. <b>Code-level cause</b> (confirmed by '
             'reading the rustims rendering path): the synthetic renderer uses deterministic fractional '
             'intensities with a hard <code>&lt;1.0</code> per-pixel threshold + rounding (no count '
             'sampling), so the dimmer member loses proportionally more sub-unit contributions before '
             'they are written. Tellingly, the real-data noise/superimpose is injected <i>after</i> that '
             'threshold, so it cannot restore the erased signal (it only adds interfering peaks → mild '
             'compression) — which is why adding real-data noise did <b>not</b> remove the over-separation. '
             '<b>Fix</b> (right place = the Rust renderer before the <code>1.0</code> filter): a per-pixel '
             'count model (e.g. Poisson + detector/background counts) for DIA MS2 — write-up + citations '
             'in <code>docs/timsim-lowinput-noise-floor.md</code> (ready to file as a rustims issue). '
             'Reproduce: <code>scripts/plot_lowabund_underquant.py</code>.</div>')

    S.append("<h3>9b. Cross-engine ID overlap on real data — 1.8 vs 2.x</h3>")
    S.append("<p>Pooled over all 12 real runs (q&lt;0.01). 2.5/2.6 identify ~23–25% more precursors "
             "than 1.8 (35k vs 28.5k); 2.6 ≈ 2.5.</p>")
    S.append(table(["overlap (precursor)", "shared", "1.8-only", "2.x-only", "Jaccard"],
                   [["1.8 vs 2.5", "26,227 (92% of 1.8)", "2,344", "8,904", "0.70"],
                    ["1.8 vs 2.6", "26,319 (92% of 1.8)", "2,252", "9,499", "0.69"],
                    ["2.5 vs 2.6", "33,345", "1,786", "2,473", "0.89"]]))
    S.append(img("real_overlap_18_vs_25.png", "Real PYE1 ID overlap, DIA-NN 1.8 vs 2.5 (precursor/peptide/protein)."))
    S.append(img("real_overlap_18_vs_26.png", "Real PYE1 ID overlap, DIA-NN 1.8 vs 2.6 — same picture as 2.5."))
    S.append('<div class="key"><b>2.x ≈ 1.8 + a large extra precursor population</b> (92% of 1.8 is '
             'recovered by 2.x, which adds ~9k more; 1.8 keeps ~2.3k unique, so not a strict superset). '
             '2.5 and 2.6 are nearly the same ID set (Jaccard 0.89). <b>At protein level the gap '
             'collapses</b> — 1.8 vs 2.5 share 4,289 with balanced uniques (1.8-only 628 / 2.5-only 672). '
             '<b>Real-vs-SIM structural difference:</b> on SIM (§5) 1.8 had 647 unique proteins vs 2.5\'s '
             '90 — the blueprint exposed 1.8\'s protein over-reporting as largely false; on real the '
             'protein-unique is symmetric (628 vs 672). So the <i>overlap structure</i> itself differs '
             'real-vs-SIM, complementing the §9 low-abundance gap.</div>')

    S.append("<h3>9c. Real-data ratio panels — DIA-NN 1.8 vs 2.5</h3>")
    S.append("<p>LFQbench-style A/B ratio panels on the real G PYE1 runs (ion level, raw "
             "Precursor.Quantity, human-anchored). Both engines show the same modest <b>global</b> "
             "compression vs nominal — E. coli +0.94 (1.8) / +0.89 (2.5) vs +1.0; yeast −1.37 / −1.33 "
             "vs −1.58 — with the human control at 0. 2.5 quantifies more precursors (e.g. ecoli "
             "11,323 vs 1.8's 9,923). This is the real-data counterpart of the §2 SIM ratio panels; "
             "the global compression is a real DIA-quant property (present in both engines), distinct "
             "from the abundance-<i>dependent</i> SIM over-separation of §9.</p>")
    S.append(img("real_ratios_diann18.png", "DIA-NN 1.8, real G PYE1 — ion-level A/B violins + MA scatter (raw quant)."))
    S.append(img("real_ratios_diann25.png", "DIA-NN 2.5, real G PYE1 — same view; slightly more compressed, more precursors."))

    # ---- caveats ----
    S.append("<h2>10. Scope &amp; caveats</h2>")
    S.append('<div class="caveat"><b>Superposition limit:</b> Stage 2 adds simulated signal '
             '<i>after</i> the real acquisition — it does NOT reproduce ion suppression, source/'
             'fragmentation competition, fill-time, or detector saturation. It tests the '
             '<b>dense-background + DIA-NN-processing</b> contribution to the low-abundance effect, '
             'not the full real-plasma phenomenon.</div>')
    S.append("""<ul>
    <li><b>Calibration is report-proxy</b> (DIA-NN report quant, not raw MS1/MS2 peak areas) — the
    raw PYE-mix <code>.d</code> would upgrade it to gold-standard. E. coli landed at 0.63× vs 0.53×
    target, within the proxy's resolution.</li>
    <li><b>Same-background A/B pairing</b> gives a perfect human anchor but likely <b>amplifies</b>
    the ratio distortion vs real separate injections (identical, perfectly-correlated interference).</li>
    <li>The earlier "lowest-decile collapse to 1:1" was background-dependent (080 only) and has been
    dropped as an over-read.</li>
    <li>Cross-engine integrity: seed from one engine, score against the union blueprint, read
    differences — not absolute self-recall.</li>
    </ul>""")
    S.append("<p><small>PlasmaBENCH · panels: P1 ratios, P3 sensitivity/bias, P5 FDR, P6 quant, "
             "ID-overlap, 1.8/2.5/2.6 quant (§8), real-vs-SIM validation (§9) · SUT: DIA-NN "
             "1.8/2.5/2.6 · quantity basis: raw Precursor.Quantity · 90-check test suite.</small></p>")
    S.append("</body></html>")
    return "\n".join(S)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/report.html", type=Path)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(build(), encoding="utf-8")
    kb = args.out.stat().st_size / 1024
    print(f"wrote {args.out}  ({kb:.0f} KB, self-contained)")


if __name__ == "__main__":
    main()
