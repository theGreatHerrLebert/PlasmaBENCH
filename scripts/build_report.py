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
                   [["protein (PG.MaxLFQ)", "+0.91", "−1.58", "4,842"],
                    ["peptide", "+1.07", "−1.65", "20,227"],
                    ["ion", "+1.08", "−1.65", "21,533"]]))
    S.append('<div class="key">Expected A/B: E. coli +1.0, yeast −1.58 (=−log₂3), human 0. '
             'The protein-level MaxLFQ rollup <b>compresses</b> ratios vs precursor level — '
             'E. coli +0.91 (protein) vs +1.08 (ion).</div>')
    S.append(img("stage1_diann25_ratios_protein.png", "DIA-NN 2.5, protein level. Plasma orange, yeast blue, E. coli green."))
    S.append(img("stage1_diann25_ratios_ion.png", "DIA-NN 2.5, ion level — denser, less compressed than protein."))
    S.append(img("stage1_diann18_ratios_protein.png", "DIA-NN 1.8, protein level — less compressed than 2.5 (E. coli +1.04); human anchor differs."))

    # ---- Stage 1: quant ----
    S.append("<h2>3. Stage 1 — quant correlation vs truth (ground-truth only)</h2>")
    S.append("<p>Simulated truth intensity vs recovered quantity, per species. Spearman 0.92–0.96 "
             "(DIA-NN 2.5) — strong rank fidelity. The paper cannot produce this panel (no truth).</p>")
    S.append(img("stage1_diann25_quant.png", "DIA-NN 2.5 precursor quant correlation. ρ 0.92–0.96 per species."))

    # ---- Stage 1: FDR ----
    S.append("<h2>4. Stage 1 — empirical FDR / recall (ground-truth only)</h2>")
    S.append('<div class="key">Headline: at the same nominal 1% q-value, DIA-NN 1.8\'s <b>true</b> '
             'error rate is ~2% at precursor and <b>~7.5% at protein</b>, while 2.5 holds precursor '
             '&lt;0.5%. Recall is comparable (ion 67% vs 66%, protein 88.6% vs 90%) — so at low '
             'abundance 2.5 detects as many precursors yet is better-calibrated. This independently '
             'reproduces the known "1.8 FDR control is poor" observation. Protein FDR is group-credit; '
             '"any-member" / "strict" definitions reported separately.</div>')
    S.append(table(["level", "1.8 FDR", "2.5 FDR", "1.8 / 2.5 recall"],
                   [["ion", "1.9%", "0.49%", "66% / 67%"],
                    ["peptide", "2.0%", "0.51%", "67% / 70%"],
                    ["protein (group-credit: any / strict)", "7.5% / 7.8%", "1.44% / 1.65%", "90% / 88.6%"]]))
    S.append(img("stage1_diann25_fdr_blanksub.png", "DIA-NN 2.5 FDR/TPR, blank-subtracted. Protein group-credit FDR over the 1% line (inference inflation)."))
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

    S.append("<h3>7b. Ratio recovery (corrected truth) — 1.8 vs 2.5</h3>")
    S.append(table(["protein A/B (expected ecoli +1.0, yeast −1.58)", "E. coli", "yeast", "human anchor (080 / 081)"],
                   [["DIA-NN 1.8", "+1.01", "−1.64", "−0.64 / −0.69 (stable)"],
                    ["DIA-NN 2.5", "+1.12", "−1.11", "−0.31 / −1.51 (unstable)"]]))
    S.append('<div class="key">On real plasma, <b>DIA-NN 1.8 recovers the spike-in ratios accurately '
             'and reproducibly</b> (yeast −1.64, E. coli +1.01 on both backgrounds, stable human anchor). '
             '<b>2.5 is unstable</b>: its quant of the real human background (the anchor) swings from '
             '−0.31 to −1.51 between two near-identical plasma backgrounds, so its YE ratios jump around '
             '— 2.5\'s QuantUMS/MaxLFQ normalization looks destabilized by the dense plasma + superimposed '
             'spike-in. <b>Preliminary</b>: may partly be the same-background pairing interacting with '
             '2.5\'s cross-run normalization; under investigation.</div>')
    S.append(img("stage2fix_18_ratios_protein.png", "DIA-NN 1.8, Stage-2 corrected, protein — accurate (yeast −1.64, E. coli +1.01)."))
    S.append(img("stage2fix_25_ratios_protein.png", "DIA-NN 2.5, Stage-2 corrected, protein — yeast compressed here (−1.11) and unstable across backgrounds (see anchor)."))
    S.append(img("stage2fix_18v25_ratiobias.png", "Ratio bias vs abundance on real plasma, 1.8 vs 2.5 (corrected truth)."))

    S.append("<h3>7c. FDR / recall on real plasma (simulated spike-in)</h3>")
    S.append(table(["", "ion FDR", "ion recall", "protein FDR", "protein recall"],
                   [["DIA-NN 2.5", "0.83%", "62%", "0.94%", "90%"],
                    ["DIA-NN 1.8", "0.88%", "59%", "1.73%", "89%"]]))
    S.append('<div class="key">On real plasma both engines sit ~0.8–0.9% at ion level; 1.8\'s protein '
             'FDR (1.73%) is still above 2.5\'s (0.94%) but the gap is far smaller than Stage 1 '
             '(7.5% vs 1.5%) — the dense background narrows the version difference.</div>')
    S.append(img("stage2fix_25_fdr.png", "DIA-NN 2.5 Stage-2 simulated-spike-in FDR/recall (corrected, human excluded)."))

    S.append("<h3>7d. Bias &amp; sensitivity vs TRUE abundance — blank vs real plasma</h3>")
    S.append(img("stage1_vs_stage2fix_trueabund.png", "Stage 1 (blank) vs Stage 2 (real plasma, corrected truth) — ratio bias & detection sensitivity vs true abundance. Human panel Stage-2-empty (no truth)."))

    # ---- caveats ----
    S.append("<h2>8. Scope &amp; caveats</h2>")
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
    S.append("<p><small>PlasmaBENCH · branch <code>plots-analysis-panels</code> · panels: "
             "P1 ratios, P3 sensitivity/bias, P5 FDR, P6 quant, ID-overlap · 85-check test suite.</small></p>")
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
