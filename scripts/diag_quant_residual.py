"""Diagnostic: per-precursor A/B ratio residual vs blueprint truth (accuracy + precision).

Companion to diag_quant_iqr.py. IQR of observed log2(A/B) measures spread but is blind
to bias (centered-but-wide vs shifted-but-tight). With exact blueprint truth we score the
residual e = observed_log2(A/B) − truth_log2(A/B) per precursor and report, per species:
  bias    = median(e)        signed accuracy
  precis. = IQR(e), MAD(e)   spread (precision)
  combined= median|e|        accuracy+precision in one number

Both sides are HUMAN-ANCHORED (observed and truth each shifted by their own human-median
log2 B/A). This is intentional: the raw observed-vs-truth intensity scales differ by an
arbitrary global convention, so unanchored `bias` just measures that convention, not quant
error (1.8 and 2.5 show a uniform per-species offset = the global scale). Anchoring removes
it and leaves species-specific bias. IQR/MAD are anchor-invariant (a constant shift), so the
precision numbers are unaffected by this choice. Stage 2 has no human blueprint → truth side
anchors to 0 (real-plasma human is genuinely 1:1). Truth ratio per ion = log2 of blueprint
truth_intensity B vs A (events × relative_abundance, the actual simulated A/B).

Usage:
    python scripts/diag_quant_residual.py \
        --engine 1.8 results/diann-1.8-AB-v3/report.tsv \
        --engine 2.5 results/diann-2.5-AB-v3/report.parquet \
        --sim-dir simulations/stage1/dia --experiment stage1 --human-scope simulated
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots.ratios import build_ratio_table, SPECIES_ORDER


def _stats(e: pd.Series) -> dict:
    e = e.dropna()
    n = len(e)
    if n < 4:
        return dict(n=n, bias=float("nan"), iqr=float("nan"),
                    mad=float("nan"), med_abs=float("nan"))
    return dict(
        n=n,
        bias=float(np.median(e)),
        iqr=float(np.subtract(*np.percentile(e, [75, 25]))),
        mad=float(np.median(np.abs(e - np.median(e)))),
        med_abs=float(np.median(np.abs(e))),
    )


def truth_ratio_by_ion(truth: pd.DataFrame, human_anchor: bool = True) -> pd.DataFrame:
    """log2(truth_B / truth_A) per ion (sequence_modified+charge), transmitted-aware truth.

    Human-anchored to match the observed side: subtracts the human-median truth ratio so the
    species-specific bias is not confounded by the global intensity-scale convention. Stage 2
    has no human truth (real plasma, background_unknown) → anchor 0 (human truly 1:1)."""
    t = truth.copy()
    t = t[t["truth_intensity"] > 0]
    t["feature"] = t["sequence_modified"] + t["charge"].astype(str)
    # Blueprint truth is per-sample (one .d per sample in these sims → one run/sample), so
    # the mean over runs is a no-op here; if a sample ever had >1 replicate run with
    # differing truth, match the observed replicate subset before comparing (not needed now).
    per = t.groupby(["feature", "sample"], as_index=False)["truth_intensity"].mean()
    wide = per.pivot(index="feature", columns="sample", values="truth_intensity")
    wide = wide.dropna(subset=["A", "B"])
    wide = wide[(wide["A"] > 0) & (wide["B"] > 0)]
    out = pd.DataFrame({"truth_log2_ba": np.log2(wide["B"] / wide["A"])}, index=wide.index)
    # species per ion (single-species only; ambiguous dropped — matches ratio table guard)
    sp = t.groupby("feature")["species"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna()
    out["species"] = out.index.map(sp)
    out = out.dropna(subset=["species"])
    if human_anchor:
        h = out.loc[out["species"] == "HUMAN", "truth_log2_ba"]
        anchor = float(h.median()) if len(h) else 0.0
        out["truth_log2_ba"] = out["truth_log2_ba"] - anchor
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", nargs=2, action="append", metavar=("LABEL", "REPORT"),
                    required=True)
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="diag")
    ap.add_argument("--human-scope", default="simulated",
                    choices=["simulated", "background_unknown"])
    ap.add_argument("--q-value-max", type=float, default=0.01)
    args = ap.parse_args()

    truth_by_engine = {}
    obs_ratio_by_engine = {}
    for label, report in args.engine:
        report = Path(report)
        manifest = build_manifest(report, args.sim_dir, args.experiment, args.human_scope)
        truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)
        obs = load_observations(report, manifest)
        wide, _ = build_ratio_table(obs, level="ion", q_value_max=args.q_value_max,
                                    human_anchor=True)
        truth_by_engine[label] = truth_ratio_by_ion(truth)
        obs_ratio_by_engine[label] = wide.set_index("feature")[["log2_ba"]]
        print(f"[{label}] {report.name}  obs_ions={len(wide):,}  "
              f"truth_ions={len(truth_by_engine[label]):,}")

    print(f"\n{'species':8} {'engine':6} {'n':>7} {'bias':>8} {'IQR':>8} "
          f"{'MAD':>8} {'med|e|':>8}")
    for sp in [s for s in SPECIES_ORDER if s != "HUMAN"]:
        for label, _ in args.engine:
            tr = truth_by_engine[label]
            ob = obs_ratio_by_engine[label]
            j = ob.join(tr, how="inner")
            j = j[j["species"] == sp]
            e = j["log2_ba"] - j["truth_log2_ba"]
            s = _stats(e)
            print(f"{sp:8} {label:6} {s['n']:>7} {s['bias']:>+8.3f} {s['iqr']:>8.3f} "
                  f"{s['mad']:>8.3f} {s['med_abs']:>8.3f}")


if __name__ == "__main__":
    main()
