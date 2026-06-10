"""Diagnostic: is DIA-NN 2.5's wide A/B ratio spread driven by its NORMALIZATION?

Re-reads the SAME report under different quantity columns (no re-running DIA-NN) and
reports per-species ratio IQR + residual-vs-truth. If 2.5's spread collapses when using
the raw Precursor.Quantity instead of Precursor.Normalised, the cross-run normalization is
the culprit and is selectable. Ms1.Area / Ms1.Normalised probe the MS1 vs MS2 quant path.

Usage:
    python scripts/diag_quant_method.py --report results/diann-2.5-s2fix-r080/report.parquet \
        --sim-dir simulations/stage2/dia-fix --experiment stage2 --human-scope background_unknown
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

QUANT_COLS = ["Precursor.Normalised", "Precursor.Quantity", "Ms1.Normalised", "Ms1.Area"]


def iqr(x):
    x = pd.Series(x).dropna()
    return float(np.subtract(*np.percentile(x, [75, 25]))) if len(x) >= 4 else float("nan")


def truth_ba(truth):
    t = truth[truth["truth_intensity"] > 0].copy()
    t["feature"] = t["sequence_modified"] + t["charge"].astype(str)
    per = t.groupby(["feature", "sample"], as_index=False)["truth_intensity"].mean()
    w = per.pivot(index="feature", columns="sample", values="truth_intensity").dropna(subset=["A", "B"])
    w = w[(w["A"] > 0) & (w["B"] > 0)]
    out = pd.DataFrame({"truth_log2_ba": np.log2(w["B"] / w["A"])}, index=w.index)
    sp = t.groupby("feature")["species"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna()
    out["species"] = out.index.map(sp)
    out = out.dropna(subset=["species"])
    h = out.loc[out["species"] == "HUMAN", "truth_log2_ba"]
    out["truth_log2_ba"] -= float(h.median()) if len(h) else 0.0
    return out["truth_log2_ba"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="diag")
    ap.add_argument("--human-scope", default="background_unknown",
                    choices=["simulated", "background_unknown"])
    ap.add_argument("--q-value-max", type=float, default=0.01)
    args = ap.parse_args()

    manifest = build_manifest(args.report, args.sim_dir, args.experiment, args.human_scope)
    truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)
    tba = truth_ba(truth)

    avail = pd.read_parquet(args.report).columns if args.report.suffix.lower() in (".parquet", ".pq") \
        else pd.read_csv(args.report, sep="\t", nrows=1).columns
    cols = [c for c in QUANT_COLS if c in set(avail)]

    # Build every column's ratio table first, so we can also report on the COMMON
    # feature set (precursors quantified under EVERY column) — otherwise a lower IQR
    # could reflect a smaller/stronger subset rather than the quantity basis. NB
    # Precursor.Normalised is a strictly-positive per-precursor rescale of
    # Precursor.Quantity, so those two share identical support; only Ms1.* differ.
    tables = {}
    for qc in cols:
        obs = load_observations(args.report, manifest, quant_col=qc)
        wide, _ = build_ratio_table(obs, level="ion", q_value_max=args.q_value_max, human_anchor=True)
        tables[qc] = wide.set_index("feature")
    common = set.intersection(*[set(t.index) for t in tables.values()]) if tables else set()

    print(f"report: {args.report}   |common feature set| = {len(common):,}")
    print(f"{'quant_col':22} {'species':7} {'n':>6} {'IQR':>7} {'bias':>7} {'med|e|':>7}"
          f"  {'IQR_common':>10} {'n_common':>8}")
    for qc in cols:
        wide = tables[qc]
        wc = wide[wide.index.isin(common)]
        for sp in [s for s in SPECIES_ORDER if s != "HUMAN"]:
            d = wide[wide["species"] == sp]
            dc = wc[wc["species"] == sp]
            e = (d["log2_ba"] - d.index.map(tba)).dropna()
            print(f"{qc:22} {sp:7} {len(d):>6} {iqr(d['log2_ba']):>7.3f} "
                  f"{(np.median(e) if len(e) else float('nan')):>+7.3f} "
                  f"{(np.median(np.abs(e)) if len(e) else float('nan')):>7.3f}"
                  f"  {iqr(dc['log2_ba']):>10.3f} {len(dc):>8}")


if __name__ == "__main__":
    main()
