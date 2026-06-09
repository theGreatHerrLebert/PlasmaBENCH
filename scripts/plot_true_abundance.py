"""Sensitivity + ratio bias vs TRUE (blueprint) abundance, comparing two engines.

Bins by ground-truth intensity (shared across versions → no reported-scale confound),
so the detection-sensitivity and ratio-bias curves are directly comparable.

Usage:
    python scripts/plot_true_abundance.py \
        --report-a results/diann-1.8-AB-v3/report.tsv --label-a "1.8" \
        --report-b results/diann-2.5-AB-v3/report.parquet --label-b "2.5" \
        --sim-dir simulations/stage1/dia --level ion \
        --out results/figures/stage1_trueabund_ion_18_vs_25.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_truth, load_observations
from plasmabench.plots.abundance import build_true_abundance, plot_true_abundance


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-a", required=True, type=Path)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--report-b", required=True, type=Path)
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--level", default="ion",
                    choices=["ion", "precursor", "peptide", "modified_peptide"])
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    # truth is the SAME blueprint for both versions → load once (manifest from report A).
    manifest = build_manifest(args.report_a, args.sim_dir, "s1")
    truth = pd.concat([load_truth(r, manifest) for r in manifest.runs], ignore_index=True)

    tables = {}
    for label, report in ((args.label_a, args.report_a), (args.label_b, args.report_b)):
        m = build_manifest(report, args.sim_dir, "s1")
        obs = load_observations(report, m)
        tables[label] = build_true_abundance(truth, obs, level=args.level, n_bins=args.n_bins)

    table = plot_true_abundance(tables, args.out, level=args.level)
    for label in tables:
        t = table[table["label"] == label]
        print(f"\n{label} ({args.level}) — sensitivity & bias by true-abundance decile:")
        for sp in ("ECOLI", "YEAST"):
            d = t[t["species"] == sp].sort_values("x")
            if len(d):
                lo, hi = d.iloc[0], d.iloc[-1]
                print(f"  {sp:6s} low: sens {100*lo['sensitivity']:.0f}% bias {lo['median_ratio']:+.2f} | "
                      f"high: sens {100*hi['sensitivity']:.0f}% bias {hi['median_ratio']:+.2f}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
