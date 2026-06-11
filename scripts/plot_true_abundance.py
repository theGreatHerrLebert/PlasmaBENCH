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
    ap.add_argument("--sim-dir-b", type=Path, default=None,
                    help="separate sim tree for report-b (cross-stage compare); default = --sim-dir")
    ap.add_argument("--level", default="ion",
                    choices=["ion", "precursor", "peptide", "modified_peptide"])
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    sim_b = args.sim_dir_b or args.sim_dir

    def truth_for(report, sim_dir):
        m = build_manifest(report, sim_dir, "s1")
        return m, pd.concat([load_truth(r, m) for r in m.runs], ignore_index=True)

    tables = {}
    for label, report, sim_dir in ((args.label_a, args.report_a, args.sim_dir),
                                    (args.label_b, args.report_b, sim_b)):
        m, truth = truth_for(report, sim_dir)
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
