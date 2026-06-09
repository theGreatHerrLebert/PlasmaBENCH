"""Abundance-resolved ratio BIAS for two engines/versions, per species.

Bins each species by abundance and plots the median log2(A/B) per decile against
the expected ratio, overlaying both reports — so abundance-dependent bias (e.g.
low-abundance compression) is directly visible and comparable.

Usage:
    python scripts/plot_ratio_bias.py \
        --report-a results/diann-1.8-AB-v3/report.tsv --label-a "1.8" \
        --report-b results/diann-2.5-AB-v3/report.parquet --label-b "2.5" \
        --sim-dir simulations/stage1/dia --level ion \
        --out results/figures/stage1_ratiobias_ion_18_vs_25.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import build_manifest, load_observations
from plasmabench.plots import ratios as R


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-a", required=True, type=Path)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--report-b", required=True, type=Path)
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--level", default="ion",
                    choices=["protein", "peptide", "modified_peptide", "ion", "precursor"])
    ap.add_argument("--convention", choices=["A_over_B", "B_over_A"], default="A_over_B")
    ap.add_argument("--aggregation", default="linear_mean",
                    choices=["linear_mean", "geometric_mean", "median"])
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--y-halfspan", type=float, default=1.0)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    conv = R.A_OVER_B if args.convention == "A_over_B" else R.B_OVER_A
    wides = {}
    for label, report in ((args.label_a, args.report_a), (args.label_b, args.report_b)):
        manifest = build_manifest(report, args.sim_dir, f"s1-{label}")
        obs = load_observations(report, manifest)
        wide, _ = R.build_ratio_table(obs, level=args.level, aggregation=args.aggregation)
        wides[label] = wide

    table = R.plot_ratio_bias(wides, args.out, conv=conv, level=args.level,
                              n_bins=args.n_bins, y_halfspan=args.y_halfspan)
    # compact print: low vs high decile median per species/label
    for label in wides:
        t = table[table["label"] == label]
        print(f"\n{label} ({args.level}):")
        for sp in R.SPECIES_ORDER:
            d = t[t["species"] == sp].sort_values("x")
            if len(d):
                print(f"  {sp:6s} low-decile {d.iloc[0]['median']:+.2f}  "
                      f"high-decile {d.iloc[-1]['median']:+.2f}  "
                      f"(expected {conv.to_display(R.EXPECTED_LOG2_BA[sp]):+.2f})")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
