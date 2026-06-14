"""LFQbench-style A/B ratio panels (violin + MA scatter) on REAL G PYE1 data, per engine.

Real data has no blueprint, so this uses the manifest-free real loader. RAW Precursor.Quantity
(benchmark standard); ion level (the precursor ratio distribution). One figure per engine.

Usage:
    python scripts/plot_real_ratios.py --report results/diann-1.8-realG/report.tsv \
        --label "DIA-NN 1.8" --out results/figures/real_ratios_diann18.png
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.loader import load_observations_real
from plasmabench.plots import ratios as R


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--label", default="DIA-NN")
    ap.add_argument("--level", default="ion")
    ap.add_argument("--quant-col", default="Precursor.Quantity")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    obs = load_observations_real(args.report, quant_col=args.quant_col)
    wide, anchor = R.build_ratio_table(obs, level=args.level, q_value_max=0.01, human_anchor=True)
    summ = R.ratio_summary(wide, R.A_OVER_B)
    print(f"{args.label}: {len(wide):,} {args.level}s, anchor {anchor:+.3f} ({args.quant_col})")
    print(summ.to_string(index=False))
    R.plot_ratio_panels(wide, args.out, conv=R.A_OVER_B, level=args.level, anchor=anchor)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
