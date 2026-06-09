"""P1 — three-species A/B ratio panels for a Stage-1 run (LFQbench-style).

Thin CLI over :mod:`plasmabench.plots.ratios`. Builds the run manifest from the SIM
tree (same as the quant CLI), loads the ``observations`` table, and draws the
per-species ratio violin + LFQbench MA scatter.

Defaults reproduce the PYE paper's axis (``log2(A/B)``: human 0, yeast −log2 3,
E. coli +1) at protein level (PG.MaxLFQ). Switch with ``--level``/``--convention``;
``--level precursor --convention B_over_A`` recovers the earlier precursor view.

Usage:
    python scripts/plot_stage1_ratios.py \
        --report results/diann-2.5-AB-v3/report.parquet \
        --sim-dir simulations/stage1/dia \
        --experiment stage1-diann-2.5 \
        --out results/figures/stage1_diann25_ratios.png
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
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sim-dir", required=True, type=Path)
    ap.add_argument("--experiment", default="stage1")
    ap.add_argument("--level", default="protein",
                    choices=["protein", "peptide", "modified_peptide", "ion", "precursor"])
    ap.add_argument("--convention", choices=["A_over_B", "B_over_A"], default="A_over_B")
    ap.add_argument("--aggregation", choices=["linear_mean", "geometric_mean", "median"],
                    default="linear_mean")
    ap.add_argument("--min-replicates", type=int, default=1)
    ap.add_argument("--q-value-max", type=float, default=0.01)
    ap.add_argument("--no-human-anchor", action="store_true")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    conv = R.A_OVER_B if args.convention == "A_over_B" else R.B_OVER_A
    manifest = build_manifest(args.report, args.sim_dir, args.experiment)
    obs = load_observations(args.report, manifest)

    wide, anchor = R.build_ratio_table(
        obs, level=args.level, aggregation=args.aggregation,
        min_replicates=args.min_replicates, q_value_max=args.q_value_max,
        human_anchor=not args.no_human_anchor,
    )
    summary = R.ratio_summary(wide, conv)
    print(f"{len(wide):,} {args.level}s quantified in both samples "
          f"(agg={args.aggregation}, anchor={anchor:+.3f} log2, convention={conv.name})")
    print(summary.to_string(index=False))

    R.plot_ratio_panels(wide, args.out, conv=conv, level=args.level, anchor=anchor)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
