"""Cross-engine ID overlap Venn (precursor / peptide / protein) for two reports.

Usage:
    python scripts/plot_id_overlap.py \
        --report-a results/diann-1.8-AB-v3/report.tsv --label-a "DIA-NN 1.8" \
        --report-b results/diann-2.5-AB-v3/report.parquet --label-b "DIA-NN 2.5" \
        --out results/figures/stage1_overlap_18_vs_25.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plasmabench.plots.overlap import id_sets, plot_overlap_venn


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-a", required=True, type=Path)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--report-b", required=True, type=Path)
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--q-value-max", type=float, default=0.01)
    ap.add_argument("--title", default="DIA-NN ID overlap (q<0.01, A+B pooled)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    sets_a = id_sets(args.report_a, args.q_value_max)
    sets_b = id_sets(args.report_b, args.q_value_max)
    table = plot_overlap_venn(sets_a, sets_b, args.label_a, args.label_b, args.out,
                              title=args.title)
    print(table.to_string(index=False))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
